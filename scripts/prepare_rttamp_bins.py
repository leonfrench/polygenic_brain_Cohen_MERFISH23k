"""RT&T-AMP-MERFISH hemisections as square bins of cells, each the pseudobulk of the cells inside it, as for the Stereo-seq and MERFISH apps.

The prepared cell dataset (--cells, from prepare_rttamp_cells.py) gives the cells, their CCF positions (µm), labels, expression
matrix rows and the gene universe. It stores ranks, not counts, so the counts of the raw h5ad in --raw are read again, with source
genes merged into the universe as its gene_mapping.parquet records. Cells are assigned by centroid to the squares of one lattice
shared by all sections, so every bin has the same size and smaller sections simply have fewer bins. --bin-um sets the side; --bins
sets it so that the widest section spans that many squares. A bin's counts are the summed transcripts of its cells, ranked within
the bin as the cells are.

Checks: every cell's counts equal n_counts in the cell dataset, each section's bin gene totals equal its cell class pseudobulk,
and --check-bins bins per section match scipy.stats.rankdata and a pairwise AUROC.

  .conda-env/bin/python scripts/prepare_rttamp_bins.py --bin-um 178
"""
from pathlib import Path
import argparse, datetime, json, math, os, platform, shutil, time
import h5py
import numpy as np
import pandas as pd
import scipy

from prepare_rttamp_cells import EXPRESSION, LEVELS, RAW, section_counts, sha256, rank_rows, validate_rows

ROOT = Path(__file__).resolve().parents[1]
PARENT = {'cell_type': 'cell_group', 'cell_group': 'cell_class'}


def main_labels(codes, labels, local, n_bins):
    """Each bin's most common label (ties: first in sorted order) and that label's fraction of the bin's cells."""
    tally = np.bincount(local*len(labels) + codes, minlength=n_bins*len(labels)).reshape(n_bins, len(labels))
    top = tally.argmax(axis=1)
    return labels[top], tally[np.arange(n_bins), top]/tally.sum(axis=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--raw', type=Path, default=RAW, help='Folder of cellxgene_WT_labels.h5ad.')
    parser.add_argument('--cells', type=Path, default=ROOT/'data/rttamp-wt', help='Prepared cell dataset (prepare_rttamp_cells.py).')
    parser.add_argument('--sections', nargs='*', help='Section (sample) numbers of the cell dataset, e.g. 1 3. Default: all of them.')
    size = parser.add_mutually_exclusive_group(required=True)
    size.add_argument('--bin-um', type=float, help='Side of the squares in µm, rounded to whole coordinate units.')
    size.add_argument('--bins', type=int, help='Squares across the longer side of the widest section; sets the side for all sections.')
    parser.add_argument('--check-bins', type=int, default=64, help='Bins per section checked against scipy.stats.rankdata.')
    parser.add_argument('--out', type=Path, help='Default: the cell dataset folder name with -bins<µm>um.')
    parser.add_argument('--seed', type=int, default=20261103)
    parser.add_argument('--force', action='store_true', help='Replace an existing output directory.')
    args = parser.parse_args()
    started = time.perf_counter()
    source = json.loads((args.cells/'manifest.json').read_text())
    entries = {s['id']: s for s in source['slices']}
    requested = args.sections or list(entries)
    if any(s not in entries for s in requested):
        raise SystemExit(f'Sections not in {args.cells}: {[s for s in requested if s not in entries]}')
    slices = [s for s in entries if s in requested]  # the cell dataset's sample order
    um_per_unit = source['coordinates']['um_per_unit']

    genes = pd.read_parquet(args.cells/'genes.parquet').sort_values('gene_index').reset_index(drop=True)
    mapping = pd.read_parquet(args.cells/'gene_mapping.parquet')
    g = len(genes)
    cells = pd.read_parquet(args.cells/'cells.parquet', columns=['cell_index', 'slice', 'h5ad_row', 'rx', 'ry', 'n_counts', *LEVELS])
    cells = cells[cells.slice.isin(slices)].sort_values('cell_index').reset_index(drop=True)
    if args.bin_um:
        side = round(args.bin_um/um_per_unit)
    else:
        extent = max(max(p.rx.max()-p.rx.min(), p.ry.max()-p.ry.min()) for _, p in cells.groupby('slice'))
        side = math.ceil(extent/args.bins)
    if side < 1:
        raise SystemExit(f'Bins must be at least {um_per_unit:g} µm')
    bin_um = side*um_per_unit
    # One lattice for all sections, from the top left of their shared (CCF) frame. Its side is a whole number of
    # units, so bin centres are exactly one side apart and neighbouring squares meet without gaps.
    origin = (math.floor(cells.rx.min()), math.floor(cells.ry.min()))
    n_rows, n_cols = int((cells.ry.max()-origin[1])//side) + 1, int((cells.rx.max()-origin[0])//side) + 1

    out = args.out or args.cells.with_name(f'{args.cells.name}-bins{bin_um:g}um')
    if out.exists() and any(out.iterdir()) and not args.force:
        raise SystemExit(f'{out} is not empty; pass --force to replace it')
    staging = out.with_name(out.name + '.staging')
    if staging.exists():
        shutil.rmtree(staging)
    for folder in ('counts', 'ranks'):
        (staging/folder).mkdir(parents=True)
    manifest_files = []
    def record(path):
        manifest_files.append(dict(path=str(path.relative_to(staging)), bytes=path.stat().st_size, sha256=sha256(path)))

    with h5py.File(args.raw/EXPRESSION, 'r') as h5:
        if not np.array_equal(h5['var/_index'].asstr()[:], mapping.source_gene.to_numpy()):
            raise ValueError('The expression matrix must list the genes of gene_mapping.parquet, in order')
        indptr = h5['X/indptr'][:].astype(np.int64)
        indices, data = h5['X/indices'][:], h5['X/data'][:]
    remap = pd.Index(genes.symbol).get_indexer(mapping.symbol)  # source genes that share a symbol are summed
    if (remap < 0).any():
        raise ValueError('Every measured gene must be in the universe')

    rng = np.random.default_rng(args.seed)
    tables, cell_bins, slice_entries, validated, offset = [], [], [], 0, 0
    for slice_id in slices:
        t = time.perf_counter()
        part = cells[cells.slice.eq(slice_id)]
        m = len(part)
        rows, cols, values = section_counts(indptr, indices, data, part.h5ad_row.to_numpy(), remap, g)
        row = ((part.ry.to_numpy()-origin[1])//side).astype(np.int64)
        col = ((part.rx.to_numpy()-origin[0])//side).astype(np.int64)
        occupied, local = np.unique(row*n_cols + col, return_inverse=True)
        n = len(occupied)
        sums = np.bincount(local[rows]*g + cols, weights=values, minlength=n*g).reshape(n, g)
        # The same cells and genes as the cell dataset: per-cell totals, then per-gene totals of the section.
        np.testing.assert_array_equal(np.bincount(rows, weights=values, minlength=m), part.n_counts.to_numpy())
        with np.load(args.cells/entries[slice_id]['pseudobulk'], allow_pickle=False) as z:
            np.testing.assert_array_equal(sums.sum(axis=0), z['cell_class_counts'].sum(axis=0, dtype=np.float64))
        if sums.max(initial=0) >= 2**32:
            raise ValueError('Bin counts exceed uint32')
        counts = sums.astype(np.uint32)
        r, c = occupied // n_cols, occupied % n_cols
        table = pd.DataFrame(dict(slice=slice_id, bin=[f'{slice_id} r{a:02d} c{b:02d}' for a, b in zip(r+1, c+1)],
                                  row=(r+1).astype(np.uint16), col=(c+1).astype(np.uint16),
                                  rx=origin[0] + (c+.5)*side, ry=origin[1] + (r+.5)*side,
                                  n_cells=np.bincount(local).astype(np.uint32), n_counts=counts.sum(axis=1, dtype=np.int64),
                                  n_genes_detected=(counts > 0).sum(axis=1).astype(np.uint16)))
        for level in LEVELS:
            labels, codes = np.unique(part[level].to_numpy(), return_inverse=True)
            table[level], table[f'{level}_fraction'] = main_labels(codes, labels, local, n)
        nz_rows, nz_cols = np.nonzero(counts)
        indptr_ranked, gene_index, rank2, zero2 = rank_rows(nz_rows, nz_cols, counts[nz_rows, nz_cols], n, g)
        check = rng.choice(n, min(n, args.check_bins), replace=False)
        validated += validate_rows(indptr_ranked, gene_index, rank2, zero2, [(b, counts[b].astype(np.float64)) for b in check], g, rng)
        files = dict(counts=f'counts/{slice_id}.npz', ranks=f'ranks/{slice_id}.npz')
        np.savez_compressed(staging/files['counts'], counts=counts)
        # The shard format of the cell batches; its rows (cell_index) are bins.
        np.savez_compressed(staging/files['ranks'], cell_index=np.arange(offset, offset+n, dtype=np.uint32), indptr=indptr_ranked.astype(np.uint32),
                            gene_index=gene_index, rank2=rank2, zero_rank2=zero2)
        with np.load(staging/files['ranks'], allow_pickle=False) as z:  # round-trip check of the serialized ranks
            np.testing.assert_array_equal(z['gene_index'], gene_index)
            np.testing.assert_array_equal(z['rank2'], rank2)
        for name in files.values():
            record(staging/name)
        cell_bins.append(pd.DataFrame(dict(cell_index=part.cell_index.to_numpy(), bin_index=(offset + local).astype(np.uint32))))
        tables.append(table)
        slice_entries.append(dict(id=slice_id, n_cells=m, n_bins=n, bin_index_start=offset, bin_index_stop=offset+n, **files,
                                  pseudobulk=entries[slice_id]['pseudobulk'], seconds=round(time.perf_counter()-t, 1)))
        offset += n
        print(json.dumps(slice_entries[-1]), flush=True)

    bins = pd.concat(tables, ignore_index=True)
    bins.insert(0, 'bin_index', np.arange(len(bins), dtype=np.uint32))
    bins['draw_order'] = np.arange(len(bins), dtype=np.uint32)  # squares do not overlap
    bins.to_parquet(staging/'bins.parquet', index=False); record(staging/'bins.parquet')
    pd.concat(cell_bins, ignore_index=True).to_parquet(staging/'cell_bins.parquet', index=False); record(staging/'cell_bins.parquet')
    parents = {level: cells.drop_duplicates(level).set_index(level)[parent].to_dict() for level, parent in PARENT.items()}
    manifest = dict(
        format=source['format'], layout='square-bins-v2', dataset=out.name,
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source=dict(path=os.path.relpath(args.cells, out), dataset=source['dataset'], created_utc=source['created_utc'],
                    expression=dict(path=EXPRESSION, sha256=source['sources']['expression']['sha256'])),
        n_bins=len(bins), n_cells=len(cells), n_genes=g, n_slices=len(slices), seed=args.seed,
        bin_um=bin_um, bin_side=side, origin_rxry=list(origin), lattice=dict(rows=n_rows, cols=n_cols),
        observations='square bins of cells; counts are the summed transcripts of the cells whose centroid lies in the bin',
        binning='one lattice of squares (bin_side units) for all sections, from origin_rxry, the top left of their shared frame; '
                'rows from the top (dorsal), columns from the left (medial), numbered from 1; only squares with cells are bins',
        coordinates=dict(columns='rx, ry', units='micrometres in the Allen CCF (ccf_z, ccf_y), bin centre', um_per_unit=um_per_unit),
        main_labels="each bin's most common label per level among its cells (ties: first in sorted order); <level>_fraction is its share",
        ranks=dict(source['ranks'], rows='bins (cell_index holds bin_index)'),
        parents=parents,
        validation=dict(cell_count_totals_match=len(cells), gene_totals_match_cell_class_pseudobulk=len(slices),
                        exact_rank_and_pairwise_auc_rows=validated, ranks_roundtripped=True),
        software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, h5py=h5py.__version__, scipy=scipy.__version__),
        slices=slice_entries, files=manifest_files, seconds=round(time.perf_counter()-started, 1))
    (staging/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    if out.exists():
        shutil.rmtree(out)
    staging.rename(out)
    print(json.dumps({k: v for k, v in manifest.items() if k not in ['files', 'slices', 'parents']}, indent=2))
    print('OUTPUT_BYTES', sum(p.stat().st_size for p in out.rglob('*') if p.is_file()), flush=True)


if __name__ == '__main__':
    main()
