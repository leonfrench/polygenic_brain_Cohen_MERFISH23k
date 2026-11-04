"""RT&T-AMP-MERFISH mouse brain hemisections (Cohen et al. 2026) -> exact sparse ranks in per-section batches.

Input in --raw (data/, as downloaded from the Zenodo deposit, record 20682863), seven coronal hemisections:
  cellxgene_WT_labels.h5ad   the whole-transcriptome ('WT') samples 1-4 of one mouse: 146,863 segmented cells x 22,930
                             genes, raw counts (CSR; rows not grouped by sample)
  cellxgene_ISO_labels.h5ad  the isoform ('ISO') samples 5-7 of two mice: 129,563 cells x 32,684 features, raw counts;
                             the same 22,930 genes (ISO_Codebook_1 and 2) and 9,754 isoforms (ISO_Codebook_3, not used)
  Both obs: sample, cell_label (cell types from the paper's clustering), ccf_z and ccf_y (Allen CCF, micrometres), annot
  (the ABC atlas parcellation_index of the cell), global_x and global_y (experimental coordinates, micrometres), volume, num_z
The two Allen Brain Cell Atlas parcellation tables and its 25 µm annotation volume are read from --atlas (data/) and the MGI
marker lists from --reference (data/reference); all are downloaded when missing.

Every cell is kept. The universe is every measured gene, coding or not, under its current MGI symbol where that is
unambiguous; source symbols that now name the same gene are summed before ranking. The files give no anterior-posterior
coordinate, so each section's coronal plane is the one whose annotation matches most of its cells' parcellations (each
cell's was read from the annotation at its registered position); sections run anterior to posterior by that plane.

  .conda-env/bin/python scripts/prepare_rttamp_cells.py --force
"""
from pathlib import Path
import argparse, colorsys, datetime, gzip, hashlib, json, platform, shutil, struct, time
import h5py
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
MGI_URLS = {'MRK_List2.rpt': 'https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt',
            'MRK_List1.rpt': 'https://www.informatics.jax.org/downloads/reports/MRK_List1.rpt'}
ATLAS_URL = 'https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/'
ATLAS = dict(parcellation='Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_acronym.csv',
             parcellation_terms='Allen-CCF-2020/20230630/parcellation_term.csv')
RAW = ROOT/'data'  # the folder of the deposit's files, also read by prepare_rttamp_bins.py
DATASET = 'Cohen-genes'
ZENODO = 'https://zenodo.org/records/20682863'
SOURCES = {'WT': 'cellxgene_WT_labels.h5ad', 'ISO': 'cellxgene_ISO_labels.h5ad'}
# Each sample's file set and mouse (the deposit's README and the paper): the WT samples 1-4 are two anterior and two posterior
# sections of one mouse, the ISO samples 5-7 posterior sections of two mice (5 and 6 from one, 7 from the other). Neither source
# says whether an ISO mouse is the WT mouse, so mice are numbered within their file set.
SAMPLES = {**dict.fromkeys(['1', '2', '3', '4'], ('WT', 'WT mouse')), '5': ('ISO', 'ISO mouse 1'), '6': ('ISO', 'ISO mouse 1'),
           '7': ('ISO', 'ISO mouse 2')}
GENE_CODEBOOKS = ('WT_Codebook_1', 'WT_Codebook_2', 'ISO_Codebook_1', 'ISO_Codebook_2')  # ISO_Codebook_3 holds the isoforms
VOLUME = 'Allen-CCF-2020/20250331/annotation_25.nii.gz'
VOLUME_URL = 'https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/image_volumes/' + VOLUME
VOXEL_UM = 25
# The paper's groups of its cell types (STAR Methods, ligand-receptor analysis; Fig. 4G,H, where 'neurons' is 'Neu'): the WT
# clustering's 41 types and the ISO clustering's Striatum_other, which the paper puts with the striatal neurons.
GROUPS = {**dict.fromkeys(['Cortex L2/3 IT', 'Cortex L2/3 RSP', 'Cortex L4 RSP', 'Cortex L4/5 IT', 'Cortex L5 IT', 'Cortex L6 IT'],
                          'Cortex IT excitatory neurons'),
          **dict.fromkeys(['Cortex 6b', 'Cortex L5 ET', 'Cortex L5 NP', 'Cortex L6 CT'], 'Cortex non-IT excitatory neurons'),
          'Cortex Subplate & Amygdalar Areas': 'Cortical subplate neurons', 'Cortex Olfactory Areas': 'Olfactory area neurons',
          **dict.fromkeys(['Hippocampus_CA1', 'Hippocampus_CA2-CA3', 'Hippocampus_DG', 'Hippocampus_IG'], 'Hippocampus neurons'),
          'Thalamus': 'Thalamus neurons', 'Hypothalamus': 'Hypothalamus neurons', 'Midbrain': 'Midbrain neurons',
          **dict.fromkeys(['Cholinergic_STR-HY', 'Cholinergic_TH'], 'Cholinergic neurons'),
          **dict.fromkeys(['Striatum_D1', 'Striatum_D2', 'Striatum_LSr', 'Striatum_OT', 'Striatum_other'], 'Striatum neurons'),
          **dict.fromkeys(['Inhibitory_Lamp5', 'Inhibitory_Pvalb', 'Inhibitory_Sst', 'Inhibitory_Vip'], 'Inhibitory neurons'),
          # Non-neuronal types are groups of their own, as is 'Other' (a class of its own in Fig. 4A, absent from the groups).
          **{t: t for t in ['Astrocytes', 'Choroid Plexus', 'Endothelial', 'Ependymal', 'Microglia', 'NFOLs', 'OPCs', 'Oligodendrocytes',
                            'Pericytes', 'SMCs', 'VLMCs', 'Other']}}
# Fig. 4A's classes: the neuron groups are 'Neurons', every other group is its own class.
CLASSES = {group: 'Neurons' if group.endswith(' neurons') else group for group in GROUPS.values()}
# Colour families of the groups, one dark-surface categorical slot each (the plot is black), in slot order; a family's
# cell types step through its hue in lightness. Glia and 'vascular and ventricular cells' are the paper's (Fig. 2D).
FAMILIES = {'Cortex IT excitatory neurons': 'Cortex IT',
            **dict.fromkeys(['Thalamus neurons', 'Hypothalamus neurons', 'Midbrain neurons', 'Cholinergic neurons'], 'Diencephalon and midbrain'),
            'Hippocampus neurons': 'Hippocampus',
            **dict.fromkeys(['Cortex non-IT excitatory neurons', 'Cortical subplate neurons', 'Olfactory area neurons'], 'Other cortical'),
            'Inhibitory neurons': 'Inhibitory', 'Striatum neurons': 'Striatum',
            **dict.fromkeys(['Astrocytes', 'Microglia', 'NFOLs', 'OPCs', 'Oligodendrocytes'], 'Glia'),
            **dict.fromkeys(['Choroid Plexus', 'Endothelial', 'Ependymal', 'Pericytes', 'SMCs', 'VLMCs'], 'Vascular and ventricular'),
            'Other': 'Other'}
# Categorical slots stepped for a dark surface (the dataviz reference palette); 'Other' is the muted grey.
PALETTE = ['#3987E5', '#D95926', '#199E70', '#C98500', '#D55181', '#008300', '#9085E9', '#E66767']
UNKNOWN = '#898781'
FAMILY_COLOURS = {**dict(zip(dict.fromkeys(f for f in FAMILIES.values() if f != 'Other'), PALETTE)), 'Other': UNKNOWN}
# OKLCH lightness of the type steps: clear of the black plot, wide enough for families of up to six types.
TYPE_BAND = (0.50, 0.80)
MIN_CONTRAST = 3.0  # against the black plot, as WCAG asks of graphical objects
# Names for the abbreviations, from the paper's text (shown as 'name (label)' in the tooltip).
NAMES = {'NFOLs': 'Newly formed oligodendrocytes', 'OPCs': 'Oligodendrocyte precursor cells', 'SMCs': 'Smooth muscle cells',
         'VLMCs': 'Vascular and leptomeningeal cells'}
CELL_LEVELS = {'cell_class': 'Cell class', 'cell_group': 'Cell group', 'cell_type': 'Cell type'}
REGION_LEVELS = {'region_division': ('division', 'Brain division'), 'region_structure': ('structure', 'Brain structure'),
                 'region_substructure': ('substructure', 'Brain substructure')}
LEVELS = [*CELL_LEVELS, *REGION_LEVELS]
UNASSIGNED = '#C9CED6'  # the atlas colour for unassigned tissue is black, the plot background


def sha256(path):
    with open(path, 'rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def download(url, path):
    if not path.exists():
        import urllib.request
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(path.suffix + '.part')
        urllib.request.urlretrieve(url, tmp)
        tmp.replace(path)
    return path


# ---------------------------------------------------------------- gene universe

def mgi_tables(reference):
    current = pd.read_csv(download(MGI_URLS['MRK_List2.rpt'], reference/'MRK_List2.rpt'), sep='\t', dtype=str, keep_default_na=False)
    current = current[current['Status'].eq('O') & current['Marker Type'].eq('Gene')]
    pc = current[current['Feature Type'].eq('protein coding gene')][['MGI Accession ID', 'Marker Symbol', 'Marker Name']]
    pc.columns = ['mgi_id', 'symbol', 'name']
    full = pd.read_csv(download(MGI_URLS['MRK_List1.rpt'], reference/'MRK_List1.rpt'), sep='\t', dtype=str, keep_default_na=False)
    withdrawn = full[full['Status'].eq('W')][['Marker Symbol', 'Current MGI Accession ID (if withdrawn)']]
    withdrawn.columns = ['old_symbol', 'mgi_id']
    return current, pc.reset_index(drop=True), withdrawn


def map_genes(names, current, pc, withdrawn):
    """Current MGI symbol first, then an unambiguous withdrawn symbol. No synonym guessing."""
    by_symbol = dict(zip(current['Marker Symbol'], current['MGI Accession ID']))
    pc_ids = dict(zip(pc.mgi_id, pc.symbol))
    old = withdrawn[withdrawn.mgi_id.ne('')].groupby('old_symbol').mgi_id.agg(lambda s: sorted(set(s)))
    rows = []
    for name in names:
        ident, method = by_symbol.get(name), 'current_symbol'
        if ident is None:
            candidates = old.get(name, [])
            method = 'withdrawn_symbol' if candidates else 'unmapped'
            ident = candidates[0] if len(candidates) == 1 else None
            if len(candidates) > 1:
                method = 'ambiguous_withdrawn_symbol'
        status = 'retained' if ident in pc_ids else ('not_protein_coding' if ident else method)
        rows.append((name, ident, pc_ids.get(ident), method, status))
    return pd.DataFrame(rows, columns=['source_gene', 'mgi_id', 'symbol', 'mapping_method', 'status'])


def gene_universe(source, codebooks, current, pc, withdrawn):
    """Each measured gene's universe symbol and the sorted universe.

    A gene takes its current MGI symbol (directly, or from an unambiguous withdrawn symbol); a symbol MGI
    lacks takes the one current symbol it matches apart from case (ETV3L -> Etv3l), else stays as measured.
    Coding and non-coding genes are kept alike: the panel already holds only protein-coding and lincRNA genes.
    """
    mapping = map_genes(list(source), current, pc, withdrawn).drop(columns=['symbol', 'status'])
    by_case = current.groupby(current['Marker Symbol'].str.lower())['MGI Accession ID'].agg(list)
    for i in mapping.index[mapping.mapping_method.eq('unmapped')]:
        ids = by_case.get(mapping.source_gene[i].lower(), [])
        if len(ids) == 1:
            mapping.loc[i, ['mgi_id', 'mapping_method']] = ids[0], 'case_insensitive_symbol'
    by_id = current.set_index('MGI Accession ID')
    mapping['symbol'] = mapping.mgi_id.map(by_id['Marker Symbol']).fillna(mapping.source_gene)
    mapping['name'] = mapping.mgi_id.map(by_id['Marker Name']).fillna('')
    mapping['feature_type'] = mapping.mgi_id.map(by_id['Feature Type']).fillna('')
    mapping['codebook'] = np.asarray(codebooks)
    symbols = np.array(sorted(set(mapping.symbol)))
    # The app matches symbols case-insensitively.
    if pd.Series(symbols).str.upper().duplicated().any():
        raise ValueError('Universe symbols must differ in more than case')
    return mapping, symbols


# ---------------------------------------------------------------- ranks and pseudobulks

def rank_rows(rows, cols, values, n_rows, n_genes):
    """Doubled average ranks per row; positive values only, zeros share (Z+1)/2."""
    rows, cols, values = (np.asarray(a, dtype=np.int64) for a in (rows, cols, values))
    # One integer key per sort; np.lexsort is much slower on unordered rows.
    order = np.argsort(rows*(int(values.max(initial=0))+1) + values)
    rows, cols, values = rows[order], cols[order], values[order]
    counts = np.bincount(rows, minlength=n_rows)
    zeros = n_genes - counts
    start = np.r_[0, np.flatnonzero((np.diff(rows) != 0) | (np.diff(values) != 0)) + 1]
    end = np.r_[start[1:], len(rows)]
    row_start = np.r_[0, np.cumsum(counts)[:-1]]
    # A 0-based [start,end) tie inside the positive block has doubled 1-based rank start+end+1.
    rank2 = np.repeat(start + end + 1 - 2*row_start[rows[start]] + 2*zeros[rows[start]], end - start)
    order = np.argsort(rows*n_genes + cols)
    indptr = np.r_[0, np.cumsum(counts)].astype(np.int64)
    return indptr, cols[order].astype(np.uint16), rank2[order].astype(np.uint16), np.where(zeros > 0, zeros+1, 0).astype(np.uint16)


def pseudobulk(rows, cols, values, codes, n_genes):
    """Summed counts per (label, gene) for the labels present; rows index cells, codes give each cell's label."""
    present, inverse = np.unique(codes, return_inverse=True)
    sums = np.bincount(inverse[rows]*n_genes + cols, weights=values, minlength=len(present)*n_genes)
    if sums.max(initial=0) >= 2**32:
        raise ValueError('Pseudobulk counts exceed uint32')
    return present, np.bincount(inverse, minlength=len(present)), sums.reshape(len(present), n_genes).astype(np.uint32)


def validate_rows(indptr, cols, rank2, zero2, dense_rows, n_genes, rng):
    checked = 0
    for row, dense in dense_rows:
        a, b = indptr[row], indptr[row+1]
        decoded = np.full(n_genes, zero2[row]/2); decoded[cols[a:b]] = rank2[a:b]/2
        np.testing.assert_array_equal(decoded, rankdata(dense, method='average'))
        selected = rng.choice(n_genes, 100, replace=False)
        other = np.setdiff1d(np.arange(n_genes), selected)
        truth = np.mean((dense[selected, None] > dense[other]) + 0.5*(dense[selected, None] == dense[other]))
        auc = (decoded[selected].sum() - 100*101/2) / (100*(n_genes-100))
        np.testing.assert_allclose(auc, truth, atol=1e-12)
        checked += 1
    return checked


# ---------------------------------------------------------------- colours

def hex_rgb(color):
    return tuple(int(color[i:i+2], 16)/255 for i in (1, 3, 5))


def separate_colors(colors):
    """Spread lightness inside groups of (near-)identical colours, keeping each hue family."""
    labels = sorted(colors)
    groups, assigned = [], {}
    for label in labels:
        rgb = np.array(hex_rgb(colors[label]))
        for group in groups:
            if np.abs(rgb-group['rgb']).max() < 12/255:
                group['labels'].append(label); break
        else:
            groups.append(dict(rgb=rgb, labels=[label]))
    for group in groups:
        m = len(group['labels'])
        for i, label in enumerate(group['labels']):
            h, l, s = colorsys.rgb_to_hls(*hex_rgb(colors[label]))
            if m > 1:
                l = float(np.clip(l + np.linspace(-1, 1, m)[i]*min(0.22, 0.09*(m-1)), 0.12, 0.88))
            r, g, b = colorsys.hls_to_rgb(h, l, s)
            assigned[label] = '#%02X%02X%02X' % tuple(round(v*255) for v in (r, g, b))
    return assigned


def oklab(color):
    rgb = np.array([int(color[i:i+2], 16)/255 for i in (1, 3, 5)])
    rgb = np.where(rgb <= 0.04045, rgb/12.92, ((rgb + 0.055)/1.055)**2.4)
    lms = np.cbrt(np.array([[0.4122214708, 0.5363325363, 0.0514459929], [0.2119034982, 0.6806995451, 0.1073969566],
                            [0.0883024619, 0.2817188376, 0.6299787005]]) @ rgb)
    return np.array([[0.2104542553, 0.7936177850, -0.0040720468], [1.9779984951, -2.4285922050, 0.4505937099],
                     [0.0259040371, 0.7827717662, -0.8086757660]]) @ lms


def from_oklab(lab):
    """sRGB hex, or None when outside the gamut."""
    lms = (np.array([[1, 0.3963377774, 0.2158037573], [1, -0.1055613458, -0.0638541728], [1, -0.0894841775, -1.2914855480]]) @ lab)**3
    rgb = np.array([[4.0767416621, -3.3077115913, 0.2309699292], [-1.2684380046, 2.6097574011, -0.3413193965],
                    [-0.0041960863, -0.7034186147, 1.7076147010]]) @ lms
    if (rgb < -1e-6).any() or (rgb > 1 + 1e-6).any():
        return None
    rgb = np.where(rgb <= 0.0031308, 12.92*rgb, 1.055*np.clip(rgb, 0, 1)**(1/2.4) - 0.055)
    return '#%02X%02X%02X' % tuple(int(round(v*255)) for v in np.clip(rgb, 0, 1))


def subtype_colours(classes, colours, band=TYPE_BAND):
    """Each class's hue at evenly spaced OKLCH lightness across the band, one step per subtype (sorted).

    A lone subtype keeps its class colour; chroma is reduced where a step leaves the sRGB gamut.
    """
    out = {}
    for parent in sorted(set(classes.values())):
        members = sorted(k for k, v in classes.items() if v == parent)
        if len(members) == 1:
            out[members[0]] = colours[parent]; continue
        _, a, b = oklab(colours[parent])
        for member, lightness in zip(members, np.linspace(*band, len(members))):
            out[member] = next(c for c in (from_oklab(np.array([lightness, a*f, b*f])) for f in np.linspace(1, 0, 101)) if c)
    return out


def contrast(colour):
    """Contrast ratio of a colour against black."""
    rgb = np.array([int(colour[i:i+2], 16)/255 for i in (1, 3, 5)])
    luminance = np.where(rgb <= 0.04045, rgb/12.92, ((rgb + 0.055)/1.055)**2.4) @ [0.2126, 0.7152, 0.0722]
    return float((luminance + 0.05)/0.05)


# ---------------------------------------------------------------- cells and annotations

def categorical(h5, column, group='obs'):
    """A column stored as categories and codes, as strings (the samples are integer categories)."""
    codes = h5[f'{group}/{column}/codes'][:]
    if (codes < 0).any():
        raise ValueError(f'{column} has missing values')
    categories = h5[f'{group}/{column}/categories']
    return (categories.asstr()[:] if categories.dtype.kind == 'O' else categories[:].astype(str))[codes]


def read_matrix(path, mapping, symbols):
    """One file's counts (CSR), each column's universe index (-1 for features outside it, the isoforms) and each row's gene total."""
    with h5py.File(path, 'r') as h5:
        if h5['X'].attrs.get('encoding-type') != 'csr_matrix':
            raise ValueError(f'{path.name}: expected a CSR count matrix')
        source = h5['var/_index'].asstr()[:]
        gene = np.isin(categorical(h5, 'codebook', 'var'), GENE_CODEBOOKS)
        indptr = h5['X/indptr'][:].astype(np.int64)
        indices, data = h5['X/indices'][:], h5['X/data'][:]
    if sorted(source[gene]) != sorted(mapping.source_gene):
        raise ValueError(f'{path.name} must measure the genes of the universe')
    remap = np.full(len(source), -1, np.int64)
    remap[gene] = pd.Index(symbols).get_indexer(mapping.set_index('source_gene').symbol.reindex(source[gene]))
    if not (np.diff(indptr) > 0).all():
        raise ValueError(f'{path.name}: every cell needs counts')
    totals = np.add.reduceat(np.where(remap[indices] >= 0, data, 0), indptr[:-1])  # exact, as no row is empty
    return indptr, indices, data, remap, totals


def section_counts(indptr, indices, data, rows, remap, g):
    """(row, gene_index, count) entries of the given matrix rows (row = position in rows), sorted by row and gene.

    Source genes that share a universe symbol are summed; features outside the universe (remap -1) are dropped.
    """
    starts, lengths = indptr[rows], indptr[rows+1] - indptr[rows]
    gather = np.repeat(starts - np.r_[0, np.cumsum(lengths)[:-1]], lengths) + np.arange(lengths.sum())
    values, cols = data[gather], remap[indices[gather]]
    if (values <= 0).any() or not np.array_equal(values, np.round(values)):
        raise ValueError('Stored counts must be positive integers')
    keep = cols >= 0
    key, inverse = np.unique(np.repeat(np.arange(len(rows)), lengths)[keep]*g + cols[keep], return_inverse=True)
    return key // g, key % g, np.bincount(inverse, weights=values[keep]).astype(np.int64)


def read_annotation(path, shape=(528, 320, 456)):
    """The 25 µm CCF parcellation volume (a parcellation_index per voxel) as [AP, DV, ML], from anterior, dorsal and left.

    A NIfTI-1 file whose first voxel axis runs posterior and second ventral, as its srow affine must confirm.
    """
    raw = gzip.decompress(Path(path).read_bytes())
    dim = struct.unpack('<8h', raw[40:56]); datatype, = struct.unpack('<h', raw[70:72]); offset, = struct.unpack('<f', raw[108:112])
    srow = np.array(struct.unpack('<12f', raw[280:328])).reshape(3, 4)
    if raw[344:347] != b'n+1' or dim[:4] != (3, *shape) or datatype != 512 or not (srow[1, 0] < 0 and srow[2, 1] < 0):
        raise ValueError('Expected the CCF annotation as a uint16 NIfTI-1 volume ordered anterior to posterior and dorsal to ventral')
    return np.frombuffer(raw, '<u2', count=int(np.prod(shape)), offset=int(offset)).reshape(shape, order='F')


def section_plane(volume, annot, y, z):
    """The coronal plane (VOXEL_UM index, anterior first) whose annotation matches most of a section's cells, and that fraction.

    Each cell's parcellation was read from the annotation at its registered position, so the plane the section was registered
    to matches nearly every cell, and the match falls off within a few planes.
    """
    yi = np.clip((y // VOXEL_UM).astype(np.int64), 0, volume.shape[1]-1)
    zi = np.clip((z // VOXEL_UM).astype(np.int64), 0, volume.shape[2]-1)
    match = (volume[:, yi, zi] == annot).mean(axis=1)
    return int(match.argmax()), float(match.max())


def majority(groups, members):
    """Each group's most common member (ties go to the first name)."""
    return pd.crosstab(np.asarray(groups), np.asarray(members)).idxmax(axis=1).to_dict()


def cell_levels(cells):
    """Cell class, group and type labels with their colours and names.

    The types of a colour family step through its hue in lightness; a group or class takes the colour of
    its most numerous type.
    """
    unknown = sorted(set(cells.cell_type) - GROUPS.keys())
    if unknown:
        raise ValueError(f'Cell types without a group: {unknown}')
    colours = subtype_colours({t: FAMILIES[g] for t, g in GROUPS.items()}, FAMILY_COLOURS)
    levels = {}
    for column, title in CELL_LEVELS.items():
        present = sorted(set(cells[column]))
        top = {k: k for k in present} if column == 'cell_type' else majority(cells[column], cells.cell_type)
        levels[column] = dict(kind='cell_type', title=title, labels={k: colours[top[k]] for k in present},
                              descriptions={k: NAMES[k] for k in present if k in NAMES})
    dark = {f'{label} ({level})': c for level, v in levels.items() for label, c in v['labels'].items() if contrast(c) < MIN_CONTRAST}
    if dark:
        raise ValueError(f'Colours below {MIN_CONTRAST}:1 contrast with the black plot: {dark}')
    return levels


def region_levels(cells, terms):
    """CCF division, structure and substructure labels with Allen colours and names."""
    # parcellation_term lists an acronym once per level it appears in; the copies agree.
    terms = terms.drop_duplicates(['acronym', 'name', 'color_hex_triplet'])
    if terms.acronym.duplicated().any():
        raise ValueError('Parcellation acronyms must have one name and colour')
    terms = terms.set_index('acronym')
    levels = {}
    for column, (_, title) in REGION_LEVELS.items():
        present = sorted(set(cells[column]))
        colors = {k: UNASSIGNED if k == 'unassigned' else terms.color_hex_triplet[k].upper() for k in present}
        descriptions = {k: 'Outside the CCF parcellation' if k == 'unassigned' else terms.name[k] for k in present}
        levels[column] = dict(kind='region', title=title, labels=separate_colors(colors), descriptions=descriptions)
    return levels


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--raw', type=Path, default=RAW, help='Folder of the deposit\'s h5ad files.')
    parser.add_argument('--atlas', type=Path, default=ROOT/'data', help='Folder of the Allen Brain Cell Atlas tables and volume, in the bucket layout.')
    parser.add_argument('--out', type=Path, default=ROOT/'data/rttamp-genes')
    parser.add_argument('--reference', type=Path, default=ROOT/'data/reference')
    parser.add_argument('--batch-cells', type=int, default=16384, help='Maximum cells per batch.')
    parser.add_argument('--batch-entries', type=int, default=12_000_000,
                        help='Maximum stored (cell, detected gene) entries per batch; keeps files well under 95 MB.')
    parser.add_argument('--seed', type=int, default=20261103)
    parser.add_argument('--force', action='store_true', help='Replace an existing output directory.')
    args = parser.parse_args()
    started = time.perf_counter()
    if args.out.exists() and any(args.out.iterdir()) and not args.force:
        raise SystemExit(f'{args.out} is not empty; pass --force to replace it')
    staging = args.out.with_name(args.out.name + '.staging')
    if staging.exists():
        shutil.rmtree(staging)
    (staging/'progressive').mkdir(parents=True)
    (staging/'pseudobulk').mkdir()

    tables = {name: download(ATLAS_URL+path, args.atlas/path) for name, path in ATLAS.items()}
    volume = read_annotation(download(VOLUME_URL, args.atlas/VOLUME))
    parcellation = pd.read_csv(tables['parcellation'], keep_default_na=False).set_index('parcellation_index')
    terms = pd.read_csv(tables['parcellation_terms'], dtype=str, keep_default_na=False)
    current, pc, withdrawn = mgi_tables(args.reference)
    paths = {name: args.raw/file for name, file in SOURCES.items()}
    with h5py.File(paths['WT'], 'r') as h5:  # the universe is the WT file's genes, which the ISO file measures too
        mapping, symbols = gene_universe(h5['var/_index'].asstr()[:], categorical(h5, 'codebook', 'var'), current, pc, withdrawn)
    g = len(symbols)
    if not 1 < g < 32768:
        raise ValueError('uint16 doubled ranks require 2 <= genes <= 32767')
    frames, matrices = [], {}
    for name, path in paths.items():
        with h5py.File(path, 'r') as h5:
            annot = h5['obs/annot'][:]  # stored as floats in the ISO file
            if not np.array_equal(annot, np.round(annot)):
                raise ValueError(f'{path.name}: parcellation indices must be whole numbers')
            frames.append(pd.DataFrame({'slice': categorical(h5, 'sample'), 'source': name, 'cell_id': h5['obs/_index'].asstr()[:],
                                        'h5ad_row': np.arange(len(annot)), 'rx': h5['obs/ccf_z'][:], 'ry': h5['obs/ccf_y'][:],
                                        'x': h5['obs/global_x'][:], 'y': h5['obs/global_y'][:], 'volume': h5['obs/volume'][:],
                                        'num_z': h5['obs/num_z'][:], 'parcellation_index': annot.astype(np.int64),
                                        'cell_type': categorical(h5, 'cell_label')}))
        matrices[name] = read_matrix(path, mapping, symbols)
        if len(matrices[name][0]) != len(frames[-1]) + 1:
            raise ValueError(f'{path.name}: every cell needs a matrix row')
    cells = pd.concat(frames, ignore_index=True)
    if set(zip(cells.slice, cells.source)) != {(s, source) for s, (source, _) in SAMPLES.items()}:
        raise ValueError(f'Expected samples {sorted(SAMPLES)} in their file sets')
    if cells.cell_id.duplicated().any() or not np.isfinite(cells[['rx', 'ry', 'x', 'y']].to_numpy()).all():
        raise ValueError('Cells need unique names and finite coordinates')
    if not parcellation.index.is_unique:
        raise ValueError('Parcellation indices must be unique')
    unknown = sorted(set(cells.parcellation_index) - set(parcellation.index))
    if unknown:
        raise ValueError(f'Parcellation indices missing from the atlas table: {unknown[:10]}')
    # Each section's registered coronal plane; sections run anterior to posterior, ties (adjacent sections) in sample order.
    planes = {s: section_plane(volume, p.parcellation_index.to_numpy(), p.ry.to_numpy(), p.rx.to_numpy()) for s, p in cells.groupby('slice')}
    slices = sorted(SAMPLES, key=lambda s: (planes[s][0], int(s)))
    rank = cells.slice.map({s: i for i, s in enumerate(slices)}).to_numpy()
    cells = cells.iloc[np.lexsort((cells.h5ad_row.to_numpy(), rank))].reset_index(drop=True)  # source order within each section
    rank = np.sort(rank)
    n = len(cells)
    cells.insert(0, 'cell_index', np.arange(n, dtype=np.uint32))
    cells.insert(cells.columns.get_loc('cell_type'), 'cell_group', cells.cell_type.map(GROUPS))
    cells.insert(cells.columns.get_loc('cell_group'), 'cell_class', cells.cell_group.map(CLASSES))
    for column, (source_level, _) in REGION_LEVELS.items():
        cells[column] = cells.parcellation_index.map(parcellation[source_level]).to_numpy()
    cells['draw_order'] = np.random.default_rng(args.seed).permutation(n).astype(np.uint32)
    cells['batch'] = ''; cells['batch_row'] = np.zeros(n, np.int32)
    cells['n_counts'] = np.zeros(n, np.int64); cells['n_genes_detected'] = np.zeros(n, np.int64)
    levels = {**cell_levels(cells), **region_levels(cells, terms)}
    bounds = {s: (int(np.searchsorted(rank, i)), int(np.searchsorted(rank, i, side='right'))) for i, s in enumerate(slices)}

    # Exact ranks, randomly ordered per-section batches and per-section pseudobulk counts.
    label_index = {level: np.unique(cells[level].to_numpy()) for level in LEVELS}
    rng = np.random.default_rng(args.seed)
    detection = np.zeros(g, np.int64); totals = np.zeros(g, np.float64); per_slice = {}
    manifest_files, slice_entries, validated = [], [], 0
    def record(path, **extra):
        item = dict(path=str(path.relative_to(staging)), bytes=path.stat().st_size, sha256=sha256(path), **extra)
        if item['bytes'] >= 95_000_000:
            raise ValueError(f'File exceeds conservative GitHub limit: {path}')
        manifest_files.append(item); return item
    for slice_id in slices:
        t = time.perf_counter()
        start, stop = bounds[slice_id]
        part = cells.iloc[start:stop]
        m = len(part)
        source, mouse = SAMPLES[slice_id]
        indptr, indices, data, remap, source_totals = matrices[source]
        rows, cols, values = section_counts(indptr, indices, data, part.h5ad_row.to_numpy(), remap, g)
        n_counts = np.bincount(rows, weights=values, minlength=m).astype(np.int64)
        if not np.array_equal(n_counts, source_totals[part.h5ad_row.to_numpy()]):  # merging symbols keeps every gene count
            raise ValueError(f'{slice_id}: cell totals changed when genes were merged')
        if (n_counts <= 0).any():
            raise ValueError(f'{slice_id}: every cell needs gene counts')
        cells.loc[part.index, 'n_counts'] = n_counts
        cells.loc[part.index, 'n_genes_detected'] = np.bincount(rows, minlength=m)
        per_slice[slice_id] = np.bincount(cols, minlength=g)
        detection += per_slice[slice_id]; totals += np.bincount(cols, weights=values, minlength=g)
        indptr_ranked, gene_index, rank2, zero2 = rank_rows(rows, cols, values, m, g)
        dense_rows = []
        for row in rng.choice(m, min(m, 32), replace=False):
            a, b = np.searchsorted(rows, [row, row+1])
            dense = np.zeros(g); dense[cols[a:b]] = values[a:b]; dense_rows.append((row, dense))
        validated += validate_rows(indptr_ranked, gene_index, rank2, zero2, dense_rows, g, rng)
        shuffled = rng.permutation(m)
        # Equal-size batches, small enough for both the cell and stored-entry limits.
        size = max(1, min(args.batch_cells, int(args.batch_entries*m/max(1, len(gene_index)))))
        size = -(-m // -(-m // size))
        batches = []
        for b, first in enumerate(range(0, m, size)):
            chosen = shuffled[first:first+size]
            lengths = indptr_ranked[chosen+1] - indptr_ranked[chosen]
            ptr = np.r_[0, np.cumsum(lengths)]
            gather = np.repeat(indptr_ranked[chosen] - ptr[:-1], lengths) + np.arange(ptr[-1])
            path = staging/'progressive'/slice_id/f'{b:05d}.npz'; path.parent.mkdir(exist_ok=True)
            global_rows = (start + chosen).astype(np.uint32)
            np.savez_compressed(path, cell_index=global_rows, indptr=ptr.astype(np.uint32),
                                gene_index=gene_index[gather], rank2=rank2[gather], zero_rank2=zero2[chosen])
            with np.load(path, allow_pickle=False) as z:  # round-trip check of the serialized batch
                np.testing.assert_array_equal(z['gene_index'], gene_index[gather])
                np.testing.assert_array_equal(z['rank2'], rank2[gather])
            cells.loc[global_rows, 'batch'] = str(path.relative_to(staging))
            cells.loc[global_rows, 'batch_row'] = np.arange(len(chosen), dtype=np.int32)
            batches.append(record(path, slice=slice_id, n_cells=len(chosen))['path'])
        arrays = {}
        for level in LEVELS:
            present, n_cells, sums = pseudobulk(rows, cols, values, np.searchsorted(label_index[level], part[level].to_numpy()), g)
            np.testing.assert_array_equal(sums.sum(axis=0, dtype=np.float64), np.bincount(cols, weights=values, minlength=g))
            arrays.update({f'{level}_index': present.astype(np.uint16), f'{level}_cells': n_cells.astype(np.uint32),
                           f'{level}_counts': sums})
        path = staging/'pseudobulk'/f'{slice_id}.npz'
        np.savez_compressed(path, **arrays)
        record(path, slice=slice_id, layout='pseudobulk')
        plane, match = planes[slice_id]
        slice_entries.append(dict(id=slice_id, sample=int(slice_id), source=source, mouse=mouse, n_cells=m,
                                  ccf_plane=plane, ccf_x_um=plane*VOXEL_UM + VOXEL_UM/2, plane_match=round(match, 4),
                                  cell_index_start=start, cell_index_stop=stop,
                                  bbox_rxry=[float(part.rx.min()), float(part.ry.min()), float(part.rx.max()), float(part.ry.max())],
                                  unassigned_cells=int(part.region_substructure.eq('unassigned').sum()),
                                  batches=batches, pseudobulk=str(path.relative_to(staging)), rank_seconds=round(time.perf_counter()-t, 2)))
        print(f'sample {slice_id} ({source}, plane {(plane + 0.5)*VOXEL_UM/1000:.2f} mm, {match:.1%} matched): {m:,} cells ranked into '
              f'{len(batches)} batches in {time.perf_counter()-t:.1f}s', flush=True)

    genes = (mapping.groupby('symbol').agg(mgi_id=('mgi_id', 'first'), name=('name', 'first'), feature_type=('feature_type', 'first'),
                                           codebook=('codebook', lambda s: '|'.join(sorted(set(s)))),
                                           source_genes=('source_gene', '|'.join))
             .reindex(symbols).reset_index())
    genes.insert(0, 'gene_index', np.arange(g, dtype=np.uint16))
    genes['n_cells_detected'] = detection; genes['total_counts'] = totals
    (staging/'annotations.json').write_text(json.dumps(levels, indent=1, ensure_ascii=False)+'\n'); record(staging/'annotations.json')
    cells.to_parquet(staging/'cells.parquet', index=False); record(staging/'cells.parquet')
    genes.to_parquet(staging/'genes.parquet', index=False); record(staging/'genes.parquet')
    mapping.to_parquet(staging/'gene_mapping.parquet', index=False); record(staging/'gene_mapping.parquet')
    pd.DataFrame({'gene_index': np.tile(np.arange(g, dtype=np.uint16), len(slices)), 'slice': np.repeat(slices, g),
                  'n_cells_detected': np.concatenate([per_slice[s] for s in slices])}).to_parquet(
        staging/'gene_detection.parquet', index=False); record(staging/'gene_detection.parquet')
    manifest = dict(
        format='sparse-doubled-average-ranks-v1', layout='per-slice-random-batches-v1', dataset=DATASET,
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        n_cells=n, n_genes=g, n_slices=len(slices), seed=args.seed, batch_cells=args.batch_cells, batch_entries=args.batch_entries,
        observations='RT&T-AMP-MERFISH segmented cells of seven coronal hemisections: the deposit\'s whole-transcriptome '
                     'samples 1-4 (cellxgene_WT_labels.h5ad) and isoform samples 5-7 (cellxgene_ISO_labels.h5ad, gene features only), '
                     'as annotated there; counts are decoded transcripts per cell',
        coordinates=dict(columns='rx, ry', units='micrometres in the Allen CCF: rx = ccf_z (medial-lateral, midline at 5,700), '
                                                  'ry = ccf_y (dorsal-ventral, dorsal up)', um_per_unit=1.0,
                         raw='x, y are global_x, global_y, the experimental coordinates (micrometres)'),
        section_order=f'anterior to posterior by ccf_x_um, the centre of the {VOXEL_UM} µm CCF coronal plane whose annotation matches '
                      'most of the section\'s cell parcellations (plane_match); ties keep sample order',
        region_assignment='ABC atlas parcellation (Allen-CCF-2020) of each cell\'s registered position, from obs annot '
                          '(its parcellation_index); unassigned = outside the parcellation',
        sources={'expression': {name: dict(path=file, sha256=sha256(paths[name]), url=ZENODO) for name, file in SOURCES.items()},
                 'paper': 'Cohen et al., Cell 189, 6189-6207 (2026), doi:10.1016/j.cell.2026.06.027',
                 **{name: dict(url=ATLAS_URL+path, sha256=sha256(tables[name])) for name, path in ATLAS.items()},
                 'annotation_volume': dict(url=VOLUME_URL, sha256=sha256(args.atlas/VOLUME)),
                 'mgi': {k: dict(url=v, sha256=sha256(args.reference/k)) for k, v in MGI_URLS.items()}},
        gene_mapping_methods=mapping.mapping_method.value_counts().to_dict(),
        feature_types=genes.feature_type.replace('', 'not in MGI').value_counts().to_dict(),
        merged_source_genes=int(len(mapping) - g),
        ranks=dict(direction='ascending: larger rank means higher expression', ties='average', scale=2,
                   zeros='all unstored genes have rank zero_rank2 / 2; stored ranks are absolute',
                   universe='every measured gene (those of the WT file, which the ISO file measures too), under its current MGI symbol '
                            'where unambiguous; the ISO isoform features are not used',
                   duplicate_genes='sum counts before ranking', auroc='(sum(rank of set) - k*(k+1)/2) / (k*(G-k))'),
        annotation_levels={k: dict(kind=v['kind'], title=v['title'], n_labels=len(v['labels'])) for k, v in levels.items()},
        pseudobulk=dict(description='per section and annotation level: summed counts per (label, gene), '
                                    '<level>_counts rows follow <level>_index into labels[level]; <level>_cells counts cells',
                        labels={level: label_index[level].tolist() for level in LEVELS}, cells={level: 'all' for level in LEVELS}),
        validation=dict(exact_rank_and_pairwise_auc_rows=validated, batches_roundtripped=True),
        software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
                      h5py=h5py.__version__),
        slices=slice_entries, files=manifest_files, seconds=round(time.perf_counter()-started, 1))
    (staging/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    if args.out.exists():
        shutil.rmtree(args.out)
    staging.rename(args.out)
    print(json.dumps({k: v for k, v in manifest.items() if k not in ['files', 'slices', 'pseudobulk']}, indent=2))
    print('OUTPUT_BYTES', sum(p.stat().st_size for p in args.out.rglob('*') if p.is_file()), flush=True)


if __name__ == '__main__':
    main()
