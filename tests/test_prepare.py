import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scipy.stats import rankdata

# Preprocessing needs requirements-prep.txt; an app-only install (requirements.txt) skips these tests.
h5py = pytest.importorskip('h5py')
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_rttamp_cells import (CLASSES, GROUPS, MIN_CONTRAST, NAMES, TYPE_BAND, UNKNOWN, categorical, cell_levels,  # noqa: E402
                                  contrast, gene_universe, map_genes, oklab, pseudobulk, rank_rows, section_counts, separate_colors)


def mgi():
    current = pd.DataFrame({'Marker Symbol': ['Exoc2', 'Ca4', 'Rmst', 'Etv3l'], 'MGI Accession ID': ['MGI:1', 'MGI:2', 'MGI:3', 'MGI:4'],
                            'Marker Name': ['exocyst complex component 2', 'carbonic anhydrase 4', 'rhabdomyosarcoma 2 associated transcript',
                                            'ets variant 3-like'],
                            'Feature Type': ['protein coding gene', 'protein coding gene', 'lncRNA gene', 'protein coding gene']})
    pc = pd.DataFrame({'mgi_id': ['MGI:1', 'MGI:2', 'MGI:4'], 'symbol': ['Exoc2', 'Ca4', 'Etv3l'], 'name': ['', '', '']})
    withdrawn = pd.DataFrame({'old_symbol': ['Gm29675', 'Car4', 'Lhfp', 'Lhfp'], 'mgi_id': ['MGI:1', 'MGI:2', 'MGI:5', 'MGI:6']})
    return current, pc, withdrawn


def test_universe_keeps_every_gene_under_its_current_symbol():
    source = ['Exoc2', 'Gm29675', 'Car4', 'Rmst', 'ETV3L', 'Lhfp', 'AC103362.1']
    mapping, symbols = gene_universe(source, ['WT_Codebook_1']*6 + ['WT_Codebook_2'], *mgi())
    table = mapping.set_index('source_gene')
    assert symbols.tolist() == ['AC103362.1', 'Ca4', 'Etv3l', 'Exoc2', 'Lhfp', 'Rmst']  # two source genes now name Exoc2
    assert table.symbol.to_dict() == {'Exoc2': 'Exoc2', 'Gm29675': 'Exoc2', 'Car4': 'Ca4', 'Rmst': 'Rmst', 'ETV3L': 'Etv3l',
                                      'Lhfp': 'Lhfp', 'AC103362.1': 'AC103362.1'}
    assert table.mapping_method.to_dict() == {'Exoc2': 'current_symbol', 'Gm29675': 'withdrawn_symbol', 'Car4': 'withdrawn_symbol',
                                              'Rmst': 'current_symbol', 'ETV3L': 'case_insensitive_symbol',
                                              'Lhfp': 'ambiguous_withdrawn_symbol', 'AC103362.1': 'unmapped'}
    # Non-coding genes stay in the universe; genes MGI lacks keep their measured symbol and have no name.
    assert table.loc['Rmst', 'feature_type'] == 'lncRNA gene' and table.loc['Car4', 'name'] == 'carbonic anhydrase 4'
    assert table.loc['AC103362.1', ['name', 'feature_type', 'codebook']].tolist() == ['', '', 'WT_Codebook_2']
    with pytest.raises(ValueError, match='case'):  # the app matches symbols case-insensitively
        gene_universe(['Xyz1', 'XYZ1'], ['WT_Codebook_1']*2, *mgi())


def test_gene_mapping_prefers_current_symbols_and_rejects_ambiguous_withdrawn():
    current = pd.DataFrame({'Marker Symbol': ['Mbp', 'Gm1'], 'MGI Accession ID': ['MGI:1', 'MGI:2']})
    pc = pd.DataFrame({'mgi_id': ['MGI:1', 'MGI:3', 'MGI:4'], 'symbol': ['Mbp', 'Pakap', 'Other'], 'name': ['', '', '']})
    withdrawn = pd.DataFrame({'old_symbol': ['Akap2', 'Split', 'Split'], 'mgi_id': ['MGI:3', 'MGI:3', 'MGI:4']})
    table = map_genes(['Mbp', 'Gm1', 'Akap2', 'Split', 'Nope'], current, pc, withdrawn).set_index('source_gene')
    assert table.loc['Mbp', 'symbol'] == 'Mbp' and table.loc['Akap2', 'symbol'] == 'Pakap'
    assert table.loc['Gm1', 'status'] == 'not_protein_coding'
    assert table.loc['Split', 'status'] == 'ambiguous_withdrawn_symbol'
    assert table.loc['Nope', 'status'] == 'unmapped'


def test_section_counts_select_rows_and_sum_genes_that_share_a_symbol():
    # Three cells x four source genes, with column indices unsorted within rows as in the source; genes 1 and 2 now name one gene.
    indptr = np.array([0, 2, 4, 6])
    indices = np.array([2, 1, 3, 0, 3, 2], np.int32)
    data = np.array([1., 3, 5, 2, 1, 4])
    rows, cols, values = section_counts(indptr, indices, data, np.array([2, 0]), np.array([0, 1, 1, 2]), 3)
    assert rows.tolist() == [0, 0, 1] and cols.tolist() == [1, 2, 1] and values.tolist() == [4, 1, 4]  # row = position in the selection
    assert values.dtype == np.int64
    with pytest.raises(ValueError, match='positive integers'):
        section_counts(indptr, indices, data*0.5, np.array([0]), np.array([0, 1, 1, 2]), 3)


def test_rank_rows_matches_rankdata_with_ties_zeros_and_empty_rows():
    rng = np.random.default_rng(3)
    g, n = 40, 25
    dense = rng.integers(0, 4, (n, g)) * (rng.random((n, g)) < .5)
    dense[0] = 0                    # no detected genes
    dense[1] = rng.integers(1, 3, g)  # no zeros
    rows, cols = np.nonzero(dense)
    shuffle = rng.permutation(len(rows))  # aggregation output is unordered
    indptr, gene_index, rank2, zero2 = rank_rows(rows[shuffle], cols[shuffle], dense[rows, cols][shuffle], n, g)
    for row in range(n):
        decoded = np.full(g, zero2[row]/2)
        a, b = indptr[row], indptr[row+1]
        assert (np.diff(gene_index[a:b].astype(int)) > 0).all()
        decoded[gene_index[a:b]] = rank2[a:b]/2
        np.testing.assert_array_equal(decoded, rankdata(dense[row], method='average'))
    assert zero2[1] == 0 and zero2[0] == g+1


def test_pseudobulk_sums_counts_per_present_label():
    rng = np.random.default_rng(5)
    dense = rng.poisson(1.5, (12, 7))
    codes = np.array([4, 4, 1, 4, 9, 1, 1, 4, 9, 9, 4, 1])  # global label indices; 0, 2, 3 ... absent
    rows, cols = np.nonzero(dense)
    present, n_cells, sums = pseudobulk(rows, cols, dense[rows, cols], codes, 7)
    assert present.tolist() == [1, 4, 9] and n_cells.tolist() == [4, 5, 3] and sums.dtype == np.uint32
    for i, label in enumerate(present):
        np.testing.assert_array_equal(sums[i], dense[codes == label].sum(axis=0))


def test_separate_colors_splits_identical_colours_only():
    colors = separate_colors({'a': '#7A944D', 'b': '#7C9650', 'c': '#5858E0'})
    assert colors['c'] == '#5858E0'
    assert colors['a'] != colors['b']


def test_groups_and_classes_follow_the_paper():
    assert len(GROUPS) == 41 and len(set(GROUPS.values())) == 23 and len(set(CLASSES.values())) == 13
    members = lambda group: sorted(t for t, g in GROUPS.items() if g == group)
    assert members('Cortex IT excitatory neurons') == ['Cortex L2/3 IT', 'Cortex L2/3 RSP', 'Cortex L4 RSP', 'Cortex L4/5 IT',
                                                       'Cortex L5 IT', 'Cortex L6 IT']
    assert members('Cortex non-IT excitatory neurons') == ['Cortex 6b', 'Cortex L5 ET', 'Cortex L5 NP', 'Cortex L6 CT']
    assert members('Striatum neurons') == ['Striatum_D1', 'Striatum_D2', 'Striatum_LSr', 'Striatum_OT']
    assert sum(c == 'Neurons' for c in CLASSES.values()) == 11  # the paper's eleven neuron groups
    assert CLASSES['Other'] == 'Other' and CLASSES['Oligodendrocytes'] == 'Oligodendrocytes'


def rttamp_cells():
    types = ['Cortex L2/3 IT']*3 + ['Cortex L6 IT'] + ['Astrocytes']*2 + ['NFOLs', 'Other']
    cells = pd.DataFrame({'cell_type': types})
    cells['cell_group'] = cells.cell_type.map(GROUPS)
    cells['cell_class'] = cells.cell_group.map(CLASSES)
    return cells


def test_cell_colours_step_through_family_hues_and_stand_out_on_black():
    levels = cell_levels(rttamp_cells())
    types = levels['cell_type']['labels']
    labs = np.array([oklab(types[t]) for t in ['Cortex L2/3 IT', 'Cortex L6 IT']])
    hue = np.degrees(np.arctan2(labs[:, 2], labs[:, 1]))
    assert np.ptp(hue) < 3 and labs[0, 0] < labs[1, 0]  # one hue, lighter in name order
    # Steps span the whole family, so a type's colour does not depend on which types are present.
    assert labs[0, 0] == pytest.approx(TYPE_BAND[0], abs=0.01) and labs[1, 0] == pytest.approx(TYPE_BAND[1], abs=0.01)
    assert types['Other'] == UNKNOWN
    # A group or class takes the colour of its most numerous type.
    assert levels['cell_group']['labels']['Cortex IT excitatory neurons'] == types['Cortex L2/3 IT']
    assert levels['cell_class']['labels'] == {'Astrocytes': types['Astrocytes'], 'NFOLs': types['NFOLs'],
                                              'Neurons': types['Cortex L2/3 IT'], 'Other': UNKNOWN}
    assert levels['cell_type']['descriptions'] == {'NFOLs': NAMES['NFOLs']}  # only the labels present
    assert all(contrast(c) >= MIN_CONTRAST for level in levels.values() for c in level['labels'].values())
    with pytest.raises(ValueError, match='without a group'):
        cell_levels(rttamp_cells().assign(cell_type='Purkinje'))


def test_categorical_reads_string_and_integer_categories(tmp_path):
    with h5py.File(tmp_path/'obs.h5', 'w') as h5:
        h5['obs/cell_label/categories'] = np.array(['Astrocytes', 'OPCs'], dtype=h5py.string_dtype())
        h5['obs/cell_label/codes'] = np.array([1, 0, 1], np.int8)
        h5['obs/sample/categories'] = np.array([1, 3])
        h5['obs/sample/codes'] = np.array([0, 1, 1], np.int8)
        h5['var/codebook/categories'] = np.array(['WT_Codebook_1'], dtype=h5py.string_dtype())
        h5['var/codebook/codes'] = np.array([-1, 0], np.int8)
    with h5py.File(tmp_path/'obs.h5', 'r') as h5:
        assert categorical(h5, 'cell_label').tolist() == ['OPCs', 'Astrocytes', 'OPCs']
        assert categorical(h5, 'sample').tolist() == ['1', '3', '3']
        with pytest.raises(ValueError, match='missing'):
            categorical(h5, 'codebook', 'var')
