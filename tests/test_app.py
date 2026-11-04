from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest
from spatial_app.enrichment import MIN_CELLS
from spatial_app.bins import MIN_BIN_CELLS
from spatial_app.layout import dot_size

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/rttamp-genes'
BINS = ROOT/'data/rttamp-genes-bins178um'
pytestmark = pytest.mark.skipif(not (DATA/'manifest.json').exists(), reason='prepared RT&T-AMP-MERFISH dataset not available')


@pytest.fixture
def manifest():
    return json.loads((DATA/'manifest.json').read_text())


def finish(app, manifest):
    """Advance the fragment timer until every batch of every section is scored."""
    batches = {b for s in manifest['slices'] for b in s['batches']}
    for _ in range(len(batches)+5):  # one batch per rerun
        if batches <= app.session_state['done']:
            return
        app.run()
        assert not app.exception
    raise AssertionError('scoring did not finish')


def test_seven_hemisections_score_and_fill_both_tables(manifest):
    app = AppTest.from_file(str(ROOT/'app.py'), default_timeout=180).run()
    assert not app.exception
    assert manifest['n_slices'] == 7 and len(app.session_state['scores']) == manifest['n_cells']
    header = next(m.value for m in app.markdown if 'hemisections' in m.value)
    assert header == '**Cohen-genes · 7 coronal hemisections, anterior to posterior**'
    panels = app.session_state['cell_plot_transport']['geometry']['panels']
    assert [p['label'] for p in panels] == [s['id'] for s in manifest['slices']] == ['1', '2', '4', '7', '5', '6', '3']
    planes = [s['ccf_x_um'] for s in manifest['slices']]  # anterior to posterior by the CCF plane each section matches
    assert planes == sorted(planes) and min(s['plane_match'] for s in manifest['slices']) > 0.9
    geometry = app.session_state['cell_plot_transport']['geometry']  # cells are drawn 4% wider than their median spacing
    assert geometry['dot_scale'] == 1.04 and geometry['dot'] == pytest.approx(1.04*dot_size(pd.read_parquet(DATA/'cells.parquet').sort_values('draw_order', ignore_index=True)))
    assert app.session_state['job']['query'].universe_size == manifest['n_genes'] == 22919
    app.run()
    finish(app, manifest)
    assert np.isfinite(app.session_state['scores']).all()
    cell_table, region_table = app.dataframe[0].value, app.dataframe[1].value
    cells = pd.read_parquet(DATA/'cells.parquet')
    for table, column in [(cell_table, 'cell_type'), (region_table, 'region_structure')]:
        counts = cells[column].value_counts()  # all cells are profiled
        assert (table['Cells'] >= MIN_CELLS).all()
        assert table.set_index(table.columns[0])['Cells'].to_dict() == counts[counts >= MIN_CELLS].to_dict()
    assert cell_table.columns[:2].tolist() == ['Cell type', 'Cell group']
    assert region_table.columns[:4].tolist() == ['Acronym', 'Full name', 'Cells', '# sections'] and '# sections' not in cell_table
    sections = cells.groupby('region_structure', observed=True).slice.nunique()  # sections contributing cells
    assert region_table.set_index('Acronym')['# sections'].to_dict() == sections.reindex(region_table['Acronym']).to_dict()
    # The default set (Tac2, Chat, Slc5a7, Chrnb3, Pou4f1, Gpr151) marks the medial habenula.
    assert app.session_state['job']['query'].target_symbols == ('CHAT', 'CHRNB3', 'GPR151', 'POU4F1', 'SLC5A7', 'TAC2')
    matching = next(e for e in app.expander if e.label.startswith('Gene matching'))
    assert matching.label == 'Gene matching · 6 targets used in a 22,919-gene universe' and not matching.warning
    assert cell_table.iloc[0, :2].tolist() == ['Cholinergic_TH', 'Cholinergic neurons'] and region_table.iloc[0, 0] == 'MH'
    assert region_table.iloc[0]['# sections'] == 2  # the habenula is in samples 4 and 7, registered to the same plane


def test_class_group_and_substructure_levels_and_colour_modes(manifest):
    app = AppTest.from_file(str(ROOT/'app.py'), default_timeout=180).run()
    app.selectbox(key='cell_level').set_value('cell_group').run()
    app.selectbox(key='region_level').set_value('region_substructure').run()
    assert not app.exception
    group_table = app.dataframe[0].value
    assert group_table.columns[:2].tolist() == ['Cell group', 'Cell class']
    classes = pd.read_parquet(DATA/'cells.parquet').drop_duplicates('cell_group').set_index('cell_group').cell_class
    assert dict(zip(group_table['Cell group'], group_table['Cell class'])) == classes.reindex(group_table['Cell group']).to_dict()
    assert len(group_table) == 23 and (group_table['Cell class'] == 'Neurons').sum() == 11
    assert app.dataframe[1].value.columns[0] == 'Acronym' and app.subheader[1].value == 'Region specificity · brain substructure'
    app.selectbox(key='cell_level').set_value('cell_class').run()
    assert app.dataframe[0].value.columns[1] == 'Cells'  # the top level has no parent column
    for mode in ['Cell type', 'Region', 'Cell type AUROC', 'Region AUROC']:
        app.radio(key='colour_mode').set_value(mode).run()
        assert not app.exception


def test_gene_matching_reports_genes_off_the_panel_and_human_input_maps_to_mouse(manifest):
    app = AppTest.from_file(str(ROOT/'app.py'), default_timeout=180).run()
    app.text_area(key='target_input').set_value('Glp1r Zic1 Snap25').run()
    app.text_area(key='background_input').set_value('Glp1r Gfap Aqp4 Mbp').run()
    assert not app.exception
    matching = next(e for e in app.expander if e.label.startswith('Gene matching'))
    assert matching.label == 'Gene matching · 1 target used in a 3-gene universe · 2 target symbols excluded'
    assert [w.value for w in matching.warning] == ['Targets not in dataset: SNAP25', 'Background genes not in dataset: MBP',
                                                   'Targets outside background: ZIC1']  # Snap25 and Mbp are too abundant for the panel
    app.text_area(key='background_input').set_value('').run()
    app.selectbox(key='input_species').set_value('Human').run()
    app.text_area(key='target_input').set_value('GLP1R ZIC1 TRPC4').run()
    assert not app.exception and app.session_state['job']['query'].target_symbols == ('GLP1R', 'TRPC4', 'ZIC1')


@pytest.mark.skipif(not (BINS/'manifest.json').exists(), reason='binned RT&T-AMP-MERFISH dataset not available')
def test_bins_are_one_more_colour_mode_and_rescored_for_a_new_gene_set():
    bins = json.loads((BINS/'manifest.json').read_text())
    mode = f"{bins['bin_um']:g}µm bin AUROC"
    app = AppTest.from_file(str(ROOT/'app.py'), default_timeout=180).run()
    assert app.radio(key='colour_mode').options == ['Gene-set AUROC', 'Cell type AUROC', 'Region AUROC', mode, 'Cell type', 'Region']
    app.radio(key='colour_mode').set_value(mode).run()
    assert not app.exception
    ((key, (values, reference)),) = app.session_state['bin_scores'].items()
    assert key[1] == bins['created_utc'] and reference is None  # keyed by the bins build too, and scored without error
    n_cells = pd.read_parquet(BINS/'bins.parquet').sort_values('draw_order').n_cells.to_numpy()
    assert len(values) == bins['n_bins'] and np.array_equal(np.isfinite(values), n_cells >= MIN_BIN_CELLS)
    assert ((values[np.isfinite(values)] >= 0) & (values[np.isfinite(values)] <= 1)).all()
    assert any(c.value.startswith(f"Each square is a {bins['bin_um']:g} µm bin") for c in app.caption)
    geometry = app.session_state['cell_plot_transport']['geometry']
    assert (geometry['items'], geometry['marker'], geometry['um_per_unit'], geometry['noun']) == ('bins', 'square', 1.0, 'sections')
    assert geometry['dot_scale'] == 1.0  # bin squares keep their side
    app.text_area(key='target_input').set_value('Gfap Aqp4 Aldh1l1').run()
    assert not app.exception and app.radio(key='colour_mode').value == 'Gene-set AUROC'  # a new gene set starts on the cells
    assert any(c.value.startswith('Each dot is one cell') for c in app.caption)
    app.radio(key='colour_mode').set_value(mode).run()
    ((_, (astrocyte, _)),) = app.session_state['bin_scores'].items()  # scores for the current query only
    assert not app.exception and not np.allclose(astrocyte, values, equal_nan=True)


@pytest.mark.skipif(not (BINS/'manifest.json').exists(), reason='binned RT&T-AMP-MERFISH dataset not available')
def test_a_failed_bin_scoring_shows_only_a_reference_and_is_not_retried(tmp_path, monkeypatch):
    bins = json.loads((BINS/'manifest.json').read_text())
    (tmp_path/'ranks').mkdir()
    for name in ['manifest.json', 'bins.parquet', *(s['ranks'] for s in bins['slices'][:-1])]:  # the last section's ranks are missing
        (tmp_path/name).symlink_to(BINS/name)
    monkeypatch.setenv('RTTAMP_APP_BINS', str(tmp_path))
    app = AppTest.from_file(str(ROOT/'app.py'), default_timeout=180).run()
    app.radio(key='colour_mode').set_value(app.radio(key='colour_mode').options[3]).run()  # the bin mode
    assert not app.exception
    ((_, (values, reference)),) = app.session_state['bin_scores'].items()
    assert reference and np.isnan(values).all()
    assert [e.value for e in app.error] == [f'Unable to score the bins (reference {reference}).']  # no paths or exception text
    app.run()  # the next fragment tick shows the same failure instead of scoring again
    assert [r for _, (_, r) in app.session_state['bin_scores'].items()] == [reference]
