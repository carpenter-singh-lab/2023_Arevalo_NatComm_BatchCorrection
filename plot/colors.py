'''Define palette for plots'''
import re

from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import plotly.express as px

from preprocessing.metadata import MAPPER

rgb_reg = re.compile(r'rgb\((\d+), (\d+), (\d+)\)')


def rgb_to_hex(code) -> str:
    r, g, b = rgb_reg.match(code).groups()
    return '#{:02x}{:02x}{:02x}'.format(int(r), int(g), int(b))


BATCH_CMAP = {
    f'{i:02d}': color
    for i, color in enumerate(px.colors.qualitative.Light24, 1)
}
# Source in order of appearance in the scenarios
SOURCES = [6, 2, 10, 3, 8]
SOURCE_CMAP = {
    f'{source:02d}': color
    for source, color in zip(SOURCES, px.colors.qualitative.Dark24)
}

METHOD_FMT = {
    'combat': 'Combat',
    'harmony': 'Harmony',
    'scvi': 'scVI',
    'scanorama': 'Scanorama',
    'mnn': 'MNN',
    'desc': 'DESC',
    'baseline': 'Baseline',
    'sphering': 'Sphering',
    'fastMNN': 'fastMNN',
    'seurat_rpca': 'Seurat RPCA',
    'seurat_cca': 'Seurat CCA'
}
METHODS = [
    'Baseline', 'Harmony', 'Scanorama', 'fastMNN', 'Seurat CCA', 'Seurat RPCA', 'MNN', 'Combat', 'DESC', 'scVI', 'Sphering'
]
METHOD_CMAP = dict(zip(METHODS, px.colors.qualitative.Light24))
METHOD_SMAP = dict(zip(METHODS, Line2D.filled_markers[1:]))

cmpds = pd.read_csv('inputs/metadata/compound.csv.gz')['Metadata_JCP2022']
cmpds = cmpds.apply(lambda x: MAPPER.get(x, x)).sort_values(
    key=lambda x: x.str.lower())
poscon_ix = ~(cmpds.str.startswith('JCP') | cmpds.isin(['DMSO', 'UNTREATED']))
POSCONS = cmpds[poscon_ix].tolist()
POSCON_MAP = dict(zip(POSCONS, px.colors.qualitative.Alphabet))
COMPOUND_COLORS = np.roll(px.colors.qualitative.Alphabet,
                          -len(POSCONS)).tolist()

MICRO_CMAP = dict(
    zip(['CV8000', 'Opera Phenix', 'ImageXpress Micro Confocal'],
        map(rgb_to_hex, px.colors.qualitative.Vivid)))

METRIC_FMT = {
    'silhouette_batch': 'Silhouette batch',
    'pcr_batch': 'PCR comparison',
    'pcr': 'PCR',
    'graph_conn': 'Graph connectivity',
    'kbet': 'KBET',
    'lisi_batch': 'LISI batch',
    'lisi_label': 'LISI label',
    'negcon_mean_map': 'mAP (controls)',
    'negcon_fraction_below_p': '%Retrieved (controls)',
    'negcon_fraction_below_corrected_p': '%Retrieved-corr (controls)',
    'nonrep_mean_map': 'mAP (nonrep)',
    'nonrep_fraction_below_p': '%Retrieved (nonrep)',
    'nonrep_fraction_below_corrected_p': '%Retrieved-corr (nonrep)',
    'nmi': 'Leiden NMI',
    'ari': 'Leiden ARI',
    'asw': 'Silhouette label',
    'il_f1': 'Isolated labels(F1)',
    'il_asw': 'Isolated labels(ASW)',
}
