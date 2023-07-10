import re
import itertools
from . import embeddings
from . import param_search
from loader import MAPPER
import pandas as pd
import numpy as np
from matplotlib.lines import Line2D

import plotly.express as px

rgb_reg = re.compile(r'rgb\((\d+), (\d+), (\d+)\)')
def rgb_to_hex(code) -> str:
    r, g, b = rgb_reg.match(code).groups()
    return '#{:02x}{:02x}{:02x}'.format(int(r), int(g), int(b))

BATCH_CMAP = {
    f'{i:02d}': color
    for i, color in enumerate(px.colors.qualitative.Light24, 1)
}
SOURCE_CMAP = {
    f'{chr(i)}': color
    for i, color in enumerate(px.colors.qualitative.Dark24, ord('A'))
}

METHODS = ['Harmony', 'Scanorama', 'Combat', 'MNN', 'DESC', 'scVI', 'Raw', 'MAD+Sphering', 'MAD', 'Sphering']
METHOD_CMAP = dict(zip(METHODS, px.colors.qualitative.D3))
METHOD_CMAP['Baseline'] = METHOD_CMAP['MAD+Sphering']

METHOD_SMAP = dict(zip(METHODS, Line2D.filled_markers[1:]))
METHOD_SMAP['Baseline'] = METHOD_SMAP['MAD+Sphering']

cmpds = pd.read_csv('inputs/metadata/compound.csv.gz')['Metadata_JCP2022']
cmpds = cmpds.apply(lambda x: MAPPER.get(x, x)).sort_values(key=lambda x: x.str.lower())
poscon_ix = ~(cmpds.str.startswith('JCP') | cmpds.isin(['DMSO', 'UNTREATED']))
POSCONS = cmpds[poscon_ix].tolist()
POSCON_MAP = dict(zip(POSCONS, px.colors.qualitative.Alphabet))
COMPOUND_COLORS = np.roll(px.colors.qualitative.Alphabet, -len(POSCONS)).tolist()


MICRO_CMAP = dict(zip(
   ['CV8000', 'Opera Phenix', 'ImageXpress Micro Confocal'],
   map(rgb_to_hex, px.colors.qualitative.Vivid)
))
