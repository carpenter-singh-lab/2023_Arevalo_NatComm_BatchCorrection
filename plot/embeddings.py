'''Plot all figures'''
import logging
import warnings

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from sklearn.preprocessing import minmax_scale

logger = logging.getLogger(__name__)


def eigen_curve(spherer_path: str) -> go.Figure:
    '''Plot eigen curve'''
    spherer = np.load(spherer_path, allow_pickle=True)
    spherer = spherer['spherer'].item()
    lambda_loc = (spherer.eigenvals > spherer.regularization).sum()
    title = f'Eigenvalues of spherized data, {spherer.regularization:0.2e}'
    data = pd.DataFrame({
        'y': list(spherer.eigenvals),
        'x': range(len(spherer.eigenvals))
    })
    fig = px.line(data, x='x', y='y', log_y=True, title=title)
    return fig.add_vline(x=lambda_loc,
                         line_width=3,
                         line_dash="dash",
                         line_color="blue")


def _jitter(trace, win_size=0.02):
    range_x = max(trace.x) - min(trace.x)
    factor = range_x * win_size
    trace.x += np.random.uniform(-factor, factor, [len(trace.x)])

    range_y = max(trace.y) - min(trace.y)
    factor = range_y * win_size
    trace.y += np.random.uniform(-factor, factor, [len(trace.y)])


def load_embeddings(embd_paths: list[str], anon=True):
    embds = []
    for path in embd_paths:
        embd = pd.read_parquet(path)
        _jitter(embd)
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            # hack to make plotly happy in ranges
            embd[['x', 'y']] = minmax_scale(embd[['x', 'y']])
        embds.append(embd)
        embd['path'] = path
    embds = pd.concat(embds)

    if anon:
        # Make batches and sources annonymous
        renamer = {
            batch: f'{i:02d}'
            for i, batch in enumerate(embds['Metadata_Batch'].unique(), 1)
        }
        embds['Metadata_Batch'] = embds['Metadata_Batch'].apply(renamer.get)

        renamer = {
                source: f'{int(source.split("_")[-1]):02d}' for source in embds['Metadata_Source'].unique()
        }
        embds['Metadata_Source'] = embds['Metadata_Source'].apply(renamer.get)

    embds = embds.rename(
        columns={
            'Metadata_Batch': 'Batch',
            'Metadata_Source': 'Source',
            'Metadata_JCP2022': 'Compound',
            'Metadata_Row': 'Row',
            'Metadata_Column': 'Column'
        })
    return embds
