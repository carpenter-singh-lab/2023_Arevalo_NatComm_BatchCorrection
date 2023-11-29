'''Plot all figures'''
import logging
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from scipy.spatial.distance import pdist
from sklearn.preprocessing import minmax_scale

from utils import preprocessing_model_name

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


def label_sort(embds: pd.DataFrame, label: str):
    medians = {}
    for label, embd in embds.groupby(label):
        error = np.median(pdist(embd[['x', 'y']].values))
        medians[label] = error
    return pd.Series(medians).sort_values().index.values


def load_embeddings(embd_paths: list[str], vis='mde', verbose=False, anon=True):
    embds = []
    for path in embd_paths:
        try:
            embd = pd.read_csv(path, dtype={'Metadata_Plate': str})
            # Capitalize first letter
            embd['Metadata_JCP2022'] = embd['Metadata_JCP2022'].str[0].str.upper() + embd['Metadata_JCP2022'].str[1:]
        except FileNotFoundError as exc:
            if verbose:
                logger.warning(f'Embeddings not ready: {exc.filename}')
            continue
        embd['model'] = preprocessing_model_name(loc.model, loc.config)
        _jitter(embd)
        embd[['x', 'y']] = minmax_scale(
            embd[['x', 'y']])  # hack to make plotly happy in ranges
        embds.append(embd)
    embds = pd.concat(embds)

    if anon:
        # Make batches and sources annonymous
        renamer = {
            batch: f'{i:02d}'
            for i, batch in enumerate(embds['Metadata_Batch'].unique(), 1)
        }
        embds['Metadata_Batch'] = embds['Metadata_Batch'].apply(renamer.get)

        renamer = {
            batch: f'{chr(i)}'
            for i, batch in enumerate(embds['Metadata_Source'].unique(),
                                      ord('A'))
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
