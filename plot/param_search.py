import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots


def plot_grid_exploration(scores: pd.DataFrame,
                          redlist: list[str] | None = None):
    if redlist:
        scores = scores[~scores['metric'].isin(redlist)]
    mean_scores = (scores.groupby(
        ['epsilon_mad', 'sphering_lambda', 'dimension', 'config_id'
         ], ).score.mean().reset_index())

    fig = make_subplots(rows=1,
                        cols=2,
                        horizontal_spacing=0.15,
                        subplot_titles=('Bio', 'Batch'))
    colorscale = 'RdBu'
    z_bio = mean_scores.query('dimension=="bio"').pivot(
        index='sphering_lambda', columns='epsilon_mad', values='score')
    z_batch = mean_scores.query('dimension=="batch"').pivot(
        index='sphering_lambda', columns='epsilon_mad', values='score')
    bio_best = mean_scores.query('dimension=="bio"').sort_values(
        'score').iloc[-1]
    batch_best = mean_scores.query('dimension=="batch"').sort_values(
        'score').iloc[-1]

    fig.add_traces([
        go.Contour(z=z_bio,
                   x=z_bio.columns,
                   y=z_bio.index,
                   colorscale=colorscale,
                   colorbar_x=0.42,
                   connectgaps=True),
        go.Scatter(x=[bio_best.epsilon_mad],
                   y=[bio_best.sphering_lambda],
                   marker_symbol='cross',
                   marker_line_color="midnightblue",
                   marker_color="lightskyblue",
                   marker_line_width=1,
                   marker_size=10,
                   showlegend=True,
                   name=f'bio_best_{bio_best.config_id}')
    ], 1, 1)
    fig.add_traces([
        go.Contour(z=z_batch,
                   x=z_bio.columns,
                   y=z_bio.index,
                   colorscale=colorscale,
                   connectgaps=True),
        go.Scatter(y=[batch_best.sphering_lambda],
                   x=[batch_best.epsilon_mad],
                   marker_symbol='x',
                   marker_line_color="midnightblue",
                   marker_color="lightskyblue",
                   marker_line_width=1,
                   marker_size=10,
                   showlegend=True,
                   name=f'batch_best_{batch_best.config_id}')
    ], 1, 2)

    fig['layout']['yaxis1'].update(title='sphering_lambda')
    fig['layout']['xaxis1'].update(title='epsilon_mad')
    fig['layout']['xaxis2'].update(title='epsilon_mad')
    fig.update_layout(legend=dict(
        yanchor="bottom",
        y=1.03,
        xanchor="center",
    ))
    fig.update_xaxes(type="log")
    fig.update_yaxes(type="log")
    return fig


def plot_agg_exploration(scores: pd.DataFrame, metric=None,
                         redlist: list[str] | None = None):
    if redlist:
        scores = scores[~scores['metric'].isin(redlist)]
    if metric:
        scores = scores[scores['metric'] == metric]
    colorscale = 'RdBu'
    mean_scores = (scores.groupby(
        ['epsilon_mad', 'sphering_lambda', 'config_id'
         ], ).score.mean().reset_index())
    z = mean_scores.pivot(index=['sphering_lambda'],
                          columns='epsilon_mad',
                          values='score')
    overall_best = mean_scores.sort_values('score').iloc[-1]
    fig = go.Figure()
    fig.add_traces([
        go.Contour(z=z,
                   x=z.columns,
                   y=z.index,
                   colorscale=colorscale,
                   connectgaps=True),
        go.Scatter(x=[overall_best.epsilon_mad],
                   y=[overall_best.sphering_lambda],
                   marker_symbol='hexagram',
                   marker_line_color="midnightblue",
                   marker_color="lightskyblue",
                   marker_line_width=1,
                   marker_size=10,
                   showlegend=True,
                   name=overall_best.config_id)
    ])
    fig.update_layout(legend=dict(
        yanchor="bottom",
        y=1.03,
        xanchor="center",
    ))

    fig.update_xaxes(type="log").update_yaxes(type="log")

    fig['layout']['yaxis'].update(title='sphering_lambda')
    fig['layout']['xaxis'].update(title='epsilon_mad')
    return fig


def rank_scores(scores: pd.DataFrame, redlist: list[str] | None = None):
    if redlist:
        scores = scores[~scores['metric'].isin(redlist)]
    mean_scores = (scores.groupby(
        ['epsilon_mad', 'sphering_lambda', 'dimension', 'config_id'
         ], ).score.mean().reset_index())
    rank = mean_scores.pivot(index='config_id',
                             columns='dimension',
                             values='score')
    rank['overall'] = (rank.batch + rank.bio) / 2

    rank = rank.sort_values('bio', ascending=False)
    rank['bio_rank'] = range(len(rank))

    rank = rank.sort_values('batch', ascending=False)
    rank['batch_rank'] = range(len(rank))

    rank = rank.sort_values('overall', ascending=False)
    rank['overall_rank'] = range(len(rank))
    rank = rank.sort_values('overall_rank')
    return rank
