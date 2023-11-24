'''
Helper functions
'''
import hashlib
import json
from collections.abc import Iterable
from pathlib import Path

METRIC_MAP = {
    'ASW_label/batch': 'batch',
    'PCR_batch': 'batch',
    'graph_conn': 'batch',
    'kBET': 'batch',
    'iLISI': 'batch',
    'cLISI': 'bio',
    'percent_repl': 'bio',
    'null_th_repl': 'bio',
    'percent_matching': 'bio',
    'fraction_positive_p': 'bio',
    'fraction_positive_q': 'bio',
    'mean_average_precision': 'bio',
    'wasserstein_distance': 'bio',
    'null_matching': 'bio',
    'NMI_cluster/label': 'bio',
    'ARI_cluster/label': 'bio',
    'ASW_label': 'bio',
    'isolated_label_F1': 'bio',
    'isolated_label_silhouette': 'bio',
    # Below metrics are not used here, require cell cycle label
    'cell_cycle_conservation': 'bio',
    'hvg_overlap': 'bio',
    'trajectory': 'bio',
}

with open('feature_set.txt', 'r', encoding='utf8') as f_in:
    FEATURE_SET = f_in.read().splitlines()


def find_feat_cols(cols: Iterable[str]):
    '''Find column names for features'''
    feat_cols = [c for c in cols if not c.startswith('Meta')]
    return feat_cols


def find_meta_cols(cols: Iterable[str]):
    '''Find column names for metadata'''
    meta_cols = [c for c in cols if c.startswith('Meta')]
    return meta_cols


def hashname(config: dict):
    '''Create hashname from a configuration id'''
    utf8_encoded = json.dumps(config, sort_keys=True).encode('utf-8')
    data_md5 = hashlib.md5(utf8_encoded).hexdigest()
    return data_md5


def preprocessing_model_name(model, config: dict):
    '''Returns a readable name for this preprocessing'''
    if model == 'sphering':
        if config['mad_norm'] and config['sphering']:
            return 'MAD+Sphering'
        if config['mad_norm']:
            return 'MAD'
        if config['sphering']:
            return 'Sphering'
        else:
            return 'Raw'
    else:
        formatter = {'combat': 'Combat', 'harmony': 'Harmony', 'scvi': 'scVI',
                     'scanorama': 'Scanorama', 'mnn': 'MNN', 'desc': 'DESC',
        }
        return formatter[model]

class PathLocator():
    '''Define output locations given a configuration'''

    def __init__(self, config_path: str | Path, model: str,
                 min_replicates_to_eval: int, output_path: str | Path):
        '''Initialize config for this experiment '''
        # Config
        with open(config_path, encoding='utf8') as freader:
            self.config = json.load(freader)
        self.hashid = hashname(self.config)
        self.model = model
        self.min_replicates_to_eval = min_replicates_to_eval

        root_dir = Path(f'{output_path}/{self.hashid}/{model}')
        score_dir = root_dir / f'replicates_{min_replicates_to_eval}'
        score_dir.mkdir(parents=True, exist_ok=True)

        self.config_path = root_dir.parent / 'config.json'
        with self.config_path.open('w', encoding='utf8') as f_out:
            json.dump(self.config, f_out)

        self.profiles_path = root_dir.parent / 'profiles.h5ad'
        self.processed_plates_path = root_dir.parent / 'processed_plates.parquet'
        self.spherer_path = root_dir.parent / 'spherer.npz'

        self.corr_dist_matching_path = root_dir / 'corr_dist_matching.csv'
        self.null_dist_matching_path = root_dir / 'null_dist_matching.csv'
        self.average_precision_path = root_dir / 'average_precision.csv'
        self.map_target_path = root_dir / 'map_target.csv'

        self.score_path = score_dir / 'sc_scores.csv'
        self.corr_dist_path = score_dir / 'corr_dist.csv'
        self.null_dist_path = score_dir / 'null_dist.csv'

        self.corrected_path = root_dir / 'corrected.h5ad'
        self.mde_path = root_dir / 'embeddings_2d.csv'
        self.pca_path = root_dir / 'pca_2d.csv'
        self.umap_path = root_dir / 'umap_2d.csv'
