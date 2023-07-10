#/bin/bash
python run_correction.py inputs/conf/scenario_3.json all outputs/scenario_3
python run_evaluation.py inputs/conf/scenario_3.json all 5 outputs/scenario_3
python eval_map.py inputs/conf/scenario_3.json all outputs/scenario_3
