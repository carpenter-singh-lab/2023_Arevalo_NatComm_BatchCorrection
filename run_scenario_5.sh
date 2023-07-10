#/bin/bash
python run_correction.py inputs/conf/scenario_5.json all outputs/scenario_5
python run_evaluation.py inputs/conf/scenario_5.json all 5 outputs/scenario_5
python eval_map.py inputs/conf/scenario_5.json all outputs/scenario_5
