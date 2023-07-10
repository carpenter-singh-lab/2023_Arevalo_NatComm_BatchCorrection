#/bin/bash
python run_correction.py inputs/conf/scenario_4.json all outputs/scenario_4
python run_evaluation.py inputs/conf/scenario_4.json all 5 outputs/scenario_4
python eval_map.py inputs/conf/scenario_4.json all outputs/scenario_4
