#/bin/bash
python run_correction.py inputs/conf/scenario_2.json all outputs/scenario_2
python run_evaluation.py inputs/conf/scenario_2.json all 5 outputs/scenario_2
python eval_map.py inputs/conf/scenario_2.json all outputs/scenario_2
