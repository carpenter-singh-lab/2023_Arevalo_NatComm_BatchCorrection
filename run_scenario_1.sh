#!/bin/bash
python run_correction.py inputs/conf/scenario_1/raw.json sphering outputs/scenario_1
python run_evaluation.py inputs/conf/scenario_1/raw.json sphering 5 outputs/scenario_1
python eval_map.py inputs/conf/scenario_1/raw.json sphering outputs/scenario_1

python run_correction.py inputs/conf/scenario_1/mad.json sphering outputs/scenario_1
python run_evaluation.py inputs/conf/scenario_1/mad.json sphering 5 outputs/scenario_1
python eval_map.py inputs/conf/scenario_1/mad.json sphering outputs/scenario_1

python run_correction.py inputs/conf/scenario_1/sphering.json sphering outputs/scenario_1
python run_evaluation.py inputs/conf/scenario_1/sphering.json sphering 5 outputs/scenario_1
python eval_map.py inputs/conf/scenario_1/sphering.json sphering outputs/scenario_1

# Running all methods for mad+sphering
python run_correction.py inputs/conf/scenario_1/mad_sphering.json all outputs/scenario_1
python run_evaluation.py inputs/conf/scenario_1/mad_sphering.json all 5 outputs/scenario_1
python eval_map.py inputs/conf/scenario_1/mad_sphering.json all outputs/scenario_1
