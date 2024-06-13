mkdir -p figures
cp outputs/scenario_1/plots/results_table.pdf figures/table_1.pdf
cp outputs/scenario_2/plots/results_table.pdf figures/table_2.pdf
cp outputs/scenario_3/plots/results_table.pdf figures/table_3.pdf
cp outputs/scenario_4/plots/full_panel.pdf figures/figure_2.pdf
cp outputs/scenario_5/plots/results_table.pdf figures/table_4.pdf
cp outputs/scenario_1/plots/full_panel.pdf figures/sup_figure_B.pdf
cp outputs/scenario_2/plots/full_panel.pdf figures/sup_figure_C.pdf
cp outputs/scenario_3/plots/full_panel.pdf figures/sup_figure_D.pdf
cp outputs/scenario_5/plots/full_panel.pdf figures/sup_figure_E.pdf
ipython plot/isolated.py
ipython plot/runtimes.py
ipython plot/boxp.py
