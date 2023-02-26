# Hutch-AMM
Hutch estimated sampling probability for AMM

The code contains:

**AMM Experiment code:**

HutchAMM_experiment.m

p_bar_modified_exp.m

Hutch_AMM_exp2.m: *Experiment with synthetic data*

realdata_experiment.m: *Experiment with practical data*

**Visualization:**

plot_visualization.m

**Algorithm:**

simple_hutchinson.m: *simple Hutch trace estimator*

hutchplusplus.m: *hutchplusplus trace estimator*

AMM_true.m: *Standard (Optimal) AMM estimator, needs predefined group partition as input* (Archived)

AMM_coarse_hutch.m: *Hutch AMM estimator, needs predefined group partition as input* (Archived)

Hutch_AMM.m: *Hutch AMM estimator, input block size and this function split columns into groups in natural sequence.* (Archived)

AMM_true_tracefun_ver2.m: *Standard (Optimal) AMM estimator, needs predefined group partition as input, fair comparison*

AMM_coarse_hutch_ver4.m: *Hutch AMM estimator, needs predefined group partition as input, fair comparison*

AMM_coarse_uni_ver2.m: *Uniform AMM estimator, fair comparison*
