# Reinforcement learning modeling — ASD early learning project

Computational modeling code and data for the ASD early-learning project.
This folder is the **modeling component** of the project's main repository
(behavioral analyses and figure generation live elsewhere in that
repository).

Early-adolescent *Tsc2* and *Shank3B* mice (male and female) were trained on an
odor-based two-alternative forced choice (2AFC) task under deterministic
(100%) or probabilistic (80–90%) reward schedules. Trial-by-trial choice data
from Session 1 were fit with a family of reinforcement learning (RL) models,
including hybrid models that alternate between an engaged RL policy and a
disengaged side-biased policy via a hidden Markov model ([dynamic noise
estimation](https://doi.org/10.1016/j.jmp.2024.102842); Li et al., 2024,
*J. Math. Psych.*).

## Layout

```
analysis/                  analysis pipeline (MATLAB)
  cohorts.m                  cohort configuration (single source of truth)
  model_specs.m              candidate model specifications (bounds, priors)
  fit_models.m               MAP fitting of all candidate models (Session 1)
  compare_params.m           Het-vs-WT parameter comparisons
  validate_models.m          model validation by simulation
  recover_params.m           parameter recovery analysis
  p_engaged_tests.m          p(engaged) KL-divergence permutation tests
  kl_divergence_permutation_test.m, get_binned_p_engaged.m   (helpers)
  csv_to_mat.m               raw CSV -> processed .mat converter (data prep)
models/
  likelihood/                negative log-likelihood + latent-state functions
  generative/                generative (simulation) function of the winning model
data/
  raw/                       raw behavior CSVs (Sessions 1–2 where available)
  processed/                 processed Session-1 .mat tables (pipeline inputs)
results/                     archived model fits and p(engaged) trajectories
```

All paths in the code are resolved relative to this folder, so it works
wherever it is placed within the parent repository.

## Cohorts

All modeling uses Session 1 only.

| id            | line     | sex    | schedule      | n (Het/WT/KO) | models fitted |
|---------------|----------|--------|---------------|---------------|---------------|
| `TSC_M_det`   | *Tsc2*   | male   | deterministic | 15/13/–       | 6             |
| `TSC_F_det`   | *Tsc2*   | female | deterministic | 10/10/–       | 5             |
| `Shank_M_det` | *Shank3B*| male   | deterministic | 14/13/13      | 6             |
| `Shank_F_det` | *Shank3B*| female | deterministic | 11/11/9       | 5             |
| `TSC_M_prob`  | *Tsc2*   | male   | probabilistic | 10/10/–       | 8             |
| `Shank_M_prob`| *Shank3B*| male   | probabilistic | 11/11/–       | 5             |

Candidate models: `RL_epsilon`, `a0b3s`, `a0b3s_hybrid`, `a0b2s_hybrid`,
`a0b1s_hybrid` (winning model), plus `a01s_hybrid` (6-model cohorts) and
`Aab1s_hybrid`, `b1s_hybrid` (`TSC_M_prob`). See `analysis/model_specs.m`
for specifications.

## Requirements

MATLAB (developed on R2022a) with toolboxes: Optimization, Global
Optimization, Parallel Computing, and Statistics and Machine Learning.

## Usage

```matlab
addpath('analysis');
fit_models('TSC_M_det');       % fit all candidate models (slow: hours)
compare_params('TSC_M_det');   % Het-vs-WT parameter comparisons
validate_models('TSC_M_det');  % simulate winning model, validate curves
recover_params('TSC_M_det');   % parameter recovery
p_engaged_tests();             % KL-divergence permutation tests
```

Fitting and parameter recovery use `GlobalSearch` (sqp `fmincon`) with
`rng('default')` per subject-model, and are deterministic apart from the
(unseeded) draw of starting values; validation simulations are stochastic by
nature (set the `seed` option for reproducible output).

The archived results in `results/` allow regenerating all comparison and
validation figures without re-fitting.

## Author

Jing-Jing Li (jl3676@berkeley.edu), Collins Lab, UC Berkeley.
