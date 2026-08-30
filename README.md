# Genetic changes linked to two different syndromic forms of autism enhance reinforcement learning in adolescent male but not female mice
### Abstract: 
Autism Spectrum Disorder (ASD) is characterized by restricted and repetitive behaviors and social differences, both of which may manifest, in part, from underlying differences in corticostriatal circuits and reinforcement learning. Here, we investigated reinforcement learning in developing mice with mutations in either Tsc2 or Shank3, both high-confidence ASD risk genes associated with major syndromic forms of ASD. Using an odor-based two-alternative forced choice (2AFC) task, we tested early adolescent mice of both sexes and found male Tsc2 and Shank3B heterozygote (Het) mice showed enhanced learning performance compared to their wild type (WT) siblings. No gain of function was observed in females. Using a novel reinforcement learning (RL) based computational model to infer learning rate as well as policy-level task engagement and disengagement, we found that the gain of function in males was driven by an enhanced positive learning rate in both Tsc2 and Shank3B Het mice. The gain of function in Het males was absent when mice were trained with a probabilistic reward schedule. These findings in two ASD mouse models reveal a convergent learning phenotype that shows similar sensitivity to sex and environmental uncertainty. These data can inform our understanding of both strengths and challenges associated with autism, while providing further evidence that sex and experience of uncertainty modulate autism-related phenotypes.
<br/>
Manuscript status: under revision
<br />

### Below is a brief description of the folders and files contained within this directory:  
## Code
- **ITL**: use ITL_analysis.ipynb to set path to metadata, output folders, and csv behavior files. This will create folders with all files seen in "results" folder. Helpers: `analysis.py`, `tail_analysis.py`, `trimming.py`.
- **mat_to_csv**: code to convert .mat behavior file to .csv files
- **behavior**: `master_ASD_RLWM_JC.m` builds a session index (`makeDataIndex_ASD.m`) and converts each raw .mat behavior log into a `behaviorDF.csv` via `extract_behavior_df` / `get_RLWM_EventTimes`. Set `root_dir` and `strainNum` at the top of the driver. Add `code/mat_to_csv/` and `code/behavior/` (and its `lib/` subfolders) to the MATLAB path before running.
<br/>`ASD_odor_summary.m` - 3 hour block performance in AB; `ASD_learning_curve.m` - quartile learning curves in AB. Summary values from both scripts were exported and plotted in GraphPad Prism.
- **exploratory_not_reported**: exploratory analyses not included in manuscript, kept for reference
- **modeling**: all reinforcement learning modeling code and data. Includes a separate `README.md` for reference
## Data
- **csvs**: individual behavior `.csv` for each animal's session, separated into 6 folders by group (Tsc2/Shank3B × sex × deterministic/probabilistic).
- **metadata**: one row per session with age, date, genotype, sex, etc.
## Results
- **ITL**: results from ITL code with individual folder for each group. Folder contains trial level summaries, cutoffs from 93, 95, 98 percentile, and a shorter .csv with the main 95 percentile data that was used for supplemental figures in the manuscript. 

<br />

## Related Repositories
- https://github.com/Wilbrecht-Lab/RP_Box: Code for custom 2AFC box and task 

<br />

## Contributors
- **Juliana Chase**: first author, behavior pipeline, ITL data, repository construction
- **Jing-Jing Li**: reinforcement learning modeling 
- **Wan Chen Lin**: initial quartile learning curve MATLAB script
- **Hongli Wang**: 3-hour block MATLAB script
- **Lung-Hao Tai**: RP-Box code construction
