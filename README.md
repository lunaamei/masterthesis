## Documentation

### Behavioural

1. **build_RDM.m**
   - **in:** `.mat` files from the multiple arrangements program
   - **out:** individual and group-averaged RDMs as `.mat` and `.png` files

2. **behavioural_correlations.m**
   - **in:** `.mat` files from the multiple arrangements program
   - Computes Pearson correlation of Gatys and Texform to natural RDMs
   - Runs one-sample t-tests on the correlations with Bonferroni correction (2 conditions × 3 tasks = 6)
   - Runs paired t-tests between conditions within task with Bonferroni correction (3 comparisons)

3. **seperability.m**
   - Computes separability index
   - Runs one-sample t-test (Bonferroni n=10: 2 conditions × 5 categories) and pairwise comparisons between conditions (Bonferroni n=10 category pairs)

---

### fMRI

1. **fmri_dprime_stats_latex.m**
   - Computes d-prime selectivity
   - Runs 3 (condition) × 5 (category) two-way repeated measures ANOVA
   - One-sample test against baseline (Bonferroni n=15)
   - Paired t-test: highest vs. second highest category (Bonferroni n=3 per ROI)
   - Paired condition tests (Bonferroni n=15: 5 categories × 3 comparisons per ROI)

2. **fmri_build_RDMs_FINAL.m**
   - Computes and saves trial-averaged RDMs for each subject

---

### Encoding Models

1. **dprime_may.m**
   - Computes d-prime selectivity
   - Runs 3 (condition) × 5 (category) two-way ANOVA
   - One-sample t-test with Bonferroni correction (n=15)
   - Paired t-test with Bonferroni correction (n=3 pairwise comparisons)

3. **RSA_may.m**
   - Correlates model RDMs to fMRI RDM
   - Runs repeated measures ANOVA with condition and category on the correlations
   - One-sample t-test with Bonferroni correction (n=6)
   - Paired t-test with Bonferroni correction (n=3 pairwise comparisons)

4. **rsa_2x3_anova.m**
   - Each correlation tested against baseline with a one-sample t-test (Bonferroni n=6: 2 models × 3 conditions per ROI)
   - Paired t-tests between models (Bonferroni n=3 per ROI)
   - Conditions compared within model (Bonferroni n=6: 3 conditions × 2 models per ROI)
