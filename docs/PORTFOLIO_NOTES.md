# RobuProt — Portfolio notes

**Context:** Undergraduate student coursework for BIOE 4266, Analytical Techniques in Bioengineering, Marmara University (spring 2025–2026). This was an exploratory educational research project, not a commissioned or clinically validated product.

## What problem did we tackle?

In LC-MS/MS proteomics, imputation and correction decisions can alter which proteins appear different between sample groups. RobuProt is a course prototype to compare those results under six preprocessing scenarios rather than relying on a single analysis.

## What was implemented?

The Python workflow downloads a public reference dataset, prepares protein intensities, visualizes missingness and PCA, calculates exploratory group comparisons in six scenarios, ranks protein-level signals with a **custom heuristic score**, and generates summary CSV/TXT outputs plus a formatted Word document.

## Authorship and contributions

**Research/report:** Pınar Nur Özdemir and Öykü Akdeniz, Marmara University, 2026.

**Programming:** Pınar Nur Özdemir led the software implementation, using AI-assisted coding tools during development. This description is based on the contributor's account; no per-commit contribution record is available.

## How should readers interpret the results?

The report records original coursework findings, not independently validated clinical biomarkers. The batch assignments in the current code are assumptions, not verified source metadata. The current statistical selection uses raw p-values rather than corrected FDR cutoffs, and the custom imputation differs from a canonical MinProb/QRILC implementation. Low-scoring signals should not be treated as *proven* false positives.

## Potential improvements

Verify biological groups and real batch metadata; implement validated imputation; account for multiple testing; add dataset-version provenance, tests, and validated out-of-sample evaluation; replace hard-coded document claims with computed results.

## Local execution check

The student developer reported a successful complete pipeline run on macOS with Python 3.10, including CSV/TXT exports, four figures and a Word report. The reported reference-scenario counts (260 candidates split into 13/215/32) matched the coursework report. SciPy numerical precision warnings were observed for some comparisons. The match verifies reproducibility of the current code path on that computer, not validity of the statistical assumptions or clinical usefulness. The original analytical code has **not** been modified in this portfolio package.
