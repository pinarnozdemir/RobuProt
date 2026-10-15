import pandas as pd

# ── Tüm sonuçları yükle ────────────────────────────────────────────────────
scores  = pd.read_csv("results/robustness_scores.csv", index_col=0)
classic = pd.read_csv("results/deresult_MinProb_noBatch.csv", index_col=0)

scenarios = ['KNN_noBatch','KNN_Batch','MinProb_noBatch',
             'MinProb_Batch','Median_noBatch','Median_Batch']
all_res = {s: pd.read_csv(f"results/deresult_{s}.csv", index_col=0)
           for s in scenarios}

sig = classic[classic['significant']==True].copy()
# Keep summary consistent with step4.py: remove contaminants/reverse hits
sig = sig[~sig.index.astype(str).str.startswith(('CON__', 'REV__'))]
high = scores[scores['category']=='High confidence'].sort_values('total_score', ascending=False)
mid  = scores[scores['category']=='Moderate confidence'].sort_values('total_score', ascending=False)
risk = scores[scores['category']=='Likely artifact']

artifacts_in_classic = sig[sig.index.isin(risk.index)]

# ── TXT dosyasına yaz ─────────────────────────────────────────────────────
with open("results/ALL_RESULTS_SUMMARY.txt", "w") as f:

    f.write("="*65 + "\n")
    f.write("PROTEOMICS ROBUSTNESS PIPELINE — COMPLETE RESULTS SUMMARY\n")
    f.write("="*65 + "\n\n")

    f.write("DATASET\n")
    f.write("-"*40 + "\n")
    f.write("Source     : StatOmics/pda — Ghent University\n")
    f.write("Disease    : Squamous cell carcinoma (OR, n=9)\n")
    f.write("Control    : Healthy controls (PD, n=9)\n")
    f.write("Raw proteins: 3,898\n")
    f.write("After filtering (<=60% missing): 3,179\n\n")

    f.write("="*65 + "\n")
    f.write("DIFFERENTIAL ANALYSIS — 6 SCENARIOS\n")
    f.write("="*65 + "\n")
    f.write(f"{'Scenario':<25} {'Imputation':<18} {'Batch':<8} {'Sig. Proteins':>15}\n")
    f.write("-"*65 + "\n")
    scenario_info = [
        ("KNN_noBatch",     "KNN (k=5)",       "No",  ),
        ("KNN_Batch",       "KNN (k=5)",       "Yes", ),
        ("MinProb_noBatch", "MinProb/QRILC",   "No",  ),
        ("MinProb_Batch",   "MinProb/QRILC",   "Yes", ),
        ("Median_noBatch",  "Median",          "No",  ),
        ("Median_Batch",    "Median",          "Yes", ),
    ]
    for s, imp, batch in scenario_info:
        n = all_res[s]['significant'].sum()
        f.write(f"{s:<25} {imp:<18} {batch:<8} {n:>15}\n")
    f.write("\n")

    f.write("="*65 + "\n")
    f.write("ROBUSTNESS SCORE SUMMARY\n")
    f.write("="*65 + "\n")
    f.write(f"Total proteins analyzed  : {len(scores)}\n")
    f.write(f"High confidence (>=7/10) : {len(high)}\n")
    f.write(f"Moderate conf.  (4-6/10) : {len(mid)}\n")
    f.write(f"Likely artifact (<4/10)  : {len(risk)}\n\n")
    f.write(f"Classic pipeline total   : {len(sig)}\n")
    f.write(f"  -> High confidence     : {len(sig[sig.index.isin(high.index)])}\n")
    f.write(f"  -> Moderate confidence : {len(sig[sig.index.isin(mid.index)])}\n")
    f.write(f"  -> Likely artifact     : {len(artifacts_in_classic)} ({len(artifacts_in_classic)/len(sig)*100:.1f}%)\n\n")

    f.write("="*65 + "\n")
    f.write("HIGH CONFIDENCE BIOMARKER CANDIDATES (score >=7/10)\n")
    f.write("="*65 + "\n")
    f.write(f"{'Protein':<12} {'Sig/6':>6} {'Missing%':>10} {'Mean_log2FC':>13} {'Score':>7} {'Category':<20}\n")
    f.write("-"*65 + "\n")
    for prot, row in high.iterrows():
        f.write(f"{prot:<12} {int(row['n_sig_of_6']):>6} {row['missing_pct']:>10.1f} "
                f"{row['mean_log2FC']:>13.3f} {int(row['total_score']):>7} {row['category']:<20}\n")
    f.write("\n")

    f.write("="*65 + "\n")
    f.write("CLASSIC PIPELINE — TOP 20 SIGNIFICANT PROTEINS\n")
    f.write("="*65 + "\n")
    f.write(f"Total significant proteins: {len(sig)}\n\n")
    f.write(f"{'Protein':<12} {'log2FC':>8} {'p-value':>10} {'FDR':>10}\n")
    f.write("-"*45 + "\n")
    for prot, row in sig.sort_values('pvalue').head(20).iterrows():
        f.write(f"{prot:<12} {row['log2FC']:>8.4f} {row['pvalue']:>10.4f} {row['fdr']:>10.4f}\n")
    f.write("\n")

    f.write("="*65 + "\n")
    f.write("LIKELY ARTIFACTS FOUND IN CLASSIC PIPELINE\n")
    f.write("="*65 + "\n")
    f.write(f"{'Protein':<12} {'log2FC':>8} {'p-value':>10} {'Rob.Score':>11}\n")
    f.write("-"*45 + "\n")
    for prot, row in artifacts_in_classic.sort_values('pvalue').iterrows():
        score_val = scores.loc[prot, 'total_score'] if prot in scores.index else 'N/A'
        f.write(f"{prot:<12} {row['log2FC']:>8.4f} {row['pvalue']:>10.4f} {str(score_val):>11}\n")



# ── CSV summary table for report consistency ───────────────────────────────
classic_total = len(sig)
high_in_classic = len(sig[sig.index.isin(high.index)])
mid_in_classic = len(sig[sig.index.isin(mid.index)])
risk_in_classic = len(artifacts_in_classic)
robust_total = high_in_classic + mid_in_classic
artifact_pct = (risk_in_classic / classic_total * 100) if classic_total else 0

summary_table = pd.DataFrame({
    'Metric': [
        'Total significant proteins',
        'High confidence (>=7/10)',
        'Moderate confidence (4-6/10)',
        'Likely artifacts removed'
    ],
    'Classic Pipeline': [
        classic_total,
        'not assessed',
        'not assessed',
        '—'
    ],
    'Robustness-Scored Pipeline': [
        robust_total,
        high_in_classic,
        mid_in_classic,
        f'{risk_in_classic} ({artifact_pct:.1f}%)'
    ]
})
summary_table.to_csv('results/robustness_summary_table.csv', index=False)

print("ALL_RESULTS_SUMMARY.txt kaydedildi -> results/ALL_RESULTS_SUMMARY.txt")
print("\nİçerik:")
print("-"*40)
with open("results/ALL_RESULTS_SUMMARY.txt") as f:
    print(f.read())
