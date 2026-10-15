import os
import numpy as np
np.random.seed(42)
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.impute import KNNImputer

C_NS     = "#AAAAAA"
C_WEAK   = "#E8A838"
C_ROBUST = "#2196A8"
C_CLASS  = "#C44E52"
BG       = "#FAFAFA"
AX_BG    = "#F7F7F7"

df_norm = pd.read_csv("data/intensity_norm.csv", index_col=0)
df_bc   = pd.read_csv("data/intensity_batch_corrected.csv", index_col=0)
meta    = pd.read_csv("data/metadata.csv", index_col=0)

ctrl_cols = [c for c in meta[meta['group']=='Control'].index if c in df_norm.columns]
dis_cols  = [c for c in meta[meta['group']=='Disease'].index  if c in df_norm.columns]

def impute_knn(df, k=5):
    imp = KNNImputer(n_neighbors=k)
    arr = imp.fit_transform(df.T).T
    return pd.DataFrame(arr, index=df.index, columns=df.columns)

def impute_minprob(df):
    result = df.copy()
    for protein in df.index:
        row = df.loc[protein]
        obs = row.dropna()
        if len(obs) < 2:
            result.loc[protein] = row.fillna(df.min().min())
            continue
        mu    = obs.quantile(0.01)
        sigma = max(obs.std() * 0.3, 0.01)
        n_mis = row.isna().sum()
        if n_mis > 0:
            result.loc[protein, row.isna()] = np.random.normal(mu, sigma, n_mis)
    return result

def impute_median(df):
    return df.T.fillna(df.median(axis=1)).T

print("Running imputation (1-2 min)...")
matrices = {
    'KNN_noBatch'    : impute_knn(df_norm),
    'KNN_Batch'      : impute_knn(df_bc),
    'MinProb_noBatch': impute_minprob(df_norm),
    'MinProb_Batch'  : impute_minprob(df_bc),
    'Median_noBatch' : impute_median(df_norm),
    'Median_Batch'   : impute_median(df_bc),
}
print("  6 scenario matrices ready")

def diff_analysis(mat, ctrl_cols, dis_cols):
    results = []
    for protein in mat.index:
        cv = mat.loc[protein, ctrl_cols].dropna().values
        dv = mat.loc[protein, dis_cols].dropna().values
        if len(cv) < 3 or len(dv) < 3:
            results.append({'protein': protein, 'log2FC': np.nan,
                            'pvalue': 1.0, 'fdr': 1.0, 'significant': False})
            continue
        log2fc = dv.mean() - cv.mean()
        _, pval = stats.ttest_ind(dv, cv, equal_var=False)
        results.append({'protein': protein, 'log2FC': round(log2fc, 4),
                        'pvalue': pval, 'fdr': np.nan, 'significant': False})
    df_res = pd.DataFrame(results).set_index('protein')
    valid  = df_res['pvalue'] < 1.0
    if valid.sum() > 0:
        _, fdrs, _, _ = multipletests(df_res.loc[valid, 'pvalue'], method='fdr_bh')
        df_res.loc[valid, 'fdr'] = fdrs
    df_res['fdr'] = df_res['fdr'].fillna(1.0)
    df_res['significant'] = (df_res['pvalue'] < 0.05) & (df_res['log2FC'].abs() > 1.0)
    return df_res

print("\nRunning differential analysis...")
all_results = {}
for scenario, mat in matrices.items():
    res = diff_analysis(mat, ctrl_cols, dis_cols)
    all_results[scenario] = res
    res.to_csv(f"results/deresult_{scenario}.csv")
    print(f"  {scenario:25s}: {res['significant'].sum():3d} significant proteins")

classic     = all_results['MinProb_noBatch']
classic_sig = set(classic[classic['significant']].index)

if os.path.exists("results/robustness_scores.csv"):
    scores_temp = pd.read_csv("results/robustness_scores.csv", index_col=0)
    robust_sig  = set(scores_temp[scores_temp['total_score'] >= 7].index)
    print(f"\nRobust proteins (score >=7): {len(robust_sig)}")
else:
    freq       = pd.DataFrame({s: r['significant'] for s, r in all_results.items()}).sum(axis=1)
    robust_sig = set(freq[freq >= 4].index)
    print(f"\nRobust proteins (>=4/6): {len(robust_sig)}")

print(f"Classic pipeline: {len(classic_sig)} significant proteins")

# ── FIGURE 3 — Volcano ────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
fig.patch.set_facecolor(BG)

# Left — classic
ax = axes[0]
ax.set_facecolor(AX_BG)
res = classic.dropna(subset=['log2FC', 'pvalue'])
nlp = -np.log10(res['pvalue'].clip(lower=1e-15))
sig = res['significant']

ax.scatter(res.loc[~sig, 'log2FC'], nlp[~sig],
           c=C_NS, s=12, alpha=0.45, linewidths=0, zorder=2, label='Not significant')
ax.scatter(res.loc[sig, 'log2FC'], nlp[sig],
           c=C_CLASS, s=22, alpha=0.85, linewidths=0, zorder=3,
           label=f'Significant (n={sig.sum()})')
for prot, row in res[sig].nsmallest(6, 'pvalue').iterrows():
    ax.annotate(prot, (row['log2FC'], -np.log10(row['pvalue'])),
                fontsize=7, color='#222222', fontweight='500',
                xytext=(6, 5), textcoords='offset points',
                arrowprops=dict(arrowstyle='-', color='#999999', lw=0.6))
ax.axhline(-np.log10(0.05), color='#888', linewidth=1, linestyle='--', alpha=0.7)
ax.axvline(-1, color='#888', linewidth=1, linestyle='--', alpha=0.7)
ax.axvline( 1, color='#888', linewidth=1, linestyle='--', alpha=0.7)
ax.set_xlabel("log₂ Fold Change  (Disease / Control)", fontsize=11, color='#333')
ax.set_ylabel("−log₁₀(p-value)", fontsize=11, color='#333')
ax.set_title("Classic Pipeline\n(MinProb imputation, no batch correction)",
             fontsize=12, fontweight='bold', color='#222', pad=12)
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['left', 'bottom']].set_color('#CCC')
ax.tick_params(colors='#555')
ax.legend(fontsize=9, framealpha=0.9, edgecolor='#DDD', facecolor='white')

# Right — robustness
ax = axes[1]
ax.set_facecolor(AX_BG)
res2 = classic.dropna(subset=['log2FC', 'pvalue']).copy()
res2['robust'] = res2.index.isin(robust_sig)
nlp2 = -np.log10(res2['pvalue'].clip(lower=1e-15))

m_ns = ~res2['significant']
m_wk = res2['significant'] & ~res2['robust']
m_rb = res2['robust']

ax.scatter(res2.loc[m_ns, 'log2FC'], nlp2[m_ns],
           c=C_NS, s=12, alpha=0.4, linewidths=0, zorder=2, label='Not significant')
ax.scatter(res2.loc[m_wk, 'log2FC'], nlp2[m_wk],
           c=C_WEAK, s=18, alpha=0.75, linewidths=0, zorder=3,
           label=f'Weak signal (n={m_wk.sum()})  [score <7]')
ax.scatter(res2.loc[m_rb, 'log2FC'], nlp2[m_rb],
           c=C_ROBUST, s=55, alpha=0.95, linewidths=0.8,
           edgecolors='white', zorder=4,
           label=f'Robust biomarker (n={m_rb.sum()})  [score ≥7/10]')

# Her mavi nokta için etiket — x pozisyonuna göre sağa veya sola yönlendir
robust_df = res2[m_rb].copy()
robust_df['nlp'] = nlp2[m_rb]
robust_df = robust_df.sort_values('nlp', ascending=False)

# Manuel offset — her protein için ayrı ayar
offsets = {
    'P46108' : (25, 15),
    'P49006' : (25, 5),
    'P20962' : (25, -15),
    'Q9GZZ1' : (25, 5),
    'P63313' : (-75, 15),
    'P52566' : (-75, 5),
    'P16949' : (-75, -10),
    'P16401' : (-75, -20),
    'O43399' : (-80, 5),
    'P02787' : (25, -15),
    'P12004' : (25, 5),
    'P55263' : (-80, -15),
    'P62310' : (-80, -25),
    'Q14195' : (25, 15),
    'Q15417' : (25, 5),
}
for prot, row in robust_df.iterrows():
    x = row['log2FC']
    y = row['nlp']
    xoff, yoff = offsets.get(prot, (20, 10))
    ax.annotate(
        prot, (x, y),
        fontsize=7.5, color='#0d4f6e', fontweight='700',
        xytext=(xoff, yoff), textcoords='offset points',
        ha='left' if xoff > 0 else 'right',
        arrowprops=dict(arrowstyle='->', color='#2196A8', lw=0.9)
    )

ax.axhline(-np.log10(0.05), color='#888', linewidth=1, linestyle='--', alpha=0.7)
ax.axvline(-1, color='#888', linewidth=1, linestyle='--', alpha=0.7)
ax.axvline( 1, color='#888', linewidth=1, linestyle='--', alpha=0.7)
ax.set_xlabel("log₂ Fold Change  (Disease / Control)", fontsize=11, color='#333')
ax.set_ylabel("−log₁₀(p-value)", fontsize=11, color='#333')
ax.set_title("Robustness-Scored Pipeline\n(High confidence = robustness score ≥7/10)",
             fontsize=12, fontweight='bold', color='#222', pad=12)
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['left', 'bottom']].set_color('#CCC')
ax.tick_params(colors='#555')
ax.legend(fontsize=8.5, framealpha=0.9, edgecolor='#DDD',
          facecolor='white', loc='lower right')

plt.suptitle("Differential Protein Analysis — Classic vs Robustness Pipeline",
             fontsize=13, fontweight='bold', color='#1a1a1a', y=1.02)
plt.tight_layout(pad=2.5)
plt.savefig("figures/fig3_volcano.png", dpi=200,
            bbox_inches="tight", facecolor=BG)
plt.close()
print("\nFigure 3 saved → figures/fig3_volcano.png")
print("\nStep 3 COMPLETE")
