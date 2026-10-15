import numpy as np
np.random.seed(42) 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.metrics import roc_auc_score
from sklearn.impute import KNNImputer

# ── Color palette ──────────────────────────────────────────────────────────
C_HIGH   = "#2196A8"
C_MID    = "#E8A838"
C_RISK   = "#C44E52"
C_HIGH_L = "#D6EEF1"
C_MID_L  = "#FDF3DC"
C_RISK_L = "#FAE0E0"
BG       = "#FAFAFA"
AX_BG    = "#F7F7F7"

meta    = pd.read_csv("data/metadata.csv", index_col=0)
df_norm = pd.read_csv("data/intensity_norm.csv", index_col=0)

ctrl_cols = [c for c in meta[meta['group']=='Control'].index if c in df_norm.columns]
dis_cols  = [c for c in meta[meta['group']=='Disease'].index  if c in df_norm.columns]
miss_pct  = df_norm.isna().mean(axis=1) * 100

scenarios   = ['KNN_noBatch','KNN_Batch','MinProb_noBatch',
               'MinProb_Batch','Median_noBatch','Median_Batch']
all_results = {s: pd.read_csv(f"results/deresult_{s}.csv", index_col=0)
               for s in scenarios}

# KNN imputed matrix for ROC
imp = KNNImputer(n_neighbors=5)
df_knn = pd.DataFrame(
    imp.fit_transform(df_norm.T).T,
    index=df_norm.index, columns=df_norm.columns
)
labels_binary = [1 if meta.loc[s,'group']=='Disease' else 0
                 for s in df_knn.columns]

# ── Robustness score ───────────────────────────────────────────────────────
print("Calculating robustness scores...")
scores_list = []

for protein in df_norm.index:
    # Skip contaminants
    if str(protein).startswith('CON__') or str(protein).startswith('REV__'):
        continue

    # Criterion 1: significant in >=4/6 scenarios
    n_sig = sum(
        all_results[s].loc[protein,'significant']
        if protein in all_results[s].index else False
        for s in scenarios
    )
    c1 = 3 if n_sig >= 4 else (1 if n_sig >= 2 else 0)

    # Criterion 2: significant after batch correction
    batch_sc = ['KNN_Batch','MinProb_Batch','Median_Batch']
    n_batch  = sum(
        all_results[s].loc[protein,'significant']
        if protein in all_results[s].index else False
        for s in batch_sc
    )
    c2 = 2 if n_batch >= 2 else (1 if n_batch == 1 else 0)

    # Criterion 3: low missing value rate
    miss = miss_pct.get(protein, 100)
    c3 = 2 if miss < 30 else (1 if miss < 50 else 0)

    # Criterion 4: fold change direction consistent
    fcs = [
        all_results[s].loc[protein,'log2FC']
        for s in scenarios
        if protein in all_results[s].index
        and not np.isnan(all_results[s].loc[protein,'log2FC'])
    ]
    if len(fcs) >= 4:
        pos = sum(f > 0 for f in fcs)
        neg = sum(f < 0 for f in fcs)
        c4 = 2 if (pos==len(fcs) or neg==len(fcs)) else (1 if max(pos,neg)>=len(fcs)-1 else 0)
    else:
        c4 = 0

    # Criterion 5: ROC-AUC > 0.75
    c5 = 0
    if protein in df_knn.index:
        try:
            auc = roc_auc_score(labels_binary, df_knn.loc[protein].values)
            c5  = 1 if max(auc, 1-auc) > 0.75 else 0
        except:
            c5 = 0

    total = c1 + c2 + c3 + c4 + c5
    cat   = ("High confidence" if total >= 7
             else "Moderate confidence" if total >= 4
             else "Likely artifact")

    scores_list.append({
        'protein'      : protein,
        'n_sig_of_6'   : n_sig,
        'c1_frequency' : c1,
        'c2_batch'     : c2,
        'c3_missing'   : c3,
        'c4_fc_dir'    : c4,
        'c5_roc'       : c5,
        'total_score'  : total,
        'category'     : cat,
        'missing_pct'  : round(miss, 1),
        'mean_log2FC'  : round(np.mean(fcs) if fcs else 0, 3),
    })

score_df = pd.DataFrame(scores_list).set_index('protein')
score_df.to_csv("results/robustness_scores.csv")

# ── Summary ────────────────────────────────────────────────────────────────
# Classic pipeline significant proteins from the reference scenario
classic_sig = set(all_results['MinProb_noBatch']
                  [all_results['MinProb_noBatch']['significant']].index)

# Remove contaminants/reverse hits so the table and Figure 4 use the same universe
classic_sig = {p for p in classic_sig
               if not str(p).startswith('CON__')
               and not str(p).startswith('REV__')}

high_conf    = set(score_df[score_df['category']=='High confidence'].index)
mid_conf     = set(score_df[score_df['category']=='Moderate confidence'].index)
risky        = set(score_df[score_df['category']=='Likely artifact'].index)

k_high = len(classic_sig & high_conf)
k_mid  = len(classic_sig & mid_conf)
k_risk = len(classic_sig & risky)

classic_total = k_high + k_mid + k_risk
robust_total  = k_high + k_mid
artifact_pct  = (k_risk / classic_total * 100) if classic_total else 0

print("=" * 52)
print("ROBUSTNESS SCORE SUMMARY")
print("=" * 52)
print(f"Total proteins analyzed    : {len(score_df)}")
print(f"High confidence  (>=7/10) : {len(high_conf)}")
print(f"Moderate conf.   (4-6/10) : {len(mid_conf)}")
print(f"Likely artifact  (<4/10)  : {len(risky)}")
print()
print(f"Classic pipeline flagged   : {classic_total} proteins")
print(f"  → High confidence        : {k_high}")
print(f"  → Moderate confidence    : {k_mid}")
print(f"  → Likely artifact        : {k_risk} ({artifact_pct:.1f}%)  ← removed")
print(f"Robustness-scored retained : {robust_total} proteins")

# Save the same values used in Figure 4, to avoid report/table mismatches
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
        k_high,
        k_mid,
        f'{k_risk} ({artifact_pct:.1f}%)'
    ]
})
summary_table.to_csv('results/robustness_summary_table.csv', index=False)

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Robustness barplot
# ══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
fig.patch.set_facecolor(BG)

# ── Left: Top 20 robustness barplot ────────────────────────────────────────
ax = axes[0]
ax.set_facecolor(AX_BG)

top20 = score_df.nlargest(20, 'total_score')
cat_colors = {
    'High confidence'     : C_HIGH,
    'Moderate confidence' : C_MID,
    'Likely artifact'     : C_RISK,
}
bar_colors = [cat_colors[c] for c in top20['category']]

bars = ax.barh(range(len(top20)), top20['total_score'],
               color=bar_colors, edgecolor='white',
               linewidth=1.0, height=0.72)

# Score labels
for bar, sc in zip(bars, top20['total_score']):
    ax.text(bar.get_width() + 0.12,
            bar.get_y() + bar.get_height()/2,
            f"{sc}/10", va='center', fontsize=9,
            color='#333333', fontweight='500')

ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20.index, fontsize=8.5, color='#222222')
ax.set_xlabel("Robustness Score (0–10)", fontsize=11, color='#333333')
ax.set_title("Top 20 Proteins — Robustness Score", fontsize=12,
             fontweight='bold', color='#222222', pad=12)
ax.axvline(7, color=C_RISK, linewidth=1.5, linestyle='--',
           alpha=0.7, label='Threshold ≥7 (high confidence)')
ax.axvline(4, color=C_MID,  linewidth=1.2, linestyle=':',
           alpha=0.7, label='Threshold ≥4 (moderate)')
ax.set_xlim(0, 12)
ax.invert_yaxis()
ax.spines[['top','right']].set_visible(False)
ax.spines[['left','bottom']].set_color('#CCCCCC')
ax.tick_params(colors='#555555')

legend_patches = [
    mpatches.Patch(color=C_HIGH, label='High confidence (≥7)'),
    mpatches.Patch(color=C_MID,  label='Moderate confidence (4–6)'),
    mpatches.Patch(color=C_RISK, label='Likely artifact (<4)'),
]
ax.legend(handles=legend_patches, fontsize=8.5, loc='lower right',
          framealpha=0.9, edgecolor='#DDDDDD', facecolor='white')

# ── Right: Classic vs Robust stacked bar ───────────────────────────────────
ax = axes[1]
ax.set_facecolor(AX_BG)

x      = [0, 1]
width  = 0.5
labels_x = ['Classic\nPipeline', 'Robustness-Scored\nPipeline']

# Classic bar
ax.bar([0], [k_high], width=width, color=C_HIGH,   edgecolor='white', linewidth=1)
ax.bar([0], [k_mid],  width=width, color=C_MID,    edgecolor='white', linewidth=1,
       bottom=[k_high])
ax.bar([0], [k_risk], width=width, color=C_RISK,   edgecolor='white', linewidth=1,
       bottom=[k_high + k_mid])

# Robust bar (artifacts removed)
ax.bar([1], [k_high], width=width, color=C_HIGH,   edgecolor='white', linewidth=1)
ax.bar([1], [k_mid],  width=width, color=C_MID,    edgecolor='white', linewidth=1,
       bottom=[k_high])

# Total labels on top
total_classic = k_high + k_mid + k_risk
total_robust  = k_high + k_mid
ax.text(0, total_classic + 3, f"n={total_classic}", ha='center',
        fontsize=12, fontweight='bold', color='#222222')
ax.text(1, total_robust  + 3, f"n={total_robust}",  ha='center',
        fontsize=12, fontweight='bold', color='#222222')
ax.set_ylim(0, total_classic + 40)


ax.set_xticks([0, 1])
ax.set_xticklabels(labels_x, fontsize=11, color='#222222')
ax.set_ylabel("Number of Significant Proteins", fontsize=11, color='#333333')
ax.set_title("Classic vs Robustness-Scored Pipeline\nSignificant Protein Comparison",
             fontsize=12, fontweight='bold', color='#222222', pad=12)
ax.set_xlim(-0.5, 1.5)
ax.spines[['top','right']].set_visible(False)
ax.spines[['left','bottom']].set_color('#CCCCCC')
ax.tick_params(colors='#555555')

legend_patches2 = [
    mpatches.Patch(color=C_HIGH, label='High confidence (≥7/10)'),
    mpatches.Patch(color=C_MID,  label='Moderate confidence (4–6/10)'),
    mpatches.Patch(color=C_RISK, label='Likely artifact (<4/10)'),
]
ax.legend(handles=legend_patches2, fontsize=9, loc='upper right',
          framealpha=0.9, edgecolor='#DDDDDD', facecolor='white')

plt.suptitle("Robustness Score Analysis", fontsize=14,
             fontweight='bold', color='#1a1a1a', y=1.02)
plt.tight_layout(pad=2.5)
plt.savefig("figures/fig4_robustness_barplot.png", dpi=200,
            bbox_inches="tight", facecolor=BG)
plt.close()
print("\nFigure 4 saved → figures/fig4_robustness_barplot.png")

# ── Final biomarker table ──────────────────────────────────────────────────
final_table = score_df[score_df['category'].isin(['High confidence','Moderate confidence'])]\
    .sort_values('total_score', ascending=False).head(20)\
    [['n_sig_of_6','missing_pct','mean_log2FC','total_score','category']]
final_table.columns = ['Significant scenarios (of 6)', 'Missing %',
                       'Mean log2FC', 'Robustness Score', 'Category']
final_table.to_csv("results/final_biomarker_table.csv")
print("Final table saved → results/final_biomarker_table.csv")
print("\nStep 4 COMPLETE")
print("\nAll outputs:")
print("  figures/fig1_pca.png")
print("  figures/fig2_missing_heatmap.png")
print("  figures/fig3_volcano.png")
print("  figures/fig4_robustness_barplot.png")
print("  results/robustness_scores.csv")
print("  results/final_biomarker_table.csv")
