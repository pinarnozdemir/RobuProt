import numpy as np
np.random.seed(42) 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer

# ── Color palette ──────────────────────────────────────────────────────────
COLORS = {
    "Control": "#4C72B0",
    "Disease": "#C44E52",
    "Batch1" : "#55A868",
    "Batch2" : "#CCB974",
    "Batch3" : "#8172B2",
}
PALETTE = ["#4C72B0", "#C44E52", "#55A868", "#CCB974", "#8172B2", "#64B5CD"]

df   = pd.read_csv("data/intensity_matrix.csv", index_col=0)
meta = pd.read_csv("data/metadata.csv", index_col=0)

# ── Filtering ──────────────────────────────────────────────────────────────
miss_pct = df.isna().mean(axis=1) * 100
df_filt  = df[miss_pct <= 60].copy()
print(f"Before filtering : {len(df)} proteins")
print(f"After filtering  : {len(df_filt)} proteins (<=60% missing)")

# ── Median normalization ───────────────────────────────────────────────────
df_norm = df_filt.subtract(df_filt.median(axis=0), axis=1)
df_norm.to_csv("data/intensity_norm.csv")

miss_protein = df_filt.isna().mean(axis=1) * 100
miss_sample  = df_filt.isna().mean(axis=0) * 100
print(f"\nMissing value summary:")
print(f"  Per protein (mean) : {miss_protein.mean():.1f}%")
print(f"  Per sample  (mean) : {miss_sample.mean():.1f}%")
print(f"  Proteins >30% miss : {(miss_protein > 30).sum()}")

# ── KNN imputation for PCA ────────────────────────────────────────────────
imp = KNNImputer(n_neighbors=5)
df_knn = pd.DataFrame(
    imp.fit_transform(df_norm.T).T,
    index=df_norm.index, columns=df_norm.columns
)

# ── Batch correction ───────────────────────────────────────────────────────
df_bc = df_norm.copy()
for batch in meta['batch'].unique():
    b_cols = [c for c in meta[meta['batch']==batch].index if c in df_bc.columns]
    df_bc[b_cols] = df_bc[b_cols].subtract(
        df_bc[b_cols].mean(axis=1) - df_norm.mean(axis=1), axis=0
    )
df_bc.to_csv("data/intensity_batch_corrected.csv")

df_bc_imp = pd.DataFrame(
    imp.fit_transform(df_bc.T).T,
    index=df_bc.index, columns=df_bc.columns
)

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 1 — PCA
# ══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.patch.set_facecolor("#FAFAFA")

def draw_pca(ax, df_input, title):
    X  = df_input.T
    Xs = StandardScaler().fit_transform(X)
    pca    = PCA(n_components=2).fit(Xs)
    coords = pca.transform(Xs)
    var1, var2 = pca.explained_variance_ratio_ * 100

    markers = {"Batch1": "o", "Batch2": "s", "Batch3": "^"}
    for idx, sample in enumerate(X.index):
        if sample not in meta.index:
            continue
        grp   = meta.loc[sample, 'group']
        batch = meta.loc[sample, 'batch']
        ax.scatter(
            coords[idx, 0], coords[idx, 1],
            c=COLORS[grp], marker=markers[batch],
            s=110, edgecolors='white', linewidths=1.2, zorder=3, alpha=0.9
        )
        ax.annotate(
            sample, (coords[idx, 0], coords[idx, 1]),
            fontsize=6.5, alpha=0.65, xytext=(4, 4),
            textcoords='offset points', color='#333333'
        )

    ax.set_xlabel(f"PC1 ({var1:.1f}% variance)", fontsize=11, color='#333333')
    ax.set_ylabel(f"PC2 ({var2:.1f}% variance)", fontsize=11, color='#333333')
    ax.set_title(title, fontsize=12, fontweight='bold', pad=12, color='#222222')
    ax.set_facecolor("#F7F7F7")
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color('#CCCCCC')
    ax.grid(True, alpha=0.3, linewidth=0.6, color='white')
    ax.tick_params(colors='#555555')

    # Legend — group
    for grp, col in COLORS.items():
        if grp in ['Control', 'Disease']:
            ax.scatter([], [], c=col, s=70, label=grp)
    # Legend — batch
    for bat, mk in markers.items():
        ax.scatter([], [], c='gray', marker=mk, s=70, label=bat)
    ax.legend(fontsize=8.5, framealpha=0.9, edgecolor='#DDDDDD',
              facecolor='white', loc='best')

draw_pca(axes[0], df_knn,    "PCA — Before batch correction")
draw_pca(axes[1], df_bc_imp, "PCA — After batch correction")

plt.suptitle("Principal Component Analysis", fontsize=14,
             fontweight='bold', color='#1a1a1a', y=1.01)
plt.tight_layout(pad=2)
plt.savefig("figures/fig1_pca.png", dpi=200, bbox_inches="tight",
            facecolor='#FAFAFA')
plt.close()
print("\nFigure 1 saved → figures/fig1_pca.png")

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Missing Value Analysis
# ══════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(14, 6))
fig.patch.set_facecolor("#FAFAFA")
gs  = gridspec.GridSpec(1, 2, width_ratios=[1, 1.4], wspace=0.35)

# Left — bar chart
ax1 = fig.add_subplot(gs[0])
bins   = [0, 10, 20, 30, 40, 50, 60, 70]
labels = ['0–10%','10–20%','20–30%','30–40%','40–50%','50–60%','60–70%']
counts = [((miss_protein>=bins[i])&(miss_protein<bins[i+1])).sum()
          for i in range(len(bins)-1)]
bar_colors = ['#2ecc71','#27ae60','#f1c40f','#e67e22','#e74c3c','#c0392b','#922b21']
bars = ax1.bar(labels, counts, color=bar_colors, edgecolor='white',
               linewidth=1.2, width=0.7)
for bar, cnt in zip(bars, counts):
    if cnt > 0:
        ax1.text(bar.get_x()+bar.get_width()/2, bar.get_height()+8,
                 str(cnt), ha='center', fontsize=9, color='#333333', fontweight='500')
ax1.axvline(x=2.5, color='#E24B4A', linestyle='--',
            linewidth=1.8, label='30% filter threshold', alpha=0.85)
ax1.set_xlabel("Missing value rate", fontsize=11, color='#333333')
ax1.set_ylabel("Number of proteins", fontsize=11, color='#333333')
ax1.set_title("Missing Value Distribution\nper Protein", fontsize=12,
              fontweight='bold', color='#222222')
ax1.set_facecolor("#F7F7F7")
ax1.spines[['top', 'right']].set_visible(False)
ax1.spines[['left', 'bottom']].set_color('#CCCCCC')
ax1.tick_params(axis='x', rotation=35, colors='#555555')
ax1.tick_params(axis='y', colors='#555555')
ax1.legend(fontsize=9, framealpha=0.9)

# Right — heatmap
ax2 = fig.add_subplot(gs[1])
top60  = miss_protein.sort_values(ascending=False).head(60).index
hm     = df_filt.loc[top60].isna().astype(int)
sorder = meta.sort_values(['batch','group']).index.tolist()
sorder = [s for s in sorder if s in hm.columns]
hm     = hm[sorder]

from matplotlib.colors import LinearSegmentedColormap
cmap_custom = LinearSegmentedColormap.from_list(
    'miss', ['#2ecc71', '#FFFFFF', '#e74c3c'], N=256
)
im = ax2.imshow(hm.values, aspect='auto', cmap='RdYlGn_r',
                vmin=0, vmax=1, interpolation='nearest')
ax2.set_xticks(range(len(sorder)))
ax2.set_xticklabels(sorder, rotation=65, fontsize=7.5,
                    ha='right', color='#333333')
ax2.set_yticks([])
ax2.set_title("Missing Value Heatmap\n(Top 60 proteins, batch-ordered)",
              fontsize=12, fontweight='bold', color='#222222')

cum = 0
batch_colors = {'Batch1':'#55A868','Batch2':'#CCB974','Batch3':'#8172B2'}
for bat in meta['batch'].unique():
    n = (meta.loc[sorder,'batch']==bat).sum()
    ax2.axvline(x=cum-0.5, color='#333333', linewidth=1.8)
    ax2.text(cum+n/2-0.5, -3.5, bat, ha='center', fontsize=8.5,
             color=batch_colors.get(bat,'#333333'), fontweight='bold')
    cum += n

cbar = plt.colorbar(im, ax=ax2, fraction=0.025, pad=0.02, shrink=0.85)
cbar.set_ticks([0, 1])
cbar.set_ticklabels(['Observed', 'Missing'], fontsize=9)
cbar.ax.tick_params(colors='#555555')

plt.suptitle("Missing Value Analysis", fontsize=14,
             fontweight='bold', color='#1a1a1a', y=1.01)
plt.savefig("figures/fig2_missing_heatmap.png", dpi=200,
            bbox_inches="tight", facecolor='#FAFAFA')
plt.close()
print("Figure 2 saved → figures/fig2_missing_heatmap.png")
print("\nStep 2 COMPLETE")
