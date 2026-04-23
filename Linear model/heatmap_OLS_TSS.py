import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from scipy.stats import zscore, variation
import seaborn as sns
import matplotlib as mpl

#jobid=Plt-G
#sbatch --export=All -t 48:00:00 --mem=150G  -e run_log/${jobid}.err -o run_log/${jobid}.out -J $jobid --wrap="python heatmap_OLS_TSS.py"

df=pd.read_excel("OLS_data_plot_TSS.xlsx", sheet_name="Sheet1", index_col=0)  # index_col=0 if you saved the index

# --- Select only expression columns ---
expr_cols = ["female_central", "female_proximal-distal", "male_central", "male_proximal-distal"]
expr_data = df[expr_cols]

# Convert to numpy array
expr_array = expr_data[expr_cols].values.astype(float)

# Row-wise z-score
zscored_array = zscore(expr_array, axis=1)

# Convert back to DataFrame
zscored = pd.DataFrame(zscored_array, index=expr_data.index, columns=expr_cols)

# --- Row color mapping based on assignment ---
palette = {
    "Condition": "orchid",
    "Annotation": "darkslategrey",
    "Interaction": "darkgoldenrod",
    "unassigned": "silver"
}
row_colors = df["assignment"].map(palette)

# Extremely important: lower raster DPI for heatmap
mpl.rcParams["figure.dpi"] = 60     # text still vector
plt.rcParams['pdf.fonttype'] = 42   # Illustrator-compatible fonts
plt.rcParams['svg.fonttype'] = 'none'

sns.set(style="white", font_scale=1.0)

g = sns.clustermap(
    zscored,
    row_colors=row_colors,
    cmap="Spectral_r",
    center=0,
    linewidths=0,
    figsize=(10, 12),
    yticklabels=False,
    cbar_kws={'label': "Z-score counts"}
)

# ------------------------------
# RASTERIZE ONLY THE HEATMAP
# ------------------------------

for im in g.ax_heatmap.get_images():
    im.set_rasterized(True)

# Optional: row colors also rasterized (recommended for 100k rows)
for im in g.ax_row_colors.get_images():
    im.set_rasterized(True)

# Text stays vector
for ax in [g.ax_heatmap, g.ax_row_colors,
           g.ax_col_dendrogram, g.ax_row_dendrogram, g.cax]:
    for txt in ax.texts:
        txt.set_rasterized(False)

# Legend vector
for label, color in palette.items():
    g.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
g.ax_col_dendrogram.legend(
    loc="center",
    ncol=3,
    bbox_to_anchor=(0.5, 1.1),
    frameon=False,
    title="Assignment"
)

# ------------------------------
# SAVE AS PDF (Illustrator friendly)
# ------------------------------

plt.savefig(
    "OLS_1kbBins_TSS.pdf",
    dpi=60,                # << the embedded PNG resolution
    bbox_inches="tight",
    metadata={}            # strip metadata → smaller file
)
plt.close()
