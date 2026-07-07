# ════════════════════════════════════════════════════════════════════
# Fig 4.4a — Immune–endothelial rescue hub: Tgfb1 & Flt1 (§4.4)
# Built on BULK RNA-seq (n=3, DESeq2) — the statistically valid direction.
# snRNA is used only to note WHERE the genes are expressed (localisation),
# not to claim a cell-type-specific direction. Requires shared-setup cell.
# ════════════════════════════════════════════════════════════════════
import pandas as pd

# Bulk LFC + padj, verified from res_age_df (aging) / res_intv_df (combi vs age):
HUB = [   # gene, tissue, aging_LFC, aging_padj, rescue_LFC, rescue_padj
    ("Tgfb1", "kidney",  0.46, 1.1e-2, -0.07, 6.2e-1),
    ("Tgfb1", "muscle",  1.37, 1.4e-10, -1.30, 8.2e-10),
    ("Flt1",  "kidney", -0.53, 6.0e-5,  0.78, 1.6e-8),
]
def star(p): return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "ns"

labels = [f"{g}\n({t})" for g, t, *_ in HUB]
aging  = [h[2] for h in HUB]; rescue = [h[4] for h in HUB]
fig, ax = plt.subplots(figsize=(7.6, 4.7)); x = np.arange(len(HUB)); w = 0.38
ax.bar(x - w/2, aging,  w, color=C_RED_D,  label="aging (age vs ctrl)")
ax.bar(x + w/2, rescue, w, color=C_TEAL_D, label="rescue (combi vs age)")
ax.axhline(0, color=C_TEXT, lw=0.8)
for xi, h in zip(x, HUB):
    ax.annotate(star(h[3]), (xi - w/2, h[2]), textcoords="offset points",
                xytext=(0, 3 if h[2] >= 0 else -11), ha="center", fontsize=8, color=C_TEXT)
    ax.annotate(star(h[5]), (xi + w/2, h[4]), textcoords="offset points",
                xytext=(0, 3 if h[4] >= 0 else -11), ha="center", fontsize=8, color=C_TEXT)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("log₂ fold change (bulk, n=3)", fontsize=8)
ax.set_title("Immune–endothelial hub: Tgfb1 & Flt1 dysregulated with aging and reversed by combi (§4.4)",
             fontsize=9, fontweight="bold", color=C_TEXT)
ax.legend(fontsize=7, frameon=False)
fig.text(0.5, 0.005,
         "Bulk RNA-seq (n=3, DESeq2); opposite-sign aging vs rescue bars = directional reversal. "
         "snRNA localisation (baseline detection): Tgfb1 in kidney IMM ~33% & muscle macrophages ~25%; "
         "Flt1 in kidney endothelium ~80–98% — genes mark the immune/endothelial compartments "
         "(direction = bulk; cell-type expression = exploratory).",
         ha="center", va="bottom", fontsize=6.3, color=C_TEXT2)
fig.tight_layout(rect=[0, 0.07, 1, 0.95])
save_dis(fig, "Fig_4_4a_immune_endothelial_hub"); plt.show()

# numbers table
tbl = pd.DataFrame({
    "aging LFC":  [h[2] for h in HUB],
    "aging sig":  [f"{h[3]:.1e} {star(h[3])}" for h in HUB],
    "rescue LFC": [h[4] for h in HUB],
    "rescue sig": [f"{h[5]:.1e} {star(h[5])}" for h in HUB],
}, index=[f"{g} ({t})" for g, t, *_ in HUB])
dis_table(tbl, "Immune–endothelial hub genes — bulk RNA-seq (n=3)",
          "Rescue = combi vs age. Significance from bulk (n=3); snRNA used only for localisation.",
          "Tab_4_4a_hub_bulk")
