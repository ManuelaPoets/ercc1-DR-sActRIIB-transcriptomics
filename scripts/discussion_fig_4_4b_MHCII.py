# ════════════════════════════════════════════════════════════════════
# Fig 4.4b — Immune remodelling: MHC-II suppression vs macrophage identity (§4.4)
# Requires the shared-setup cell (atlas_stats, dis_table, save_dis, CONDS, colours).
# Kidney IMM + muscle macrophages. snRNA detection rates (exploratory, n=1).
# ════════════════════════════════════════════════════════════════════
detK, _, _ = atlas_stats(ATLAS_KIDNEY, ["H2-Aa", "H2-Ab1", "H2-Eb1", "Adgre1", "Csf1r"], "kidney")
detM, _, _ = atlas_stats(ATLAS_MUSCLE, ["H2-Aa", "H2-Ab1", "H2-Eb1", "Adgre1", "Csf1r"], "muscle")
MHC = ["H2-Aa", "H2-Ab1", "H2-Eb1"]; IDENT = ["Adgre1", "Csf1r"]
MHC_COLORS = [C_RED_D, C_CORAL_D, C_AMBER_D]; ID_COLORS = [C_TEAL_D, C_BLUE_D]
panels = [("Kidney immune cells (IMM)", detK, "IMM"), ("Muscle macrophages", detM, "Macrophages")]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4)); x = np.arange(len(CONDS))
for ax, (title, det, cell) in zip(axes, panels):
    for gn, col in zip(MHC, MHC_COLORS):
        ax.plot(x, det(cell, gn), marker="o", color=col, lw=1.6, label=f"{gn} (MHC-II)")
    for gn, col in zip(IDENT, ID_COLORS):
        ax.plot(x, det(cell, gn), marker="s", ls="--", color=col, lw=1.3, label=f"{gn} (identity)")
    ax.set_title(title, fontsize=9, fontweight="bold", color=C_TEXT)
    ax.set_ylabel("detection rate (% nuclei)", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(CONDS, fontsize=8)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
fig.suptitle("Immune remodelling: MHC-II antigen presentation suppressed while macrophage identity persists (§4.4)",
             fontsize=10.5, fontweight="bold")
fig.text(0.5, 0.01,
         "snRNA detection rates, exploratory (n=1 per condition). MHC-II markers (solid) fall with aging and further under "
         "combi in both tissues; macrophage-identity markers (dashed) stay flat or rise — antigen presentation is "
         "suppressed, not macrophage abundance.", ha="center", va="bottom", fontsize=6.6, color=C_TEXT2)
fig.tight_layout(rect=[0, 0.06, 1, 0.95])
save_dis(fig, "Fig_4_4b_MHCII_immune_remodelling"); plt.show()

# numbers table
import pandas as pd
rows = {}
for cell, det, tis in [("IMM", detK, "kidney"), ("Macrophages", detM, "muscle")]:
    for gn in MHC + IDENT:
        rows[f"{tis}:{gn}"] = det(cell, gn)
tbl = pd.DataFrame(rows, index=CONDS).T
dis_table(tbl, "MHC-II & macrophage-identity detection rate (%)",
          "Exploratory snRNA (n=1). Rows: tissue:gene; columns: condition.",
          "Tab_4_4b_MHCII_detection")
