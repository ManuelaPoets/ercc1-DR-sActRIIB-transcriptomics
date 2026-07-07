# ════════════════════════════════════════════════════════════════════
# Fig 4.5 — Proximal tubule (PT-2): damage, identity loss & depletion (§4.5)
# Kidney snRNA atlas. Detection rate = % nuclei with expression > 0 within
# celltype x condition. Exploratory (n=1). Bulk (n=3) context in the caption.
# ════════════════════════════════════════════════════════════════════
import h5py, gc
from scipy.sparse import csr_matrix

ATLAS_KIDNEY = f"{DATA_BASE}/.../annotated-aging-soupxed.h5ad"   # >>> SET to your Drive path <<<
DIS_OUT = OUTPUT_DIR_DIS
CONDS = ["ctrl", "age", "DR", "sAct", "combi"]

def load_kidney_pt(path, genes):
    """detection rate (%) per (celltype, gene, condition) + cell counts per celltype."""
    f = h5py.File(path, "r")
    g = f["raw"]["X"]; shape = tuple(g.attrs["shape"])
    X = csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).tocsc()
    var = np.array([x.decode() if isinstance(x, bytes) else x for x in f["raw"]["var"]["_index"][:]])
    gidx = {gn: int(np.where(var == gn)[0][0]) for gn in genes}
    def cat(c):
        o = f["obs"][c]
        cats = np.array([z.decode() if isinstance(z, bytes) else z for z in o["categories"][:]])
        return cats[o["codes"][:]]
    ct = cat("celltype"); cond = cat("sample"); f.close()
    cols = {gn: np.asarray(X[:, gidx[gn]].todense()).ravel() for gn in genes}
    del X; gc.collect()
    def det(cell, gn):
        return [100*np.mean(cols[gn][(ct == cell) & (cond == c)] > 0)
                if ((ct == cell) & (cond == c)).sum() else np.nan for c in CONDS]
    def count(cell):
        return [int(((ct == cell) & (cond == c)).sum()) for c in CONDS]
    return det, count

det, count = load_kidney_pt(ATLAS_KIDNEY, ["Havcr1", "Cdkn1a", "Slc34a1", "Vcam1"])
COND_COLORS = ['#888888', C_RED_D, C_GREEN_D, C_PURP_D, C_TEAL_D]

fig, axes = plt.subplots(1, 3, figsize=(13, 4.3)); x = np.arange(len(CONDS))
# A — PT-2 damage / senescence / identity
ax = axes[0]
ax.plot(x, det("PT-2", "Havcr1"),  marker="o", color=C_RED_D,  lw=1.6, label="Havcr1 (injury)")
ax.plot(x, det("PT-2", "Cdkn1a"),  marker="s", color=C_PURP_D, lw=1.6, label="Cdkn1a (senescence)")
ax.plot(x, det("PT-2", "Slc34a1"), marker="^", color=C_BLUE_D, lw=1.6, label="Slc34a1 (PT identity)")
ax.set_title("PT-2: injury, senescence & identity", fontsize=9, fontweight="bold", color=C_TEXT)
ax.set_ylabel("detection rate (% nuclei)", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(CONDS, fontsize=8); ax.legend(fontsize=7, frameon=False)
# B — Vcam1 failed-repair, PT-1 vs PT-2
ax = axes[1]
ax.plot(x, det("PT-1", "Vcam1"), marker="o", color=C_GRAY_D,  lw=1.6, label="PT-1")
ax.plot(x, det("PT-2", "Vcam1"), marker="o", color=C_CORAL_D, lw=1.6, label="PT-2")
ax.set_title("Vcam1 (failed-repair state)", fontsize=9, fontweight="bold", color=C_TEXT)
ax.set_ylabel("detection rate (% nuclei)", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(CONDS, fontsize=8); ax.legend(fontsize=7, frameon=False)
# C — PT-2 nuclei recovered per condition
ax = axes[2]
cc = count("PT-2"); ax.bar(x, cc, color=COND_COLORS)
for xi, v in zip(x, cc):
    ax.annotate(str(v), (xi, v), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=7, color=C_TEXT)
ax.set_title("PT-2 nuclei recovered", fontsize=9, fontweight="bold", color=C_TEXT)
ax.set_ylabel("n nuclei", fontsize=8); ax.set_xticks(x); ax.set_xticklabels(CONDS, fontsize=8)

fig.suptitle("Proximal tubule PT-2: a selective, layered response to rescue (§4.5)",
             fontsize=11, fontweight="bold")
fig.text(0.5, 0.005,
         "snRNA detection rates & counts: exploratory (n=1 per condition). Bulk RNA-seq (n=3): Havcr1, Cdkn1a, Vcam1 "
         "significantly up and Slc34a1 down with aging (padj<1e-9), and NONE significantly rescued (combi vs age, ns) — "
         "bulk confirms the persistence; snRNA localises it to the depleted PT-2 subpopulation.",
         ha="center", va="bottom", fontsize=6.6, color=C_TEXT2)
fig.tight_layout(rect=[0, 0.07, 1, 0.95])
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(DIS_OUT, f"Fig_4_5_PT2_damage_identity.{ext}"),
                dpi=DPI, bbox_inches="tight", facecolor="white")
print("Saved Fig_4_5_PT2_damage_identity to", DIS_OUT)
plt.show()
