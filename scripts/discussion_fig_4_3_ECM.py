# ════════════════════════════════════════════════════════════════════
# Fig 4.3 — ECM remodelling: collagen & protease balance in stromal cells
# Discussion §4.3.  Kidney fibroblasts (FIB) and muscle FAPs.
# Computed straight from the snRNA atlases (version-independent — snRNA was
# NOT re-processed for muscle v2). Metric = mean log-normalised expression
# (raw/X) within celltype x condition — reproduces the §4.3 values exactly.
# ════════════════════════════════════════════════════════════════════
import h5py, gc
from scipy.sparse import csr_matrix

# >>> SET these to your Drive h5ad locations <<<
ATLAS_KIDNEY = f"{DATA_BASE}/.../annotated-aging-soupxed.h5ad"    # kidney, 56,686 nuclei
ATLAS_MUSCLE = f"{DATA_BASE}/.../annotated-muscle-soupxed.h5ad"   # muscle, 24,419 nuclei
DIS_OUT      = OUTPUT_DIR_DIS                                     # change here for a different folder

# muscle sample IDs -> condition (CLAUDE.md: rgzj1=ctrl, 2=age, 3=sAct, 4=DR, 5=combi)
MUSCLE_SMAP = {"rgzj1": "ctrl", "rgzj2": "age", "rgzj3": "sAct", "rgzj4": "DR", "rgzj5": "combi"}

def mean_lognorm(path, celltype_value, genes, conditions, sample_map=None):
    """Mean log-norm expression per gene per condition within one cell type,
    read directly from raw/X via h5py (no scanpy/anndata required)."""
    f = h5py.File(path, "r")
    g = f["raw"]["X"]; shape = tuple(g.attrs["shape"])
    X = csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).tocsc()
    var = np.array([x.decode() if isinstance(x, bytes) else x
                    for x in f["raw"]["var"]["_index"][:]])
    gidx = {gn: int(np.where(var == gn)[0][0]) for gn in genes}
    def cat(col):
        o = f["obs"][col]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in o["categories"][:]])
        return cats[o["codes"][:]]
    ct = cat("celltype"); samp = cat("sample")
    f.close()
    cond = np.array([sample_map.get(s, s) for s in samp]) if sample_map else samp
    res = {}
    for gn in genes:
        col = np.asarray(X[:, gidx[gn]].todense()).ravel()
        res[gn] = {c: float(col[(ct == celltype_value) & (cond == c)].mean()) for c in conditions}
    del X; gc.collect()
    return res

CONDS = ["ctrl", "age", "DR", "sAct", "combi"]
GENES = ["Col3a1", "Col4a1", "Mmp2", "Timp2"]
kid = mean_lognorm(ATLAS_KIDNEY, "FIB",  GENES, CONDS)
mus = mean_lognorm(ATLAS_MUSCLE, "FAPs", GENES, CONDS, sample_map=MUSCLE_SMAP)

# ── Bulk significance (n=3, DESeq2) for the RESCUE contrast (combi vs age) ──
# Hard-coded from the canonical bulk tables (res_intv_df), verified locally.
# Stars: *** p<.001, ** p<.01, * p<.05, ns = not significant.
#   KIDNEY combi-vs-age padj: Col3a1 7.0e-5 ***, Col4a1 0.97 ns, Mmp2 2.8e-5 ***, Timp2 0.82 ns
#   MUSCLE combi-vs-age padj: Col3a1 8.2e-17 ***, Col4a1 0.044 *, Mmp2 1.4e-6 ***, Timp2 0.034 *
BULK_RESCUE_SIG = {
    "kidney": {"Col3a1": "***", "Col4a1": "ns", "Mmp2": "***", "Timp2": "ns"},
    "muscle": {"Col3a1": "***", "Col4a1": "*",  "Mmp2": "***", "Timp2": "*"},
}

# ── Plot: 2 rows (kidney FIB / muscle FAP) x 2 cols (collagen / MMP) ──
GENE_COLORS = {"Col3a1": C_CORAL_D, "Col4a1": C_BLUE_D,   # red = interstitial/fibrotic, blue = basement membrane
               "Mmp2":   C_AMBER_D, "Timp2":  C_TEAL_D}   # amber = protease, teal = inhibitor
panels = [("Kidney fibroblasts (FIB)", "kidney", kid), ("Muscle FAPs", "muscle", mus)]
pairs  = [("Col3a1", "Col4a1"), ("Mmp2", "Timp2")]

fig, axes = plt.subplots(2, 2, figsize=(10, 7.4))
x = np.arange(len(CONDS)); w = 0.38
for r, (row_title, tkey, data) in enumerate(panels):
    for c, (num, den) in enumerate(pairs):
        ax = axes[r, c]
        sig = BULK_RESCUE_SIG[tkey]
        ax.bar(x - w/2, [data[num][cd] for cd in CONDS], w, color=GENE_COLORS[num],
               label=f"{num}  (bulk {sig[num]})")
        ax.bar(x + w/2, [data[den][cd] for cd in CONDS], w, color=GENE_COLORS[den],
               label=f"{den}  (bulk {sig[den]})")
        ax.set_xticks(x); ax.set_xticklabels(CONDS, fontsize=8)
        ax.set_ylabel("mean log-norm expr.", fontsize=8)
        ax.set_title(f"{row_title} — {num}/{den}", fontsize=9, fontweight="bold", color=C_TEXT)
        ax.legend(fontsize=7, frameon=False, loc="upper left")
        ratio = [data[num][cd] / data[den][cd] for cd in CONDS]
        ax2 = ax.twinx()
        ax2.plot(x, ratio, color=C_TEXT, marker="o", ms=4, lw=1.3, zorder=5)
        ax2.axhline(ratio[0], ls="--", lw=0.8, color=C_TEXT2)   # ctrl reference
        ax2.set_ylabel(f"{num}/{den} ratio", fontsize=8, color=C_TEXT)
        ax2.set_ylim(0, max(ratio) * 1.35)
        for xi, rv in zip(x, ratio):
            ax2.annotate(f"{rv:.2f}", (xi, rv), textcoords="offset points",
                         xytext=(0, 5), ha="center", fontsize=6.5, color=C_TEXT)

fig.suptitle("ECM remodelling: collagen and protease balance in stromal cells (§4.3)",
             fontsize=11, fontweight="bold")
fig.text(0.5, 0.015,
         "Bars/ratios: mean log-norm expression in the snRNA atlas (exploratory, n=1 per condition). "
         "'(bulk ***/**/*/ns)' = significance of that gene in the bulk RNA-seq rescue contrast\n"
         "(combi vs age; n=3, DESeq2 padj). Significance derives from the bulk data; the single-nucleus "
         "panels localise the change to fibroblasts/FAPs.",
         ha="center", va="bottom", fontsize=6.8, color=C_TEXT2)
fig.tight_layout(rect=[0, 0.055, 1, 0.96])
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(DIS_OUT, f"Fig_4_3_ECM_collagen_MMP_balance.{ext}"),
                dpi=DPI, bbox_inches="tight", facecolor="white")
print("Saved Fig_4_3_ECM_collagen_MMP_balance to", DIS_OUT)
plt.show()
