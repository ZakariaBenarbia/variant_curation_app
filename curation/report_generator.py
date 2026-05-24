# =============================================================================
# GÉNÉRATEUR DE RAPPORTS CLINIQUES
# =============================================================================
# Ce module génère des rapports HTML détaillés pour l'analyse
# de variants génomiques avec visualisations et tableaux interactifs.
#
# Fonctionnalités principales:
# - Génération de graphiques (VAF, rôles de gènes, origine)
# - Création de tableaux HTML avec filtres interactifs
# - Intégration de liens externes (GeneCards, COSMIC)
# - Export en format HTML responsive et moderne
#
# Visualisations générées:
# - Distribution VAF (histogramme + zoom)
# - Rôles de gènes (barres horizontales)
# - Origine des mutations (camembert + barres)
#
# Sources de données:
# - Données de variants scorées par l'IA
# - Métadonnées COSMIC et GeneCards
# - Statistiques cliniques et filtres dynamiques

import pandas as pd
import datetime
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np


# ═══════════════════════════════════════════════════════════════════════════════
#  PALETTE & STYLE
# ═══════════════════════════════════════════════════════════════════════════════
PALETTE = {
    "bg":      "#F9FAFB", "panel":  "#FFFFFF", "grid":   "#E5E7EB",
    "text":    "#1F2937", "sub":    "#6B7280",
    "blue":    "#3B82F6", "purple": "#A855F7", "amber":  "#F59E0B",
    "red":     "#DC2626", "teal":   "#0D9488", "green":  "#16A34A",
    "orange":  "#EA8C55",
}

plt.rcParams.update({
    "figure.facecolor": PALETTE["bg"],  "axes.facecolor":   PALETTE["panel"],
    "axes.edgecolor":   PALETTE["grid"], "axes.labelcolor": PALETTE["text"],
    "xtick.color":      PALETTE["sub"], "ytick.color":      PALETTE["sub"],
    "text.color":       PALETTE["text"], "grid.color":       PALETTE["grid"],
    "grid.linestyle":   "--",           "grid.alpha":       0.7,
    "font.family":      "DejaVu Sans",  "axes.spines.top":  False,
    "axes.spines.right": False,
})


# ═══════════════════════════════════════════════════════════════════════════════
#  LINK HELPERS  — all clickable columns are defined here
# ═══════════════════════════════════════════════════════════════════════════════
_LINK = (
    "color:#0A4D8C;text-decoration:none;font-weight:600;"
    "border-bottom:1px dashed #0A4D8C;transition:color .15s"
)


def link_gene(gene):
    """GeneCards link for a gene symbol."""
    gene = str(gene).strip()
    if not gene or gene.lower() in ("nan", "intergenic", ""):
        return gene
    url = f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}"
    return f'<a href="{url}" target="_blank" style="{_LINK}" title="GeneCards: {gene}">{gene}</a>'


def format_cosmic(cosmic_id):
    """Display COSMIC ID as plain text (non-clickable)."""
    val = str(cosmic_id).strip()
    if not val or val.lower() in ("nan", "n/a", "", "unreported"):
        return "—"
    return val


# ═══════════════════════════════════════════════════════════════════════════════
#  PRIORITY BADGE
# ═══════════════════════════════════════════════════════════════════════════════
def priority_badge(row, origin_col, score_col):
    """Generate badge using shared classification logic"""
    priority_class = classify_variant_priority(row, score_col, origin_col)
    
    badge_classes = {
        0: ("low", "Signification clinique inconnue"),
        1: ("moderate", "Signification clinique potentielle"),
        2: ("critical", "Signification clinique très forte")
    }
    
    cls, label = badge_classes[priority_class]
    return f'<span class="badge {cls}">{label}</span>'


# ═══════════════════════════════════════════════════════════════════════════════
#  CLASSIFICATION HELPER (same logic as priority_badge and cards)
# ═══════════════════════════════════════════════════════════════════════════════
def classify_variant_priority(row, score_col="priority_score", origin_col="origin_improved"):
    """
    Classifie un variant en 3 niveaux de priorité clinique.
    Même logique que priority_badge() et classify_for_count().
    Retourne: 0=Inconnue, 1=Potentielle, 2=Très forte
    """
    mut_status = row.get("mut_status", "")
    gene = str(row.get("gene", "")).strip()
    origin = str(row.get(origin_col, "")).lower()
    score = row.get(score_col, 0) or 0
    
    if gene == "Intergenic":
        return 0  # Inconnue
    elif pd.notna(mut_status) and "Confirmed somatic variant" in str(mut_status):
        return 2  # Très forte
    elif pd.isna(mut_status) or mut_status == "":
        return 1  # Potentielle
    elif "loh" in origin or score > 0.6:
        return 2  # Très forte
    elif score > 0.3:
        return 1  # Potentielle
    else:
        return 0  # Inconnue


# ═══════════════════════════════════════════════════════════════════════════════
#  CHARTS
# ═══════════════════════════════════════════════════════════════════════════════
def make_vaf_plot(df, path):
    """Histogramme de distribution VAF - version originale"""
    fig = plt.figure(figsize=(13, 5.5), facecolor=PALETTE["bg"])
    fig.suptitle("Distribution VAF", fontsize=14, fontweight="bold",
                 color=PALETTE["text"], x=0.02, ha="left", y=1.01)

    if "vaf" not in df.columns:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "Colonne VAF non trouvée", ha="center", va="center",
                fontsize=13, color=PALETTE["sub"])
        ax.set_axis_off()
        fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=PALETTE["bg"])
        plt.close(fig)
        return

    vaf = df["vaf"].dropna().astype(float)
    bins = np.linspace(0, 1, 60)

    ax = fig.add_axes([0.06, 0.14, 0.54, 0.72])
    n, edges, patches = ax.hist(vaf, bins=bins, edgecolor=PALETTE["bg"], linewidth=0.4)
    for patch, left in zip(patches, edges[:-1]):
        patch.set_facecolor(
            PALETTE["red"] if left >= 0.85 else PALETTE["amber"] if left >= 0.3 else PALETTE["blue"]
        )
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda x, _: f"{int(x):,}" if x >= 1 else "")
    )
    ax.set_xlabel("Fréquence Allélique du Variant (VAF)", fontsize=11)
    ax.set_ylabel("Nombre de variants (échelle log)", fontsize=11)
    ax.set_xlim(0, 1)
    ax.grid(True, axis="y", alpha=0.5)
    for vline, label, col in [
        (0.5, "Hétérozygote (0,5)", PALETTE["teal"]),
        (0.85, "Seuil LOH (0,85)", PALETTE["red"]),
    ]:
        ax.axvline(vline, color=col, linestyle="--", linewidth=1.4, alpha=0.9, label=label)
    ax.legend(fontsize=9, framealpha=0.6, facecolor=PALETTE["panel"], edgecolor=PALETTE["grid"])

    ax2 = fig.add_axes([0.66, 0.14, 0.31, 0.72])
    bins2 = np.linspace(0, 0.5, 40)
    _, _, patches2 = ax2.hist(vaf[vaf <= 0.5], bins=bins2, edgecolor=PALETTE["bg"], linewidth=0.4)
    for patch, left in zip(patches2, bins2[:-1]):
        patch.set_facecolor(PALETTE["amber"] if left >= 0.3 else PALETTE["blue"])
    ax2.set_xlabel("VAF (zoom : 0 - 0,5)", fontsize=10)
    ax2.set_ylabel("Comptage (linéaire)", fontsize=10)
    ax2.set_xlim(0, 0.5)
    ax2.grid(True, axis="y", alpha=0.5)
    ax2.set_title("Détail VAF faible", fontsize=10, color=PALETTE["sub"], pad=6)

    fig.legend(
        handles=[
            mpatches.Patch(color=PALETTE["blue"], label="VAF < 0,30"),
            mpatches.Patch(color=PALETTE["amber"], label="0,30 - 0,85"),
            mpatches.Patch(color=PALETTE["red"], label="VAF ≥ 0,85  (zone LOH)"),
        ],
        loc="lower center", ncol=3, fontsize=9, framealpha=0.6,
        facecolor=PALETTE["panel"], edgecolor=PALETTE["grid"],
        bbox_to_anchor=(0.5, -0.04)
    )
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=PALETTE["bg"])
    plt.close(fig)


def make_gene_role_chart(df, path):
    """Répartition des variants par catégorie fonctionnelle - Signification très forte uniquement"""
    fig, ax = plt.subplots(figsize=(10, 5), facecolor=PALETTE["bg"])
    fig.patch.set_facecolor(PALETTE["bg"])
    col = next((c for c in ["gene_role", "gene_roles"] if c in df.columns), None)

    if not col:
        ax.text(0.5, 0.5, "Colonne gene_role non trouvée", ha="center", va="center",
                fontsize=13, color=PALETTE["sub"])
        ax.set_axis_off()
    else:
        # Filtrer uniquement les variants à signification très forte (classe 2)
        score_col = next((c for c in ["improved_score", "priority_score", "score"] if c in df.columns), None)
        origin_col = "origin_improved" if "origin_improved" in df.columns else "origin"
        df["priority_class"] = df.apply(lambda row: classify_variant_priority(row, score_col, origin_col), axis=1)
        df_critical = df[df["priority_class"] == 2]  # Très forte uniquement
        
        all_counts = df_critical[col].value_counts()
        unknown = all_counts.get("Unknown", 0)
        counts = all_counts.drop(labels=["Unknown"], errors="ignore")

        if counts.empty:
            ax.text(0.5, 0.5,
                    "Aucun rôle de gène connu trouvé\n(tous les variants classés comme Inconnu)",
                    ha="center", va="center", fontsize=12,
                    color=PALETTE["sub"], wrap=True)
            ax.set_axis_off()
        else:
            def get_role_color(role):
                role = str(role).strip().lower()
                # Check for dual/combined roles first
                if ("oncogene" in role and ("tsg" in role or "suppressor" in role)) or "dual" in role:
                    return PALETTE["blue"]
                if "oncogene" in role:
                    return PALETTE["purple"]
                if "suppressor" in role or "tsg" in role:
                    return PALETTE["orange"]
                return PALETTE["teal"]
            colours = [get_role_color(r) for r in counts.index]
            bars = ax.barh(counts.index[::-1], counts.values[::-1],
                           color=colours[::-1], height=0.55,
                           edgecolor=PALETTE["bg"], linewidth=0.8)
            total = counts.sum()
            for bar in bars:
                width = bar.get_width()
                ax.text(width + counts.max() * 0.015,
                        bar.get_y() + bar.get_height() / 2,
                        f"{int(width):,}  ({width/total*100:.1f}%)",
                        va="center", ha="left", fontsize=10, color=PALETTE["text"])
            ax.set_xlabel("Nombre de variants", fontsize=11)
            ax.set_xlim(0, counts.max() * 1.28)
            ax.grid(True, axis="x", alpha=0.5)
            ax.set_title("Rôles des gènes - Signification très forte uniquement",
                         fontsize=13, fontweight="bold", color=PALETTE["text"], pad=12)

        if unknown:
            fig.text(0.98, 0.02,
                     f"† {int(unknown):,} variants avec rôle de gène inconnu exclus",
                     ha="right", va="bottom", fontsize=8,
                     color=PALETTE["sub"], style="italic")

    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=PALETTE["bg"])
    plt.close(fig)


def make_origin_chart(df, path, origin_col):
    """Répartition des variants par origine - tous les variants"""
    fig, (ax_donut, ax_bar) = plt.subplots(
        1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1.2, 1]},
        facecolor=PALETTE["bg"]
    )
    fig.patch.set_facecolor(PALETTE["bg"])
    fig.suptitle("Classification de l'origine des gènes", fontsize=14, fontweight="bold",
                 color=PALETTE["text"], x=0.02, ha="left")

    def origin_colour(label):
        l = label.lower()
        if "strictly somatic" in l:
            return "#B91C1C"
        if "somatic" in l and "loh" in l:
            return PALETTE["red"]
        if "somatic" in l:
            return PALETTE["orange"]
        if "dual" in l or "ambiguous" in l:
            return PALETTE["amber"]
        if "strictly germline" in l:
            return "#1E3A8A"
        if "germline" in l:
            return PALETTE["blue"]
        return PALETTE["teal"]

    if origin_col not in df.columns:
        for ax in (ax_donut, ax_bar):
            ax.text(0.5, 0.5, "Colonne origin non trouvée",
                    ha="center", va="center", fontsize=12, color=PALETTE["sub"])
            ax.set_axis_off()
    else:
        # Use the full origin distribution for all variants (exclude only Unknown and Associated)
        all_counts = df[origin_col].fillna("Unknown").value_counts()
        excl_unk = all_counts.get("Unknown", 0)
        excl_assoc = all_counts.get("Associated", 0)
        counts = all_counts.drop(labels=["Unknown", "Associated"], errors="ignore")
        total = counts.sum()
        colours = [origin_colour(s) for s in counts.index]

        if excl_unk or excl_assoc:
            fig.text(0.98, 0.02,
                     f"† {int(excl_unk):,} inconnu + {int(excl_assoc):,} associé exclus",
                     ha="right", va="bottom", fontsize=8,
                     color=PALETTE["sub"], style="italic")

        wedges, _, autotexts = ax_donut.pie(
            counts.values, colors=colours, autopct="%1.1f%%", pctdistance=0.78,
            startangle=90, wedgeprops={"edgecolor": PALETTE["bg"], "linewidth": 2.5}
        )
        ax_donut.add_patch(plt.Circle((0, 0), 0.52, fc=PALETTE["bg"], linewidth=0))
        ax_donut.text(0, 0.10, f"{total:,}", ha="center", va="center",
                      fontsize=14, fontweight="bold", color=PALETTE["text"])
        ax_donut.text(0, -0.18, "Total", ha="center", va="center",
                      fontsize=9, color=PALETTE["sub"])
        for i, (at, wedge) in enumerate(zip(autotexts, wedges)):
            pct = counts.values[i] / total * 100
            if pct >= 8:
                at.set_color("white")
                at.set_fontsize(9)
                at.set_fontweight("bold")
            else:
                angle = np.deg2rad((wedge.theta1 + wedge.theta2) / 2)
                at.set_position((1.25 * np.cos(angle), 1.25 * np.sin(angle)))
                at.set_color(PALETTE["text"])
                at.set_fontsize(8.5)
                at.set_fontweight("bold")
        ax_donut.legend(
            handles=[
                mpatches.Patch(
                    color=colours[i],
                    label=f"{counts.index[i]} — {counts.values[i]:,} ({counts.values[i]/total*100:.1f}%)"
                )
                for i in range(len(counts))
            ],
            loc="lower center", bbox_to_anchor=(0.5, -0.32), ncol=2, fontsize=8.5,
            framealpha=0.7, facecolor=PALETTE["panel"], edgecolor=PALETTE["grid"]
        )
        ax_donut.set_aspect("equal")

        cs = counts.sort_values(ascending=True)
        cum = 0
        segments = []
        y_levels = [1.5, -1.5, 2.2, -2.2, 0.8, -0.8]
        for val, lbl, col in zip(cs.values, cs.index, [origin_colour(s) for s in cs.index]):
            ax_bar.barh(0, val, left=cum, height=0.45,
                        color=col, edgecolor=PALETTE["bg"], linewidth=1.2)
            segments.append((cum + val / 2, val, lbl, val / total * 100))
            cum += val
        ax_bar.set_xlim(0, total)
        ax_bar.set_ylim(-2.4, 2.4)
        for i, (cx, val, lbl, pct) in enumerate(segments):
            y = y_levels[i % len(y_levels)]
            ax_bar.annotate(
                f"{lbl.split('(')[0].strip()}\n{val:,} ({pct:.1f}%)",
                xy=(cx, 0.25 if y > 0 else -0.25), xytext=(cx, y),
                ha="center", va="center", fontsize=8, color=PALETTE["text"],
                arrowprops=dict(arrowstyle="-", color=PALETTE["grid"], lw=0.8),
                bbox=dict(boxstyle="round,pad=0.3", fc="white",
                          ec=PALETTE["grid"], lw=0.6, alpha=0.9)
            )
        ax_bar.set_xlabel("Nombre de variants", fontsize=10)
        ax_bar.set_title("Répartition proportionnelle", fontsize=11, color=PALETTE["text"], pad=10)
        ax_bar.set_yticks([])
        ax_bar.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
        ax_bar.grid(True, axis="x", alpha=0.4)
        for sp in ["left", "top", "right"]:
            ax_bar.spines[sp].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=PALETTE["bg"])
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════
def main(input_file=None, output_file=None, vaf_plot_path=None,
         gene_role_chart_path=None, loh_pie_path=None, sample_id="sample"):
    try:
        input_file = snakemake.input.scores
        output_file = snakemake.output.report
        vaf_plot_path = snakemake.output.vaf_plot
        gene_role_chart_path = snakemake.output.gene_role_bar_chart
        loh_pie_path = snakemake.output.loh_pie_chart
        sample_id = snakemake.wildcards.sample
    except NameError:
        pass

    print(f"[generate_report] Sample: {sample_id}")
    for lbl, val in [("Input", input_file), ("Report", output_file),
                     ("VAF", vaf_plot_path), ("Roles", gene_role_chart_path),
                     ("Origin", loh_pie_path)]:
        print(f"  {lbl}: {val}")

    for p in [output_file, vaf_plot_path, gene_role_chart_path, loh_pie_path]:
        if p:
            os.makedirs(os.path.dirname(p), exist_ok=True)

    df = pd.read_csv(input_file)
    print(f"[generate_report] Columns: {list(df.columns)}")

    if "vaf" in df.columns:
        if df["vaf"].dtype == object:
            df["vaf"] = df["vaf"].str.replace("%", "", regex=False).astype(float)
        if df["vaf"].max() > 1.0:
            df["vaf"] = df["vaf"] / 100.0

    origin_col = "origin_improved" if "origin_improved" in df.columns else "origin"
    gene_role_col = next((c for c in ["gene_role", "gene_roles"] if c in df.columns), None)
    score_col = next((c for c in ["improved_score", "final_score", "score", "priority_score"] if c in df.columns), None)
    if score_col:
        df = df.sort_values(by=score_col, ascending=False)

    df["_gene_link"] = df["gene"].apply(link_gene) if "gene" in df.columns else "—"
    df["_cosmic_fmt"] = df["cosmic_id"].apply(format_cosmic) if "cosmic_id" in df.columns else "—"
    df["_priority"] = df.apply(lambda r: priority_badge(r, origin_col, score_col), axis=1)

    for label, fn, args in [
        ("VAF plot", make_vaf_plot, (df, vaf_plot_path)),
        ("Gene role", make_gene_role_chart, (df, gene_role_chart_path)),
        ("Origin chart", make_origin_chart, (df, loh_pie_path, origin_col)),
    ]:
        print(f"[generate_report] Generating {label}...")
        fn(*args)
        print("  ✓ Done")

    # AI Score column excluded from display — used internally for sorting/classification only
    col_order = [
        "_priority", "_gene_link", gene_role_col, origin_col,
        "tumor_type", "loh_status",
        "chrom", "pos", "ref", "alt", "vaf", "depth",
        "_cosmic_fmt", "mut_status",
    ]
    display_cols = [c for c in col_order if c and c in df.columns]
    rename = {
        "_priority": "Priority",
        "_gene_link": "Gene",
        gene_role_col: "Gene Role",
        origin_col: "Origin",
        "tumor_type": "Tumor Type",
        "loh_status": "LOH Status",
        "chrom": "Chr",
        "pos": "Position",
        "ref": "Ref",
        "alt": "Alt",
        "vaf": "VAF",
        "depth": "Depth",
        "_cosmic_fmt": "COSMIC ID",
        "mut_status": "Mutation Status",
    }
    df_display = df[display_cols].rename(columns={k: v for k, v in rename.items() if k in display_cols})

    html_table = df_display.to_html(
        classes="vtable", index=False, escape=False, table_id="vtable"
    )

    n_total = len(df)
    
    # Count variants using shared classify_variant_priority function
    # This ensures perfect consistency between cards and badges
    classifications = df.apply(lambda row: classify_variant_priority(row, score_col, origin_col), axis=1)
    n_critical = sum(classifications == 2)
    n_moderate = sum(classifications == 1)
    n_loh = (len(df[df["loh_status"].str.contains("Likely LOH", na=False, case=False)])
             if "loh_status" in df.columns else 0)
    n_somatic = (len(df[df[origin_col].str.contains("Somatic", na=False, case=False)])
                 if origin_col in df.columns else 0)
    # Cartes rôles COSMIC : prefer improved_score when present
    score_roles = next((c for c in ["improved_score", "priority_score", "score"] if c in df.columns), None)
    if gene_role_col and gene_role_col in df.columns and score_roles:
        role_crit = df.apply(
            lambda row: classify_variant_priority(row, score_roles, origin_col), axis=1
        )
        roles = df.loc[role_crit == 2, gene_role_col].astype(str).str.strip()
        n_onco = int((roles == "Oncogene").sum())
        n_tsg = int((roles == "Tumor Suppressor").sum())
        n_dual = int((roles == "Oncogene/TSG").sum())
    else:
        n_onco = n_tsg = n_dual = 0

    now = datetime.datetime.now()

    vaf_col_idx = (df_display.columns.get_loc("VAF")
                   if "VAF" in df_display.columns else -1)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Clinical Genomic Report — {sample_id}</title>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;500;600;700&family=Spectral:ital,wght@0,600;0,700;1,400&display=swap" rel="stylesheet">
<style>
*{{margin:0;padding:0;box-sizing:border-box}}
:root{{
  --ink:#0D1117;--ink2:#4B5563;--ink3:#9CA3AF;
  --blue:#0A4D8C;--blue-lt:#EBF3FB;
  --red:#DC2626;--red-lt:#FEE2E2;
  --amber:#D97706;--amber-lt:#FEF3C7;
  --gray-lt:#F9FAFB;--border:#E5E7EB;--border2:#D1D5DB;
  --purple:#7C3AED;--teal:#0D9488;--orange:#EA580C;
  --radius:10px;--shadow:0 1px 3px rgba(0,0,0,.08),0 4px 12px rgba(0,0,0,.06);
}}
body{{font-family:'IBM Plex Sans',sans-serif;background:#F0F4F8;color:var(--ink);line-height:1.6;padding:2rem 1.5rem;min-height:100vh}}
a{{color:var(--blue);text-decoration:none;border-bottom:1px dashed var(--blue)}}
a:hover{{color:#073A6A;border-bottom-color:#073A6A}}
.wrap{{max-width:none;margin:0 auto;display:flex;flex-direction:column;gap:1.5rem}}
header{{
  background:linear-gradient(135deg,#0A4D8C 0%,#073A6A 100%);
  color:#fff;border-radius:var(--radius);
  box-shadow:0 8px 30px rgba(10,77,140,.25);
  overflow:hidden;position:relative;
}}
header::after{{
  content:'';position:absolute;top:-60px;right:-60px;
  width:320px;height:320px;border-radius:50%;
  background:radial-gradient(circle,rgba(255,255,255,.07) 0%,transparent 70%);
  pointer-events:none;
}}
.h-topbar{{
  display:flex;align-items:center;justify-content:space-between;
  padding:1rem 2rem;border-bottom:1px solid rgba(255,255,255,.1);
  background:rgba(0,0,0,.12);
}}
.h-logo{{display:flex;align-items:center;gap:.75rem;text-decoration:none;border:none}}
.h-logo-icon{{
  width:36px;height:36px;border-radius:8px;
  background:rgba(255,255,255,.18);border:1px solid rgba(255,255,255,.25);
  display:flex;align-items:center;justify-content:center;
  font-family:'IBM Plex Mono',monospace;font-size:1rem;font-weight:700;
  color:#fff;flex-shrink:0;
}}
.h-logo-text{{line-height:1.2}}
.h-logo-name{{font-weight:700;font-size:1rem;color:#fff;letter-spacing:-.01em}}
.h-logo-tagline{{font-size:.72rem;opacity:.7;font-weight:300;color:#fff}}
.h-badge{{
  font-size:.7rem;font-family:'IBM Plex Mono',monospace;font-weight:600;
  background:rgba(255,255,255,.12);border:1px solid rgba(255,255,255,.2);
  padding:.3rem .8rem;border-radius:20px;opacity:.85;letter-spacing:.04em;
}}
.h-body{{
  display:grid;grid-template-columns:1fr auto;gap:2rem;
  align-items:center;padding:2rem 2rem 2.2rem;
}}
.h-eyebrow{{
  font-size:.7rem;text-transform:uppercase;letter-spacing:.12em;
  opacity:.6;font-weight:600;margin-bottom:.5rem;
}}
.h-title{{font-family:'Spectral',serif;font-size:2rem;font-weight:700;letter-spacing:-.02em;line-height:1.2;margin-bottom:.4rem}}
.h-sub{{opacity:.75;font-size:.9rem;font-weight:300}}
.h-meta{{display:grid;grid-template-columns:1fr 1fr;gap:.6rem;min-width:280px}}
.meta-pill{{
  background:rgba(255,255,255,.1);border:1px solid rgba(255,255,255,.18);
  border-radius:7px;padding:.6rem 1rem;backdrop-filter:blur(6px);
}}
.meta-label{{font-size:.65rem;text-transform:uppercase;letter-spacing:.07em;opacity:.7;margin-bottom:.15rem}}
.meta-val{{font-family:'IBM Plex Mono',monospace;font-size:.9rem;font-weight:600}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:1rem}}
.card{{
  background:#fff;border-radius:var(--radius);padding:1.25rem 1.5rem;
  box-shadow:var(--shadow);border-top:3px solid var(--blue);transition:transform .2s
}}
.card:hover{{transform:translateY(-3px)}}
.card.red{{border-top-color:var(--red)}} .card.amber{{border-top-color:var(--amber)}}
.card.teal{{border-top-color:var(--teal)}} .card.blue2{{border-top-color:#3B82F6}}
.card.purple{{border-top-color:var(--purple)}} .card.orange{{border-top-color:var(--orange)}}
.card.dual{{border-top-color:var(--blue)}}
.card-num{{font-family:'IBM Plex Mono',monospace;font-size:2.2rem;font-weight:700;line-height:1;color:var(--blue)}}
.card.red .card-num{{color:var(--red)}}
.card.amber .card-num{{color:var(--amber)}}
.card.teal .card-num{{color:var(--teal)}}
.card.blue2 .card-num{{color:#3B82F6}}
.card.purple .card-num{{color:var(--purple)}}
.card.orange .card-num{{color:var(--orange)}}
.card.dual .card-num{{color:var(--blue)}}
.card-lbl{{font-size:.75rem;color:var(--ink2);text-transform:uppercase;letter-spacing:.05em;margin-top:.4rem;font-weight:500}}
.note{{
  background:var(--blue-lt);border-left:4px solid var(--blue);
  border-radius:0 var(--radius) var(--radius) 0;padding:1.25rem 1.75rem
}}
.note-title{{font-weight:700;color:var(--blue);font-size:.8rem;text-transform:uppercase;letter-spacing:.06em;margin-bottom:.4rem}}
.note-body{{font-size:.9rem;color:var(--ink);line-height:1.7}}
section{{background:#fff;border-radius:var(--radius);padding:2rem 2.5rem;box-shadow:var(--shadow)}}
.sec-title{{font-family:'Spectral',serif;font-size:1.5rem;font-weight:600;color:#073A6A;
  margin-bottom:1.5rem;padding-bottom:.75rem;border-bottom:2px solid var(--border)}}
.filters{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));
  gap:1rem;padding:1.25rem;background:var(--gray-lt);border-radius:8px;
  border:1px solid var(--border);margin-bottom:1.5rem}}
.filter-group label{{display:block;font-size:.75rem;font-weight:600;color:var(--ink2);
  text-transform:uppercase;letter-spacing:.05em;margin-bottom:.35rem}}
.filter-group input,.filter-group select{{
  width:100%;padding:.5rem .75rem;border:1px solid var(--border2);
  border-radius:6px;background:#fff;font-size:.875rem;color:var(--ink);
  font-family:'IBM Plex Sans',sans-serif;transition:border-color .15s
}}
.filter-group input:focus,.filter-group select:focus{{
  outline:none;border-color:var(--blue);box-shadow:0 0 0 3px rgba(10,77,140,.1)
}}
.btn{{
  display:inline-flex;align-items:center;gap:.4rem;
  padding:.48rem 1rem;border-radius:6px;cursor:pointer;font-weight:600;
  font-size:.8rem;font-family:'IBM Plex Sans',sans-serif;
  transition:all .18s;white-space:nowrap;
}}
.btn-clear{{
  background:#fff;color:var(--ink2);
  border:1.5px solid var(--border2);
}}
.btn-clear:hover{{background:#F3F4F6;border-color:#9CA3AF;color:var(--ink)}}
.btn-excel{{
  background:#16A34A;color:#fff;border:1.5px solid #16A34A;
  box-shadow:0 1px 3px rgba(22,163,74,.25);
}}
.btn-excel:hover{{background:#15803D;border-color:#15803D;box-shadow:0 3px 8px rgba(22,163,74,.35)}}
#filter-count{{font-size:.75rem;color:var(--ink3);font-weight:500;white-space:nowrap}}
.tbl-wrap{{overflow-x:visible;border-radius:8px;border:1px solid var(--border)}}
.vtable{{width:100%;border-collapse:separate;border-spacing:0;font-size:.85rem}}
.vtable thead th{{
  background:#F3F4F6;color:var(--ink2);font-weight:600;font-size:.72rem;
  text-transform:uppercase;letter-spacing:.06em;padding:.9rem 1.1rem;
  text-align:left;border-bottom:2px solid var(--border2);
  position:sticky;top:0;z-index:10;white-space:nowrap
}}
.vtable tbody tr{{background:#fff;border-bottom:1px solid var(--border);transition:background .15s}}
.vtable tbody tr:nth-child(even){{background:#FAFAFA}}
.vtable tbody tr:hover{{background:var(--blue-lt)}}
.vtable tbody td{{padding:1rem 1.1rem;vertical-align:middle}}
.vtable tbody td:nth-child(2){{font-family:'IBM Plex Mono',monospace;font-weight:600;font-size:.9rem}}
.vtable tbody td:nth-child(7),.vtable tbody td:nth-child(8),
.vtable tbody td:nth-child(9),.vtable tbody td:nth-child(10){{font-family:'IBM Plex Mono',monospace;font-size:.8rem}}
.badge{{display:inline-block;padding:.35rem .8rem;border-radius:5px;
  font-size:.68rem;font-weight:700;letter-spacing:.04em;text-transform:uppercase;border:1.5px solid}}
.badge.critical{{background:var(--red-lt);color:var(--red);border-color:var(--red)}}
.badge.moderate{{background:var(--amber-lt);color:var(--amber);border-color:var(--amber)}}
.badge.low{{background:#F3F4F6;color:#6B7280;border-color:var(--border2)}}
footer{{
  background:#fff;border-radius:var(--radius);padding:2rem 2.5rem;
  box-shadow:var(--shadow);display:grid;
  grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:2rem
}}
.ft-title{{font-weight:700;font-size:.8rem;color:#073A6A;text-transform:uppercase;
  letter-spacing:.05em;margin-bottom:.5rem}}
.ft-body{{font-size:.82rem;color:var(--ink2);line-height:1.6}}
.ft-legal{{grid-column:1/-1;border-top:1px solid var(--border);padding-top:1.25rem;
  font-size:.75rem;color:var(--ink3);line-height:1.6}}
@media(max-width:768px){{
  .h-body{{grid-template-columns:1fr;}}
  .h-meta{{grid-template-columns:1fr 1fr;min-width:unset}}
  .h-topbar{{flex-wrap:wrap;gap:.5rem}}
  section,footer{{padding:1.5rem}}
}}
@media print{{body{{background:#fff}} .vtable tbody tr:hover{{background:#fff}}}}
</style>
</head>
<body>
<div class="wrap">
<header>
  <div class="h-topbar">
    <div class="h-logo">
      <div class="h-logo-icon">V</div>
      <div class="h-logo-text">
        <div class="h-logo-name">VarCurate</div>
        <div class="h-logo-tagline">Curation clinique</div>
      </div>
    </div>
    <div class="h-badge">AI-Powered Pipeline v3.0 · GRCh38</div>
  </div>
  <div class="h-body">
    <div class="h-left">
      <div class="h-eyebrow">Genomic Analysis Report</div>
      <div class="h-title">Clinical Genomic Report</div>
      <div class="h-sub">Analyse des variants alléliques · Genomic Variant Analysis</div>
    </div>
    <div class="h-meta">
      <div class="meta-pill"><div class="meta-label">Sample ID</div><div class="meta-val">{sample_id}</div></div>
      <div class="meta-pill"><div class="meta-label">Analysis Date</div><div class="meta-val">{now.strftime('%Y-%m-%d')}</div></div>
      <div class="meta-pill"><div class="meta-label">Pipeline</div><div class="meta-val">AI-Powered v3.0</div></div>
      <div class="meta-pill"><div class="meta-label">Reference</div><div class="meta-val">GRCh38</div></div>
    </div>
  </div>
</header>
<div class="cards">
  <div class="card"><div class="card-num">{n_total}</div><div class="card-lbl">Total des variantes</div></div>
  <div class="card red"><div class="card-num">{n_critical}</div><div class="card-lbl">Signification très forte</div></div>
  <div class="card amber"><div class="card-num">{n_moderate}</div><div class="card-lbl">Signification potentielle</div></div>
  <div class="card teal"><div class="card-num">{n_loh}</div><div class="card-lbl">perte de l'hétérozygotie (LOH)</div></div>
  <div class="card blue2"><div class="card-num">{n_somatic}</div><div class="card-lbl">Mutations somatiques</div></div>
  <div class="card purple"><div class="card-num">{n_onco}</div><div class="card-lbl">Oncogènes</div></div>
  <div class="card orange"><div class="card-num">{n_tsg}</div><div class="card-lbl">Suppresseurs de tumeurs</div></div>
  <div class="card dual"><div class="card-num">{n_dual}</div><div class="card-lbl">Oncogène / TSG</div></div>
</div>
<section>
  <h2 class="sec-title">Detailed Variant Analysis</h2>
  <div class="filters">
    <div class="filter-group">
      <label>Priority</label>
      <select id="f-priority">
        <option value="">All</option>
        <option value="très forte">Signification très forte</option>
        <option value="potentielle">Signification potentielle</option>
        <option value="inconnue">Signification inconnue</option>
      </select>
    </div>
    <div class="filter-group">
      <label>Gene</label>
      <input id="f-gene" type="text" placeholder="e.g. TP53">
    </div>
    <div class="filter-group">
      <label>Origin</label>
      <select id="f-origin">
        <option value="">All</option>
        <option value="Somatic">Somatic</option>
        <option value="Strictly Somatic">Strictly Somatic</option>
        <option value="Germline">Germline</option>
        <option value="Dual">Dual</option>
      </select>
    </div>
    <div class="filter-group">
      <label>Min VAF</label>
      <input id="f-vaf" type="number" placeholder="0.00" min="0" max="1" step="0.01">
    </div>
    <div class="filter-group" style="display:flex;align-items:flex-end;gap:.5rem">
      <button class="btn btn-clear" id="btn-clear">
        <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round"><line x1="18" y1="6" x2="6" y2="18"/><line x1="6" y1="6" x2="18" y2="18"/></svg>
        Clear
      </button>
      <button class="btn btn-excel" id="btn-excel" title="Export visible rows to Excel">
        <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
        Export
      </button>
    </div>
  </div>
  <div style="margin-bottom:.75rem">
    <span id="filter-count" style="font-size:.8rem;color:var(--ink3);font-weight:500"></span>
  </div>
  <div class="tbl-wrap">{html_table}</div>
</section>
<footer>
  <div>
    <div class="ft-title">Reference Database</div>
    <div class="ft-body">Classification based on COSMIC Cancer Gene Census and Cosmic Complete Targeted Screen Mutant.</div>
  </div>
  <div>
    <div class="ft-title">Quality Control</div>
    <div class="ft-body">Minimum depth >=10x, base quality >=20, mapping quality >=30.</div>
  </div>
  <div>
    <div class="ft-title">Clickable Resources</div>
    <div class="ft-body">
      <strong>Gene</strong> -> GeneCards &nbsp;·&nbsp;
      <strong>COSMIC ID</strong> -> COSMIC mutation identifier
    </div>
  </div>
  <div class="ft-legal">
    <strong>For Research Use Only.</strong>
    This report is generated for academic and thesis purposes.
    Clinical decisions should not be made based solely on this analysis without proper validation
    and consultation with qualified healthcare professionals. &nbsp;|&nbsp;
    Generated {now.strftime('%Y-%m-%d %H:%M:%S')}
  </div>
</footer>
</div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js"></script>
<script>
(function(){{
  const tbl = document.getElementById('vtable');
  const rows = Array.from(tbl.tBodies[0].rows);
  const fPri = document.getElementById('f-priority');
  const fGene = document.getElementById('f-gene');
  const fOrig = document.getElementById('f-origin');
  const fVaf = document.getElementById('f-vaf');
  const count = document.getElementById('filter-count');

  function filter(){{
    let vis = 0;
    rows.forEach(r=>{{
      const td = r.cells;
      const pri = td[0]?.textContent.trim().toLowerCase() || '';
      const gene = td[1]?.textContent.trim().toLowerCase() || '';
      const orig = td[3]?.textContent.trim() || '';
      const vaf = parseFloat(td[{vaf_col_idx}]?.textContent) || 0;
      const ok = (!fPri.value || pri.includes(fPri.value.toLowerCase()))
              && (!fGene.value || gene.includes(fGene.value.toLowerCase()))
              && (!fOrig.value || orig.toLowerCase().includes(fOrig.value.toLowerCase()))
              && (!fVaf.value || vaf >= parseFloat(fVaf.value));
      r.style.display = ok ? '' : 'none';
      if(ok) vis++;
    }});
    count.textContent = `Showing ${{vis}} of ${{rows.length}} variants`;
  }}

  [fPri, fOrig].forEach(el => el.addEventListener('change', filter));
  [fGene, fVaf].forEach(el => el.addEventListener('input', filter));
  document.getElementById('btn-clear').addEventListener('click',()=>{{
    fPri.value=''; fGene.value=''; fOrig.value=''; fVaf.value=''; filter();
  }});
  count.textContent = `Showing ${{rows.length}} of ${{rows.length}} variants`;

  document.getElementById('btn-excel').addEventListener('click', () => {{
    // Collect headers from thead
    const headers = Array.from(tbl.tHead.rows[0].cells).map(th => th.textContent.trim());

    // Collect only visible rows, using plain text (strips HTML/badges)
    const data = [headers];
    rows.forEach(r => {{
      if (r.style.display === 'none') return;
      data.push(Array.from(r.cells).map(td => td.textContent.trim()));
    }});

    const ws = XLSX.utils.aoa_to_sheet(data);

    // Auto-size columns
    const colWidths = headers.map((_, ci) =>
      Math.min(60, Math.max(10, ...data.map(row => (row[ci] || '').length)))
    );
    ws['!cols'] = colWidths.map(w => ({{ wch: w }}));

    // Style header row bold (basic)
    headers.forEach((_, ci) => {{
      const cellRef = XLSX.utils.encode_cell({{r: 0, c: ci}});
      if (ws[cellRef]) ws[cellRef].s = {{ font: {{ bold: true }} }};
    }});

    const wb = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(wb, ws, 'Variants');

    const filename = `variants_{sample_id}_${{new Date().toISOString().slice(0,10)}}.xlsx`;
    XLSX.writeFile(wb, filename);
  }});
}})();
</script>
</body>
</html>"""

    try:
        with open(output_file, "w", encoding="utf-8") as handle:
            handle.write(html)
        print(f"[generate_report] ✓ Report → {output_file}")
        print(f"[generate_report] ✓ VAF    → {vaf_plot_path}")
        print(f"[generate_report] ✓ Roles  → {gene_role_chart_path}")
        print(f"[generate_report] ✓ Origin → {loh_pie_path}")
    except Exception as exc:
        print(f"[generate_report] ERROR: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    return output_file


try:
    if hasattr(snakemake, "log") and snakemake.log:
        sys.stdout = open(str(snakemake.log[0]), "w")
        sys.stderr = sys.stdout
    main()
except NameError:
    pass
except Exception as e:
    print(f"Snakemake execution error: {e}")
    raise
