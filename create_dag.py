#!/usr/bin/env python3
"""
Generate a publication-quality DAG figure for the Qualimap-HPV16 QC pipeline.
Requires: graphviz (Python package) + dot (system binary)
"""

import graphviz

dot = graphviz.Digraph(
    "Qualimap_HPV16_DAG",
    format="png",
    graph_attr={
        "rankdir": "TB",
        "fontname": "Helvetica",
        "fontsize": "14",
        "label": "Qualimap-HPV16 v1.1 — Pipeline DAG\nHPV16 Sequencing Quality Control",
        "labelloc": "t",
        "labeljust": "c",
        "bgcolor": "white",
        "dpi": "200",
        "pad": "0.5",
        "nodesep": "0.5",
        "ranksep": "0.7",
    },
)

# --- Styles ---
input_style = {
    "shape": "folder",
    "style": "filled",
    "fillcolor": "#E8F5E9",
    "color": "#388E3C",
    "fontname": "Helvetica",
    "fontsize": "10",
}
extract_style = {
    "shape": "box",
    "style": "filled,rounded",
    "fillcolor": "#E3F2FD",
    "color": "#1565C0",
    "fontname": "Helvetica-Bold",
    "fontsize": "10",
    "penwidth": "2",
}
qc_style = {
    "shape": "box",
    "style": "filled,rounded",
    "fillcolor": "#FFF3E0",
    "color": "#E65100",
    "fontname": "Helvetica-Bold",
    "fontsize": "10",
    "penwidth": "2",
}
hpv16_style = {
    "shape": "box",
    "style": "filled,rounded",
    "fillcolor": "#FCE4EC",
    "color": "#C62828",
    "fontname": "Helvetica-Bold",
    "fontsize": "10",
    "penwidth": "2",
}
classify_style = {
    "shape": "box",
    "style": "filled,rounded",
    "fillcolor": "#EDE7F6",
    "color": "#4527A0",
    "fontname": "Helvetica-Bold",
    "fontsize": "10",
    "penwidth": "2",
}
aggregate_style = {
    "shape": "box",
    "style": "filled,rounded",
    "fillcolor": "#E0F2F1",
    "color": "#00695C",
    "fontname": "Helvetica-Bold",
    "fontsize": "10",
    "penwidth": "2",
}
output_style = {
    "shape": "folder",
    "style": "filled",
    "fillcolor": "#FFEBEE",
    "color": "#C62828",
    "fontname": "Helvetica",
    "fontsize": "10",
}

edge_main = {"color": "#1565C0", "penwidth": "1.8"}
edge_qc = {"color": "#E65100", "penwidth": "1.8"}
edge_hpv = {"color": "#C62828", "penwidth": "1.8"}
edge_cls = {"color": "#4527A0", "penwidth": "1.8", "style": "dashed"}
edge_agg = {"color": "#00695C", "penwidth": "1.8"}

# --- Input ---
dot.node("input_bam", "Input BAM files\n(HPV-aligned reads)", **input_style)
dot.node("feature_file", "Feature file\n(hpv16_gene_annot.bed)", **input_style)
dot.node("rkmh_data", "rkmh reference data", **input_style)

# --- extractHPVandUNMAP ---
dot.node(
    "extract",
    "extractHPVandUNMAP\n─────────────────────\nExtract HPV-mapped\n+ unmapped reads",
    **extract_style,
)

# --- QC processes ---
dot.node(
    "qualimap",
    "qualimap\n─────────────────────\nQualimap bamqc\n(alignment QC metrics)",
    **qc_style,
)
dot.node(
    "flagstat",
    "flagstat\n─────────────────────\nsamtools flagstat\n(alignment flag stats)",
    **qc_style,
)

# --- HPV16 check ---
dot.node(
    "hpv16check",
    "hpv16check\n─────────────────────\nPer-gene read counting\n(E6,E7,E1,E2,E5a,L2,L1)\n+ FPKM normalization",
    **hpv16_style,
)

# --- Classification ---
dot.node(
    "bam2fq",
    "bam2fq\n─────────────────────\nBAM → FASTQ conversion",
    **classify_style,
)
dot.node(
    "rkmh",
    "rkmh\n─────────────────────\nSub-lineage classification\n(kmer-based)",
    **classify_style,
)

# --- Depth ---
dot.node(
    "depth",
    "hpv16ReadDepth\n─────────────────────\nDepth statistics\n(mean, Q1, Q3, median, std)\n+ coverage analysis",
    **hpv16_style,
)

# --- Aggregation ---
dot.node(
    "multiqc",
    "multiqc\n─────────────────────\nAggregate Qualimap\n+ flagstat reports",
    **aggregate_style,
)
dot.node(
    "summary",
    "summaryTable\n─────────────────────\nMerge all QC metrics\ninto summary_table.txt",
    **aggregate_style,
)

# --- Outputs ---
dot.node("multiqc_out", "multiqc_report.html\n+ multiqc_data/", **output_style)
dot.node("summary_out", "summary_table.txt\n(all metrics per sample)", **output_style)
dot.node("individual_out", "individual_reports/\n(per-sample BAMs,\nstats, depth)", **output_style)

# --- Edges: Input → extract ---
dot.edge("input_bam", "extract", **edge_main)

# --- extract → 3 branches ---
dot.edge("extract", "qualimap", label="  hpv_bam  ", **edge_qc)
dot.edge("extract", "flagstat", label="  hpv_bam2  ", **edge_qc)
dot.edge("extract", "hpv16check", label="  hpv_bam3  ", **edge_hpv)
dot.edge("feature_file", "qualimap", style="dotted", color="#E65100")

# --- hpv16check → branches ---
dot.edge("hpv16check", "bam2fq", label="  prim_hpv16.bam  ", **edge_cls)
dot.edge("hpv16check", "depth", label="  prim_hpv16.bam  ", **edge_hpv)
dot.edge("bam2fq", "rkmh", label="  .fastq  ", **edge_cls)
dot.edge("rkmh_data", "rkmh", style="dotted", color="#4527A0")

# --- QC → multiqc ---
dot.edge("qualimap", "multiqc", label="  qualimap_results  ", **edge_agg)
dot.edge("flagstat", "multiqc", label="  flagstat_results  ", **edge_agg)

# --- multiqc → summary ---
dot.edge("multiqc", "summary", label="  QC tables  ", **edge_agg)
dot.edge("multiqc", "multiqc_out", **edge_agg)

# --- hpv16check, rkmh, depth → summary ---
dot.edge("hpv16check", "summary", label="  hpv16_results  ", **edge_hpv)
dot.edge("rkmh", "summary", label="  rkmh_out  ", **edge_cls)
dot.edge("depth", "summary", label="  depth_results  ", **edge_hpv)

# --- summary → output ---
dot.edge("summary", "summary_out", **edge_agg)

# --- individual outputs ---
dot.edge("extract", "individual_out", style="dotted", color="#9E9E9E")
dot.edge("hpv16check", "individual_out", style="dotted", color="#9E9E9E")
dot.edge("depth", "individual_out", style="dotted", color="#9E9E9E")

# --- Legend ---
with dot.subgraph(name="cluster_legend") as legend:
    legend.attr(
        label="Legend",
        style="dashed",
        color="gray",
        fontname="Helvetica-Bold",
        fontsize="11",
    )
    legend.node("leg1", "Read extraction", shape="box", style="filled,rounded",
                fillcolor="#E3F2FD", color="#1565C0", fontsize="9", fontname="Helvetica")
    legend.node("leg2", "QC metrics", shape="box", style="filled,rounded",
                fillcolor="#FFF3E0", color="#E65100", fontsize="9", fontname="Helvetica")
    legend.node("leg3", "HPV16 analysis", shape="box", style="filled,rounded",
                fillcolor="#FCE4EC", color="#C62828", fontsize="9", fontname="Helvetica")
    legend.node("leg4", "Classification", shape="box", style="filled,rounded",
                fillcolor="#EDE7F6", color="#4527A0", fontsize="9", fontname="Helvetica")
    legend.node("leg5", "Aggregation", shape="box", style="filled,rounded",
                fillcolor="#E0F2F1", color="#00695C", fontsize="9", fontname="Helvetica")
    legend.edge("leg1", "leg2", style="invis")
    legend.edge("leg2", "leg3", style="invis")
    legend.edge("leg3", "leg4", style="invis")
    legend.edge("leg4", "leg5", style="invis")

output_path = "/tmp/Qualimap/dag"
dot.render(output_path, cleanup=True)
print(f"DAG saved to {output_path}.png")
