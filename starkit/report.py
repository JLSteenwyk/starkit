"""
HTML report generation module.

Produces a self-contained HTML report with inline CSS, JS, and SVG gene
diagrams for each predicted Starship element.
"""

import logging
import os
from typing import List

import jinja2

from starkit.models import EvidenceLevel, StarKITRun, StarshipResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# SVG gene diagram
# ---------------------------------------------------------------------------

def generate_svg_diagram(
    starship_result: StarshipResult,
    width: int = 800,
    height: int = 80,
) -> str:
    """Generate an SVG string showing a linear gene map for one Starship.

    Draws the Starship region as a horizontal track with arrows for the
    captain and cargo genes and small markers for TIR boundaries.

    Parameters
    ----------
    starship_result : StarshipResult
        The Starship prediction to visualise.
    width : int
        SVG viewport width in pixels (default 800).
    height : int
        SVG viewport height in pixels (default 80).

    Returns
    -------
    str
        An SVG element as a raw XML string.
    """
    pad_left = 10
    pad_right = 10
    track_y = 35
    track_h = 20
    arrow_head = 8
    tir_w = 6
    tir_h = 26

    region_start = starship_result.start
    region_end = starship_result.end
    region_len = region_end - region_start
    if region_len <= 0:
        region_len = 1

    draw_w = width - pad_left - pad_right

    def x_of(pos: int) -> float:
        """Map a genomic position to an SVG x coordinate."""
        return pad_left + (pos - region_start) / region_len * draw_w

    lines: List[str] = []
    lines.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" '
        f'style="background:#fff;font-family:Oxygen,Garamond,Georgia,serif;'
        f'font-size:10px;">'
    )

    # Background track (light grey line representing the Starship region)
    lines.append(
        f'  <rect x="{pad_left}" y="{track_y}" '
        f'width="{draw_w}" height="{track_h}" '
        f'rx="3" fill="#f0f0f0" stroke="#ccc" stroke-width="0.5"/>'
    )

    # -- helper: draw a gene arrow ------------------------------------------
    def _gene_arrow(
        g_start: int,
        g_end: int,
        strand: int,
        fill: str,
        label: str = "",
        data_tooltip: str = "",
        css_class: str = "",
    ) -> None:
        x1 = x_of(g_start)
        x2 = x_of(g_end)
        gene_w = abs(x2 - x1)
        if gene_w < 1:
            gene_w = 1
            x2 = x1 + 1

        # Ensure x1 < x2 for drawing
        lx = min(x1, x2)
        rx = max(x1, x2)
        cy = track_y + track_h / 2
        top = track_y + 2
        bot = track_y + track_h - 2
        ah = min(arrow_head, gene_w * 0.4)

        if strand >= 0:
            # Arrow pointing right
            pts = (
                f"{lx},{top} "
                f"{rx - ah},{top} "
                f"{rx},{cy} "
                f"{rx - ah},{bot} "
                f"{lx},{bot}"
            )
        else:
            # Arrow pointing left
            pts = (
                f"{rx},{top} "
                f"{lx + ah},{top} "
                f"{lx},{cy} "
                f"{lx + ah},{bot} "
                f"{rx},{bot}"
            )

        extra_attrs = ""
        if data_tooltip:
            extra_attrs += f' data-tooltip="{data_tooltip}"'
        if css_class:
            extra_attrs += f' class="{css_class}"'

        lines.append(
            f'  <polygon points="{pts}" '
            f'fill="{fill}" stroke="#333" stroke-width="0.4" opacity="0.9"'
            f'{extra_attrs}/>'
        )
        # Label: inside the arrow if it fits, above the arrow otherwise
        if label:
            tx = (lx + rx) / 2
            if gene_w > 30:
                lines.append(
                    f'  <text x="{tx:.1f}" y="{cy + 3:.1f}" '
                    f'text-anchor="middle" fill="#fff" font-size="8">'
                    f'{label}</text>'
                )
            else:
                lines.append(
                    f'  <text x="{tx:.1f}" y="{top - 4:.1f}" '
                    f'text-anchor="middle" fill="{fill}" font-size="8">'
                    f'{label}</text>'
                )

    # Draw cargo genes first (below captain)
    TAG_COLORS = {
        "nlr": "#dc3545", "fre": "#fd7e14", "d37": "#6f42c1",
        "plp": "#17a2b8", "secondary_metabolism": "#28a745",
        "drug_resistance": "#e83e8c", "tyr": "#c77c11", "myb": "#20c997",
    }
    for cargo in starship_result.cargo_genes:
        tag_color = TAG_COLORS.get(cargo.tag) if cargo.tag else None
        fill = tag_color if tag_color else ("#4B77BE" if cargo.strand >= 0 else "#6BAADB")
        strand_str = "+" if cargo.strand >= 0 else "-"
        tag_str = f"|tag: {cargo.tag}" if cargo.tag else ""
        tooltip = (
            f"{cargo.gene_id}|{cargo.product}|"
            f"{cargo.start:,}-{cargo.end:,}|{strand_str}{tag_str}"
        )
        _gene_arrow(cargo.start, cargo.end, cargo.strand, fill,
                     data_tooltip=tooltip, css_class="gene-hover")

    # Draw captain gene on top
    captain = starship_result.captain
    captain_strand = "+" if captain.strand >= 0 else "-"
    captain_tooltip = (
        f"{captain.protein_id}|Captain ({starship_result.captain_family})|"
        f"{captain.start:,}-{captain.end:,}|{captain_strand}|"
        f"e-value: {captain.evalue:.2e}"
    )
    _gene_arrow(
        captain.start,
        captain.end,
        captain.strand,
        "#c77c11",
        label="Captain",
        data_tooltip=captain_tooltip,
        css_class="gene-hover",
    )

    # Draw additional captains (from merged duplicates) in a darker amber
    for extra_captain in starship_result.additional_captains:
        ec_strand = "+" if extra_captain.strand >= 0 else "-"
        ec_tooltip = (
            f"{extra_captain.protein_id}|Additional captain ({extra_captain.hmm_name})|"
            f"{extra_captain.start:,}-{extra_captain.end:,}|{ec_strand}|"
            f"e-value: {extra_captain.evalue:.2e}"
        )
        _gene_arrow(
            extra_captain.start,
            extra_captain.end,
            extra_captain.strand,
            "#a06010",
            label="Captain*",
            data_tooltip=ec_tooltip,
            css_class="gene-hover",
        )

    # DR markers (blue rectangles at element boundaries, same style as TIRs)
    dr_color = "#0b81d5"
    if starship_result.tsd:
        # Left DR
        dr_lx = x_of(region_start)
        dr_ly = track_y + (track_h - tir_h) / 2
        dr_tip_l = f"Left DR|{starship_result.tsd}|{region_start:,}"
        lines.append(
            f'  <rect x="{dr_lx - tir_w / 2:.1f}" y="{dr_ly:.1f}" '
            f'width="{tir_w}" height="{tir_h}" '
            f'fill="{dr_color}" opacity="0.85" rx="1" '
            f'class="gene-hover" data-tooltip="{dr_tip_l}"/>'
        )
        lines.append(
            f'  <text x="{dr_lx:.1f}" y="{dr_ly - 2:.1f}" '
            f'text-anchor="middle" fill="{dr_color}" font-size="8">'
            f'DR</text>'
        )
        # Right DR
        dr_rx = x_of(region_end)
        dr_tip_r = f"Right DR|{starship_result.tsd}|{region_end:,}"
        lines.append(
            f'  <rect x="{dr_rx - tir_w / 2:.1f}" y="{dr_ly:.1f}" '
            f'width="{tir_w}" height="{tir_h}" '
            f'fill="{dr_color}" opacity="0.85" rx="1" '
            f'class="gene-hover" data-tooltip="{dr_tip_r}"/>'
        )
        lines.append(
            f'  <text x="{dr_rx:.1f}" y="{dr_ly - 2:.1f}" '
            f'text-anchor="middle" fill="{dr_color}" font-size="8">'
            f'DR</text>'
        )

    # Scale bar
    scale_y = height - 8
    # Pick a nice round scale bar length
    desired_bp = region_len / 4
    magnitude = 10 ** int(f"{desired_bp:.0e}".split("e+")[-1])
    scale_bp = max(magnitude, 1)
    scale_px = scale_bp / region_len * draw_w
    if scale_px < 20:
        scale_bp = magnitude * 5
        scale_px = scale_bp / region_len * draw_w

    sx = pad_left + 5
    lines.append(
        f'  <line x1="{sx:.1f}" y1="{scale_y}" '
        f'x2="{sx + scale_px:.1f}" y2="{scale_y}" '
        f'stroke="#666" stroke-width="1.5"/>'
    )
    # Tick marks at ends
    for tick_x in (sx, sx + scale_px):
        lines.append(
            f'  <line x1="{tick_x:.1f}" y1="{scale_y - 3}" '
            f'x2="{tick_x:.1f}" y2="{scale_y + 3}" '
            f'stroke="#666" stroke-width="1.5"/>'
        )

    if scale_bp >= 1000:
        scale_label = f"{scale_bp / 1000:.0f} kb"
    else:
        scale_label = f"{scale_bp} bp"

    lines.append(
        f'  <text x="{sx + scale_px / 2:.1f}" y="{scale_y - 4}" '
        f'text-anchor="middle" fill="#666" font-size="9">'
        f'{scale_label}</text>'
    )

    lines.append("</svg>")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Genome map
# ---------------------------------------------------------------------------

FAMILY_COLORS = {
    "Hephaestus": "#c77c11", "Enterprise": "#0b81d5", "Prometheus": "#28a745",
    "Galactica": "#6f42c1", "Tardis": "#dc3545", "Phoenix": "#fd7e14",
    "Arwing": "#20c997", "Voyager": "#e83e8c", "Serenity": "#6c757d",
    "Moya": "#17a2b8", "unclassified": "#999",
}


def generate_genome_map(starkit_run, width=900):
    """Generate an SVG genome map showing contig bars with Starship positions."""
    contig_lengths = starkit_run.genome_stats.get("contig_lengths", {})
    if not contig_lengths:
        return ""

    starships = starkit_run.starships
    if not starships and not contig_lengths:
        return ""

    # Sort contigs by length descending
    total_contigs = len(contig_lengths)
    sorted_contigs = sorted(contig_lengths.items(), key=lambda x: -x[1])

    # Group starships by contig (needed for filtering)
    ships_by_contig = {}
    for s in starkit_run.starships:
        ships_by_contig.setdefault(s.contig_id, []).append(s)

    pad_left = 140
    pad_right = 20
    bar_h = 14
    bar_gap = 6
    top_pad = 10
    subtitle = ""

    if total_contigs > 30:
        # Show only contigs containing Starships
        contigs_with_ships = {ctg for ctg in ships_by_contig}
        sorted_contigs = [
            (ctg, length) for ctg, length in sorted_contigs
            if ctg in contigs_with_ships
        ]
        subtitle = (
            f"Showing {len(sorted_contigs)} of {total_contigs:,} contigs "
            f"(those containing predicted Starships)"
        )
    else:
        sorted_contigs = sorted_contigs[:30]

    if not sorted_contigs:
        return ""

    height = top_pad + len(sorted_contigs) * (bar_h + bar_gap) + 30
    if subtitle:
        height += 16
    draw_w = width - pad_left - pad_right
    genome_max = sorted_contigs[0][1] if sorted_contigs else 1

    def x_of(pos, ctg_len):
        scale = draw_w * (ctg_len / genome_max)
        return pad_left + (pos / ctg_len) * scale if ctg_len > 0 else pad_left

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" '
        f'style="background:#fff;font-family:Oxygen,Garamond,Georgia,serif;font-size:10px;">'
    ]

    # Add subtitle if filtered
    if subtitle:
        lines.append(
            f'  <text x="{width / 2}" y="{top_pad + 4}" '
            f'text-anchor="middle" fill="#888" font-size="9" font-style="italic">'
            f'{subtitle}</text>'
        )
        top_pad += 16

    for i, (ctg_id, ctg_len) in enumerate(sorted_contigs):
        y = top_pad + i * (bar_h + bar_gap)
        bar_w = draw_w * (ctg_len / genome_max)

        # Contig label
        label = ctg_id if len(ctg_id) <= 18 else ctg_id[:15] + "..."
        lines.append(
            f'  <text x="{pad_left - 6}" y="{y + bar_h / 2 + 3}" '
            f'text-anchor="end" fill="#444" font-size="9">{label}</text>'
        )

        # Contig bar
        lines.append(
            f'  <rect x="{pad_left}" y="{y}" width="{bar_w:.1f}" height="{bar_h}" '
            f'rx="2" fill="#eee" stroke="#ccc" stroke-width="0.5"/>'
        )

        # Starships on this contig
        for s in ships_by_contig.get(ctg_id, []):
            sx = x_of(s.start, ctg_len)
            ex = x_of(s.end, ctg_len)
            sw = max(2, ex - sx)
            color = FAMILY_COLORS.get(s.captain_family, "#999")
            tip = (f"{s.starship_id}|{s.captain_family} ({s.size:,}bp)|"
                   f"{s.contig_id}:{s.start:,}-{s.end:,}")
            lines.append(
                f'  <rect x="{sx:.1f}" y="{y + 1}" width="{sw:.1f}" height="{bar_h - 2}" '
                f'rx="1" fill="{color}" opacity="0.85" '
                f'class="gene-hover" data-tooltip="{tip}"/>'
            )

    # Scale bar
    scale_y = height - 12
    scale_bp = 10 ** int(f"{genome_max / 4:.0e}".split("e+")[-1])
    scale_px = scale_bp / genome_max * draw_w
    if scale_px < 30:
        scale_bp *= 5
        scale_px = scale_bp / genome_max * draw_w
    sx = pad_left
    lines.append(
        f'  <line x1="{sx}" y1="{scale_y}" x2="{sx + scale_px:.1f}" y2="{scale_y}" '
        f'stroke="#666" stroke-width="1.5"/>'
    )
    for tx in (sx, sx + scale_px):
        lines.append(
            f'  <line x1="{tx:.1f}" y1="{scale_y - 3}" x2="{tx:.1f}" y2="{scale_y + 3}" '
            f'stroke="#666" stroke-width="1.5"/>'
        )
    scale_label = f"{scale_bp / 1e6:.1f} Mb" if scale_bp >= 1e6 else f"{scale_bp / 1000:.0f} kb"
    lines.append(
        f'  <text x="{sx + scale_px / 2:.1f}" y="{scale_y - 5}" '
        f'text-anchor="middle" fill="#666" font-size="9">{scale_label}</text>'
    )

    lines.append("</svg>")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Jinja2 HTML template (inline string)
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>StarKIT Report</title>
<style>
/* ---- Base ---- */
*, *::before, *::after { box-sizing: border-box; }
body {
    font-family: Oxygen, Garamond, Georgia, serif;
    color: #333;
    background: #fafafa;
    margin: 0;
    padding: 0;
    line-height: 1.5;
}
a { color: #0b81d5; text-decoration: none; }
a:hover { text-decoration: underline; }

/* ---- Layout ---- */
.container {
    max-width: 1100px;
    margin: 0 auto;
    padding: 20px 30px 60px;
}

/* ---- Header ---- */
.header {
    background: #444;
    color: #fff;
    padding: 24px 30px;
    margin-bottom: 30px;
    border-radius: 6px;
}
.header pre.banner {
    font-family: monospace;
    font-size: 13px;
    line-height: 1.25;
    margin: 0 0 14px;
    white-space: pre;
    color: #ddd;
}
.header h1 {
    margin: 0 0 6px;
    font-size: 22px;
    font-weight: 600;
}
.header .meta {
    display: flex;
    flex-wrap: wrap;
    gap: 8px 28px;
    font-size: 14px;
    color: #ccc;
}
.header .meta strong { color: #fff; }

/* ---- Stats cards ---- */
.stats-row {
    display: flex;
    flex-wrap: wrap;
    gap: 14px;
    margin-bottom: 28px;
}
.stat-card {
    flex: 1 1 150px;
    background: #fff;
    border: 1px solid #ddd;
    border-radius: 6px;
    padding: 14px 18px;
    text-align: center;
}
.stat-card .value {
    font-size: 26px;
    font-weight: 700;
    color: #444;
}
.stat-card .label {
    font-size: 12px;
    color: #888;
    text-transform: uppercase;
    letter-spacing: 0.04em;
}

/* ---- Evidence badges ---- */
.badge {
    display: inline-block;
    padding: 2px 10px;
    border-radius: 10px;
    font-size: 12px;
    font-weight: 600;
    color: #fff;
    text-transform: uppercase;
}
.badge-high   { background: #0b81d5; }
.badge-medium { background: #c77c11; }
.badge-low    { background: #999; }
.badge-known  { background: #6c757d; }
.badge-new    { background: #28a745; }
.badge-putative-novel { background: #e83e8c; }

/* ---- Summary table ---- */
.summary-section h2 {
    color: #444;
    margin-bottom: 10px;
}
table.summary {
    width: 100%;
    border-collapse: collapse;
    background: #fff;
    border: 1px solid #ddd;
    border-radius: 6px;
    overflow: hidden;
    margin-bottom: 34px;
    font-size: 14px;
}
table.summary thead th {
    background: #EEE;
    color: #444;
    padding: 10px 12px;
    text-align: left;
    font-weight: 600;
    cursor: pointer;
    user-select: none;
    white-space: nowrap;
    border-bottom: 2px solid #ddd;
}
table.summary thead th:hover {
    background: #e2e2e2;
}
table.summary thead th .sort-arrow {
    font-size: 10px;
    margin-left: 4px;
    color: #aaa;
}
table.summary tbody tr {
    border-bottom: 1px solid #eee;
}
table.summary tbody tr:hover {
    background: #f7f9fc;
}
table.summary tbody td {
    padding: 8px 12px;
    vertical-align: middle;
}
table.summary tbody td.num {
    text-align: right;
    font-variant-numeric: tabular-nums;
}

/* ---- Detail sections ---- */
.detail-section {
    background: #fff;
    border: 1px solid #ddd;
    border-radius: 6px;
    margin-bottom: 18px;
    overflow: hidden;
}
.detail-header {
    background: #EEE;
    padding: 12px 18px;
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: space-between;
    user-select: none;
}
.detail-header:hover {
    background: #e4e4e4;
}
.detail-header h3 {
    margin: 0;
    font-size: 16px;
    color: #444;
}
.detail-header .toggle-icon {
    font-size: 18px;
    color: #888;
    transition: transform 0.2s;
}
.detail-header.open .toggle-icon {
    transform: rotate(90deg);
}
.detail-body {
    display: none;
    padding: 18px;
}
.detail-body.open {
    display: block;
}
.detail-body .svg-wrap {
    overflow-x: auto;
    margin-bottom: 16px;
}

/* Detail tables */
table.detail-table {
    width: 100%;
    border-collapse: collapse;
    margin-bottom: 16px;
    font-size: 13px;
}
table.detail-table th {
    background: #EEE;
    text-align: left;
    padding: 6px 10px;
    font-weight: 600;
    color: #444;
    border-bottom: 1px solid #ddd;
}
table.detail-table td {
    padding: 6px 10px;
    border-bottom: 1px solid #eee;
    vertical-align: top;
}
table.detail-table tr:last-child td {
    border-bottom: none;
}
table.detail-table td.num {
    text-align: right;
    font-variant-numeric: tabular-nums;
}
.mono { font-family: monospace; font-size: 12px; word-break: break-all; }

/* ---- Genome map ---- */
.genome-map { margin-bottom: 28px; }
.genome-map h2 { color: #444; margin-bottom: 10px; }

/* ---- Gene tooltip ---- */
.gene-hover { cursor: pointer; }
.gene-hover:hover { opacity: 0.7 !important; filter: brightness(1.1); }
#gene-tooltip {
    display: none;
    position: fixed;
    background: #333;
    color: #fff;
    padding: 8px 12px;
    border-radius: 5px;
    font-size: 12px;
    line-height: 1.5;
    max-width: 340px;
    pointer-events: none;
    z-index: 1000;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
    font-family: Oxygen, Garamond, Georgia, serif;
}
#gene-tooltip .tt-label {
    font-weight: 600;
    color: #c77c11;
    margin-right: 4px;
}
#gene-tooltip .tt-row {
    white-space: nowrap;
}
#gene-tooltip .tt-product {
    font-style: italic;
    color: #ccc;
}

/* ---- Parameters ---- */
.params-section {
    margin-bottom: 28px;
}
.params-section h2 {
    color: #444;
    margin-bottom: 10px;
}
table.params {
    border-collapse: collapse;
    font-size: 14px;
    background: #fff;
    border: 1px solid #ddd;
    border-radius: 6px;
    overflow: hidden;
}
table.params th {
    background: #EEE;
    text-align: left;
    padding: 6px 14px;
    font-weight: 600;
    color: #444;
    border-bottom: 1px solid #ddd;
}
table.params td {
    padding: 6px 14px;
    border-bottom: 1px solid #eee;
}

/* ---- Footer ---- */
.footer {
    text-align: center;
    font-size: 12px;
    color: #aaa;
    margin-top: 40px;
    padding-top: 16px;
    border-top: 1px solid #eee;
}
</style>
</head>
<body>
<div class="container">

<!-- ===== Header ===== -->
<div class="header">
<pre class="banner">   _____ _             _  _______ _______
  / ____| |           | |/ /_   _|__   __|
 | (___ | |_ __ _ _ __| ' /  | |    | |
  \___ \| __/ _` | '__|  <   | |    | |
  ____) | || (_| | |  | . \ _| |_   | |
 |_____/ \__\__,_|_|  |_|\_\_____|  |_|</pre>
<h1>StarKIT Report</h1>
<div class="meta">
    <span><strong>Input:</strong> {{ input_file }}</span>
    <span><strong>Version:</strong> {{ version }}</span>
    <span><strong>Contigs:</strong> {{ genome_stats.contig_count | default(0) }}</span>
    <span><strong>Genome size:</strong> {{ "{:,}".format(genome_stats.total_length | default(0)) }} bp</span>
    <span><strong>GC:</strong> {{ "%.1f" | format((genome_stats.gc_content | default(0)) * 100) }}%</span>
    {% if runtime %}<span><strong>Runtime:</strong> {{ "%.1f" | format(runtime) }}s</span>{% endif %}
</div>
</div>

<!-- ===== Stats cards ===== -->
<div class="stats-row">
    <div class="stat-card">
        <div class="value">{{ starships | length }}</div>
        <div class="label">Starships found</div>
    </div>
    <div class="stat-card">
        <div class="value">{{ counts.high }}</div>
        <div class="label"><span class="badge badge-high">High</span> evidence</div>
    </div>
    <div class="stat-card">
        <div class="value">{{ counts.medium }}</div>
        <div class="label"><span class="badge badge-medium">Medium</span> evidence</div>
    </div>
    <div class="stat-card">
        <div class="value">{{ counts.low }}</div>
        <div class="label"><span class="badge badge-low">Low</span> evidence</div>
    </div>
    <div class="stat-card">
        <div class="value">{{ counts.novel }}</div>
        <div class="label"><span class="badge badge-new">New</span> to Starbase</div>
    </div>
</div>

<!-- ===== Genome Map ===== -->
{% if genome_map_svg %}
<div class="genome-map">
<h2>Genome Map</h2>
<div style="overflow-x:auto;">
{{ genome_map_svg }}
</div>
</div>
{% endif %}

<!-- ===== Parameters ===== -->
{% if parameters %}
<div class="params-section">
<h2>Run Parameters</h2>
<table class="params">
{% for key, val in parameters.items() %}
<tr><th>{{ key }}</th><td>{{ val }}</td></tr>
{% endfor %}
</table>
</div>
{% endif %}

<!-- ===== Summary table ===== -->
<div class="summary-section">
<h2>Summary</h2>
<table class="summary" id="summaryTable">
<thead>
<tr>
    <th data-type="str">ID <span class="sort-arrow"></span></th>
    <th data-type="str">Contig <span class="sort-arrow"></span></th>
    <th data-type="num">Start <span class="sort-arrow"></span></th>
    <th data-type="num">End <span class="sort-arrow"></span></th>
    <th data-type="num">Size (bp) <span class="sort-arrow"></span></th>
    <th data-type="str">Captain <span class="sort-arrow"></span></th>
    <th data-type="str">Family <span class="sort-arrow"></span></th>
    <th data-type="str">Evidence <span class="sort-arrow"></span></th>
    <th data-type="num">Confidence <span class="sort-arrow"></span></th>
    <th data-type="str">Novelty <span class="sort-arrow"></span></th>
    <th data-type="num">Cargo <span class="sort-arrow"></span></th>
</tr>
</thead>
<tbody>
{% for s in starships %}
<tr>
    <td><a href="#detail-{{ loop.index }}">{{ s.starship_id }}</a></td>
    <td>{{ s.contig_id }}</td>
    <td class="num">{{ "{:,}".format(s.start) }}</td>
    <td class="num">{{ "{:,}".format(s.end) }}</td>
    <td class="num">{{ "{:,}".format(s.size) }}</td>
    <td>{{ s.captain.protein_id }}</td>
    <td>{{ s.captain_family }}</td>
    <td>
        <span class="badge badge-{{ s.evidence_level.value }}">{{ s.evidence_level.value }}</span>
    </td>
    <td class="num">{{ "%.3f" | format(s.confidence_score) }}</td>
    <td>
        {% if s.homology_identity >= 0.80 and s.homology_coverage >= 0.50 %}
        <span class="badge badge-known">Known</span>
        {% else %}
        <span class="badge badge-new">New</span>
        {% endif %}
    </td>
    <td class="num">{{ s.cargo_genes | length }}</td>
</tr>
{% endfor %}
</tbody>
</table>
</div>

<!-- ===== Detail sections ===== -->
{% for s in starships %}
<div class="detail-section" id="detail-{{ loop.index }}">
<div class="detail-header" onclick="toggleDetail(this)">
    <h3>
        {{ s.starship_id }}
        &mdash;
        <span class="badge badge-{{ s.evidence_level.value }}">{{ s.evidence_level.value }}</span>
        {% if s.classification_status == "putative_novel" %}
        <span class="badge badge-putative-novel">putative novel</span>
        {% endif %}
        {{ s.captain_family }}
        ({{ "{:,}".format(s.size) }} bp)
    </h3>
    <span class="toggle-icon">&#9654;</span>
</div>
<div class="detail-body">

<!-- SVG diagram -->
<div class="svg-wrap">
{{ svg_diagrams[loop.index0] }}
</div>

<!-- Starship details -->
<table class="detail-table">
<tr><th>Starship ID</th><td>{{ s.starship_id }}</td></tr>
<tr><th>Contig</th><td>{{ s.contig_id }}</td></tr>
<tr><th>Coordinates</th><td>{{ "{:,}".format(s.start) }} &ndash; {{ "{:,}".format(s.end) }} ({{ "{:,}".format(s.size) }} bp)</td></tr>
<tr><th>Captain protein</th><td>{{ s.captain.protein_id }}</td></tr>
<tr><th>Captain e-value</th><td class="mono">{{ "%.2e" | format(s.captain.evalue) }}</td></tr>
<tr><th>Captain HMM</th><td>{{ s.captain.hmm_name }}</td></tr>
<tr><th>Family</th><td>{{ s.captain_family }}</td></tr>
<tr><th>Family score</th><td>{{ "%.1f" | format(s.family_score) }}
    {% if s.classification_status == "putative_novel" %}
    &mdash; <span class="badge badge-putative-novel">putative novel family</span> (score below threshold)
    {% elif s.classification_status == "unclassified" %}
    &mdash; unclassified (no family HMM hit)
    {% endif %}
</td></tr>
<tr>
    <th>DR sequence</th>
    <td class="mono">{{ s.tsd if s.tsd else "not detected" }}</td>
</tr>
<tr>
    <th>Left TIR</th>
    <td class="mono">{% if s.tir_left %}{{ s.tir_left.sequence }} ({{ "{:,}".format(s.tir_left.start) }}&ndash;{{ "{:,}".format(s.tir_left.end) }}){% else %}not detected{% endif %}</td>
</tr>
<tr>
    <th>Right TIR</th>
    <td class="mono">{% if s.tir_right %}{{ s.tir_right.sequence }} ({{ "{:,}".format(s.tir_right.start) }}&ndash;{{ "{:,}".format(s.tir_right.end) }}){% else %}not detected{% endif %}</td>
</tr>
<tr><th>Confidence score</th><td>{{ "%.4f" | format(s.confidence_score) }}</td></tr>
<tr><th>Boundary method</th><td>{{ s.boundary_method }}</td></tr>
<tr>
    <th>Starbase novelty</th>
    <td>
        {% if s.homology_identity >= 0.80 and s.homology_coverage >= 0.50 %}
        <span class="badge badge-known">Known</span> &mdash; matches Starbase reference ({{ "%.0f" | format(s.homology_identity * 100) }}% identity, {{ "%.0f" | format(s.homology_coverage * 100) }}% coverage)
        {% else %}
        <span class="badge badge-new">New</span> &mdash; no significant match in Starbase
        {% endif %}
    </td>
</tr>
<tr><th>Truncated</th><td>{{ "Yes" if s.truncated else "No" }}</td></tr>
<tr><th>Captain orientation flag</th><td>{{ "Yes" if s.captain_orientation_flag else "No" }}</td></tr>
<tr><th>Captain truncated flag</th><td>{{ "Yes" if s.captain_truncated_flag else "No" }}</td></tr>
{% if s.nested_in %}
<tr><th>Nested inside</th><td>{{ s.nested_in }}</td></tr>
{% endif %}
{% if s.adjacent_to %}
<tr><th>Adjacent to</th><td>{{ s.adjacent_to }}</td></tr>
{% endif %}
{% if s.additional_captains %}
<tr><th>Additional captains</th><td>{{ s.additional_captains | map(attribute='protein_id') | join(', ') }}</td></tr>
{% endif %}
</table>

<!-- Cargo genes -->
{% if s.cargo_genes %}
<h4 style="color:#444;margin:16px 0 8px;">Cargo Genes ({{ s.cargo_genes | length }})</h4>
<table class="detail-table">
<thead>
<tr>
    <th>Gene ID</th>
    <th>Product</th>
    <th>Start</th>
    <th>End</th>
    <th>Strand</th>
    <th>Tag</th>
</tr>
</thead>
<tbody>
{% for cg in s.cargo_genes %}
<tr>
    <td class="mono">{{ cg.gene_id }}</td>
    <td>{{ cg.product }}</td>
    <td class="num">{{ "{:,}".format(cg.start) }}</td>
    <td class="num">{{ "{:,}".format(cg.end) }}</td>
    <td>{{ "+" if cg.strand >= 0 else "-" }}</td>
    <td>{{ cg.tag if cg.tag else "-" }}</td>
</tr>
{% endfor %}
</tbody>
</table>
{% else %}
<p style="color:#888;font-size:13px;">No cargo genes annotated.</p>
{% endif %}

</div>
</div>
{% endfor %}

<div class="footer">
    Generated by StarKIT {{ version }}{% if runtime %} in {{ "%.1f" | format(runtime) }}s{% endif %}
</div>

</div><!-- /.container -->

<!-- Tooltip element -->
<div id="gene-tooltip"></div>

<script>
/* ---- Table sorting ---- */
(function() {
    var table = document.getElementById("summaryTable");
    if (!table) return;
    var headers = table.querySelectorAll("thead th");
    var tbody = table.querySelector("tbody");
    var sortCol = -1;
    var sortAsc = true;

    headers.forEach(function(th, idx) {
        th.addEventListener("click", function() {
            if (sortCol === idx) {
                sortAsc = !sortAsc;
            } else {
                sortCol = idx;
                sortAsc = true;
            }

            // Update arrows
            headers.forEach(function(h) {
                h.querySelector(".sort-arrow").textContent = "";
            });
            th.querySelector(".sort-arrow").textContent = sortAsc ? " \u25B2" : " \u25BC";

            var rows = Array.from(tbody.querySelectorAll("tr"));
            var isNum = th.getAttribute("data-type") === "num";

            rows.sort(function(a, b) {
                var aText = a.children[idx].textContent.trim();
                var bText = b.children[idx].textContent.trim();
                if (isNum) {
                    var aVal = parseFloat(aText.replace(/,/g, "")) || 0;
                    var bVal = parseFloat(bText.replace(/,/g, "")) || 0;
                    return sortAsc ? aVal - bVal : bVal - aVal;
                } else {
                    return sortAsc
                        ? aText.localeCompare(bText)
                        : bText.localeCompare(aText);
                }
            });

            rows.forEach(function(row) { tbody.appendChild(row); });
        });
    });
})();

/* ---- Expand / collapse ---- */
function toggleDetail(header) {
    var body = header.nextElementSibling;
    var isOpen = body.classList.contains("open");
    if (isOpen) {
        body.classList.remove("open");
        header.classList.remove("open");
    } else {
        body.classList.add("open");
        header.classList.add("open");
    }
}

/* ---- Gene tooltips ---- */
(function() {
    var tooltip = document.getElementById("gene-tooltip");
    if (!tooltip) return;

    function formatTooltip(data) {
        // data format: "geneID|product|coords|strand" or with extra "|e-value: X"
        var parts = data.split("|");
        var html = "";
        if (parts[0]) html += '<div class="tt-row"><span class="tt-label">ID:</span> ' + parts[0] + '</div>';
        if (parts[1]) html += '<div class="tt-row"><span class="tt-label">Product:</span> <span class="tt-product">' + parts[1] + '</span></div>';
        if (parts[2]) html += '<div class="tt-row"><span class="tt-label">Coordinates:</span> ' + parts[2] + '</div>';
        if (parts[3]) html += '<div class="tt-row"><span class="tt-label">Strand:</span> ' + parts[3] + '</div>';
        if (parts[4]) html += '<div class="tt-row"><span class="tt-label">' + parts[4] + '</span></div>';
        return html;
    }

    document.addEventListener("mouseover", function(e) {
        var el = e.target.closest(".gene-hover");
        if (!el) return;
        var data = el.getAttribute("data-tooltip");
        if (!data) return;
        tooltip.innerHTML = formatTooltip(data);
        tooltip.style.display = "block";
    });

    document.addEventListener("mousemove", function(e) {
        if (tooltip.style.display === "block") {
            var x = e.clientX + 14;
            var y = e.clientY + 14;
            // Keep tooltip on screen
            var tw = tooltip.offsetWidth;
            var th = tooltip.offsetHeight;
            if (x + tw > window.innerWidth - 10) x = e.clientX - tw - 10;
            if (y + th > window.innerHeight - 10) y = e.clientY - th - 10;
            tooltip.style.left = x + "px";
            tooltip.style.top = y + "px";
        }
    });

    document.addEventListener("mouseout", function(e) {
        var el = e.target.closest(".gene-hover");
        if (el) {
            tooltip.style.display = "none";
        }
    });
})();
</script>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def generate_report(
    starkit_run: StarKITRun,
    output_prefix: str,
    runtime: float = None,
) -> str:
    """Generate a self-contained HTML report for a StarKIT run.

    Parameters
    ----------
    starkit_run : StarKITRun
        The completed analysis run containing genome stats, parameters,
        and Starship results.
    output_prefix : str
        File path prefix. The report is written to ``{output_prefix}.html``.

    Returns
    -------
    str
        The path to the written HTML file.
    """
    starships: List[StarshipResult] = starkit_run.starships

    # Compute evidence-level counts and novelty
    counts = {"high": 0, "medium": 0, "low": 0, "novel": 0}
    for s in starships:
        level_key = s.evidence_level.value
        counts[level_key] = counts.get(level_key, 0) + 1
        if not (s.homology_identity >= 0.80 and s.homology_coverage >= 0.50):
            counts["novel"] += 1

    # Pre-render SVG diagrams
    svg_diagrams = [generate_svg_diagram(s) for s in starships]
    genome_map_svg = generate_genome_map(starkit_run)

    # Render the template
    template = jinja2.Template(HTML_TEMPLATE)
    html = template.render(
        input_file=os.path.basename(starkit_run.input_file),
        version=starkit_run.version,
        genome_stats=starkit_run.genome_stats,
        parameters=starkit_run.parameters,
        starships=starships,
        counts=counts,
        svg_diagrams=svg_diagrams,
        genome_map_svg=genome_map_svg,
        runtime=runtime,
    )

    output_path = f"{output_prefix}.html"
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(html)

    logger.info(f"HTML report written to {output_path}")
    return output_path
