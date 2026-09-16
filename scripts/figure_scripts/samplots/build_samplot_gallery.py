#!/usr/bin/env python3
#
# Build a single self-contained HTML page (no external files/CDNs - every
# image is base64-embedded) reviewing every samplot trio PNG samplot_trio.py
# wrote for one family/label (denovo or inherited). Written directly into
# the same folder as the PNGs as index.html, so it opens standalone in any
# browser with no server and no dependency on the PNGs staying alongside it.
#
# Usage:
#   build_samplot_gallery.py --pngdir <dir> --family <FAMILY> --label denovo|inherited [--out <path>]

import argparse
import base64
import glob
import html
import os
import re


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pngdir", required=True)
    p.add_argument("--family", required=True)
    p.add_argument("--label", required=True, choices=["denovo", "inherited"])
    p.add_argument("--out", default=None)
    args = p.parse_args()

    out_path = args.out or os.path.join(args.pngdir, "index.html")

    pattern = re.compile(
        rf"^{re.escape(args.family)}_{re.escape(args.label)}_(.+)_chr([0-9XYM]+)_(\d+)_(\d+)_(DEL|DUP|INS|INV)\.png$"
    )

    items = []
    for path in sorted(glob.glob(os.path.join(args.pngdir, "*.png"))):
        fname = os.path.basename(path)
        m = pattern.match(fname)
        if not m:
            continue
        genes, chrom, start, end, svtype = m.groups()
        items.append({
            "genes": genes.replace("_", "/"),
            "chrom": f"chr{chrom}",
            "start": int(start),
            "end": int(end),
            "svtype": svtype,
            "img": base64.b64encode(open(path, "rb").read()).decode(),
        })

    types_present = sorted(set(it["svtype"] for it in items))

    cards_html = []
    for it in items:
        size = it["end"] - it["start"]
        size_str = f"{size:,} bp" if size > 0 else "breakend"
        genes_esc = html.escape(it["genes"])
        title = f"{genes_esc} &middot; {it['chrom']}:{it['start']:,}&ndash;{it['end']:,} &middot; {it['svtype']}"
        cards_html.append(f'''
    <article class="card" data-svtype="{it['svtype']}" data-search="{html.escape(it['genes'].lower())} {it['chrom']} {it['svtype'].lower()}">
      <header class="card-head">
        <div class="card-gene">{genes_esc}</div>
        <span class="pill pill-{it['svtype']}">{it['svtype']}</span>
      </header>
      <div class="card-coords">{it['chrom']}:{it['start']:,}&ndash;{it['end']:,} <span class="dim">&middot; {size_str}</span></div>
      <button class="thumb-btn" data-title="{title}" aria-label="Open full size plot for {genes_esc}">
        <img class="thumb" src="data:image/png;base64,{it['img']}" alt="Samplot trio plot: {genes_esc} {it['chrom']}:{it['start']}-{it['end']} {it['svtype']}" loading="lazy" />
      </button>
    </article>''')

    filter_buttons = "\n".join(
        f'        <button class="filter-btn" data-filter-type="{t}">{t}</button>'
        for t in types_present
    )

    label_title = "De novo" if args.label == "denovo" else "Inherited"

    if items:
        body_main = f'''
  <div class="grid" id="grid">
    {''.join(cards_html)}
  </div>
  <div class="empty-state" id="empty-state" hidden>No variants match this filter.</div>'''
    else:
        body_main = f'''
  <div class="empty-state">No {label_title.lower()} R135 green-gene candidates survived filtering for family {html.escape(args.family)}.</div>'''

    html_doc = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>Family {html.escape(args.family)} &middot; {label_title} Samplot Review</title>
<style>
:root {{
  --bg: #ffffff;
  --surface: #ffffff;
  --surface-2: #f7f7f7;
  --border: #e0e0e0;
  --text: #1a1a1a;
  --text-muted: #6b6b6b;
  --accent: #2563eb;
  --accent-soft: #eaf1fd;
  --del: #a8412a;
  --del-soft: #f5e5e0;
  --dup: #2f5f97;
  --dup-soft: #e4edf6;
  --ins: #8a6d1f;
  --ins-soft: #f1ecd9;
  --inv: #6b4f8c;
  --inv-soft: #ebe3f1;
}}
@media (prefers-color-scheme: dark) {{
  :root {{
    --bg: #18191b;
    --surface: #1f2123;
    --surface-2: #232527;
    --border: #333537;
    --text: #e7e7e7;
    --text-muted: #9a9a9a;
    --accent: #5b9bf5;
    --accent-soft: #1d2b40;
    --del: #d9856d;
    --del-soft: #33241e;
    --dup: #85aede;
    --dup-soft: #1e2c3c;
    --ins: #cbab5c;
    --ins-soft: #302a1a;
    --inv: #b09bd0;
    --inv-soft: #291f34;
  }}
}}
* {{ box-sizing: border-box; }}
html, body {{
  margin: 0;
  padding: 0;
  background: var(--bg);
  color: var(--text);
  font-family: system-ui, -apple-system, "Segoe UI", Roboto, Arial, sans-serif;
}}
body {{ overflow-x: hidden; }}
header.topbar {{
  position: sticky;
  top: 0;
  z-index: 10;
  background: var(--bg);
  border-bottom: 1px solid var(--border);
}}
.topbar-inner {{
  max-width: 1400px;
  margin: 0 auto;
  padding: 20px 28px 16px;
  display: flex;
  flex-direction: column;
  gap: 14px;
}}
.title-row {{ display: flex; align-items: baseline; gap: 12px; flex-wrap: wrap; }}
h1 {{ font-size: 1.25rem; font-weight: 600; margin: 0; text-wrap: balance; }}
.subtitle {{ color: var(--text-muted); font-size: 0.9rem; line-height: 1.5; margin: 0; max-width: 70ch; }}
.count-badge {{
  font-size: 0.8rem;
  color: var(--text-muted);
  background: var(--surface-2);
  border: 1px solid var(--border);
  padding: 2px 10px;
  border-radius: 4px;
  font-variant-numeric: tabular-nums;
}}
.controls {{ display: flex; gap: 10px; flex-wrap: wrap; align-items: center; }}
.filter-group {{ display: flex; gap: 6px; flex-wrap: wrap; }}
.filter-btn {{
  font-size: 0.82rem;
  border: 1px solid var(--border);
  background: var(--surface);
  color: var(--text-muted);
  padding: 5px 12px;
  border-radius: 5px;
  cursor: pointer;
  font-family: inherit;
}}
.filter-btn:hover {{ border-color: var(--accent); color: var(--text); }}
.filter-btn.active {{ background: var(--accent-soft); border-color: var(--accent); color: var(--accent); }}
.filter-btn:focus-visible, .thumb-btn:focus-visible, .search-input:focus-visible, .lightbox-close:focus-visible {{
  outline: 2px solid var(--accent);
  outline-offset: 2px;
}}
.search-input {{
  font-size: 0.85rem;
  border: 1px solid var(--border);
  background: var(--surface);
  color: var(--text);
  padding: 6px 12px;
  border-radius: 5px;
  min-width: 200px;
  flex: 1 1 200px;
  max-width: 320px;
  font-family: inherit;
}}
.search-input::placeholder {{ color: var(--text-muted); }}
main {{ max-width: 1400px; margin: 0 auto; padding: 24px 28px 64px; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(380px, 1fr)); gap: 18px; }}
.card {{
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: 8px;
  overflow: hidden;
  display: flex;
  flex-direction: column;
}}
.card[hidden] {{ display: none; }}
.card-head {{ display: flex; align-items: center; justify-content: space-between; gap: 10px; padding: 12px 14px 0; }}
.card-gene {{ font-weight: 600; font-size: 0.95rem; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
.pill {{ font-size: 0.72rem; font-weight: 600; letter-spacing: 0.02em; padding: 2px 8px; border-radius: 4px; flex-shrink: 0; }}
.pill-DEL {{ background: var(--del-soft); color: var(--del); }}
.pill-DUP {{ background: var(--dup-soft); color: var(--dup); }}
.pill-INS {{ background: var(--ins-soft); color: var(--ins); }}
.pill-INV {{ background: var(--inv-soft); color: var(--inv); }}
.card-coords {{ font-size: 0.78rem; color: var(--text-muted); padding: 4px 14px 10px; font-variant-numeric: tabular-nums; }}
.card-coords .dim {{ opacity: 0.8; }}
.thumb-btn {{
  border: none;
  background: var(--surface);
  padding: 0;
  cursor: zoom-in;
  display: block;
  width: 100%;
  border-top: 1px solid var(--border);
  border-bottom: 1px solid var(--border);
}}
.thumb {{ display: block; width: 100%; height: auto; }}
.empty-state {{ text-align: center; color: var(--text-muted); padding: 60px 20px; font-size: 0.9rem; }}
.lightbox {{
  position: fixed;
  inset: 0;
  background: rgba(0, 0, 0, 0.75);
  display: none;
  align-items: center;
  justify-content: center;
  padding: 32px;
  z-index: 100;
}}
.lightbox.open {{ display: flex; }}
.lightbox-panel {{ background: var(--surface); border-radius: 8px; max-width: min(1100px, 100%); max-height: 100%; overflow: auto; }}
.lightbox-head {{ display: flex; align-items: center; justify-content: space-between; gap: 12px; padding: 14px 18px; border-bottom: 1px solid var(--border); }}
.lightbox-title {{ font-size: 0.92rem; font-weight: 600; }}
.lightbox-close {{
  border: 1px solid var(--border);
  background: var(--surface);
  color: var(--text-muted);
  width: 30px;
  height: 30px;
  border-radius: 5px;
  cursor: pointer;
  font-size: 1rem;
  line-height: 1;
}}
.lightbox-close:hover {{ color: var(--text); border-color: var(--accent); }}
.lightbox-panel img {{ display: block; width: 100%; height: auto; }}
@media (prefers-reduced-motion: reduce) {{ * {{ transition: none !important; }} }}
</style>
</head>
<body>
<header class="topbar">
  <div class="topbar-inner">
    <div class="title-row">
      <h1>Family {html.escape(args.family)} &middot; {label_title} R135 Candidates</h1>
      <span class="count-badge" id="count-badge">{len(items)} variant{'s' if len(items) != 1 else ''}</span>
    </div>
    <p class="subtitle">Samplot trio evidence (proband / mother / father) for every SV surviving tier-1 filtering, AnnotSV annotation, and R135 green-gene panel matching.</p>
'''

    if items:
        html_doc += f'''    <div class="controls">
      <div class="filter-group" id="type-filters">
        <button class="filter-btn active" data-filter-type="all">All</button>
{filter_buttons}
      </div>
      <input class="search-input" id="search" type="search" placeholder="Filter by gene&hellip;" aria-label="Filter by gene name" />
    </div>
'''

    html_doc += '''  </div>
</header>

<main>'''
    html_doc += body_main
    html_doc += '''
</main>
'''

    if items:
        html_doc += '''
<div class="lightbox" id="lightbox">
  <div class="lightbox-panel">
    <div class="lightbox-head">
      <span class="lightbox-title" id="lightbox-title"></span>
      <button class="lightbox-close" id="lightbox-close" aria-label="Close">&times;</button>
    </div>
    <img id="lightbox-img" src="" alt="" />
  </div>
</div>

<script>
const grid = document.getElementById('grid');
const cards = Array.from(grid.querySelectorAll('.card'));
const emptyState = document.getElementById('empty-state');
const countBadge = document.getElementById('count-badge');
const typeButtons = Array.from(document.querySelectorAll('#type-filters .filter-btn'));
const searchInput = document.getElementById('search');

let activeType = 'all';

function applyFilters() {
  const q = searchInput.value.trim().toLowerCase();
  let visible = 0;
  for (const card of cards) {
    const matchesType = activeType === 'all' || card.dataset.svtype === activeType;
    const matchesSearch = !q || card.dataset.search.includes(q);
    const show = matchesType && matchesSearch;
    card.hidden = !show;
    if (show) visible++;
  }
  emptyState.hidden = visible !== 0;
  countBadge.textContent = visible === cards.length ? `${cards.length} variants` : `${visible} / ${cards.length} variants`;
}

typeButtons.forEach(btn => {
  btn.addEventListener('click', () => {
    typeButtons.forEach(b => b.classList.remove('active'));
    btn.classList.add('active');
    activeType = btn.dataset.filterType;
    applyFilters();
  });
});

searchInput.addEventListener('input', applyFilters);

const lightbox = document.getElementById('lightbox');
const lightboxImg = document.getElementById('lightbox-img');
const lightboxTitle = document.getElementById('lightbox-title');
const lightboxClose = document.getElementById('lightbox-close');

function openLightbox(btn) {
  const src = btn.querySelector('img').src;
  const title = btn.dataset.title;
  lightboxImg.src = src;
  lightboxImg.alt = title;
  lightboxTitle.textContent = title;
  lightbox.classList.add('open');
}
function closeLightbox() {
  lightbox.classList.remove('open');
  lightboxImg.src = '';
}

grid.querySelectorAll('.thumb-btn').forEach(btn => {
  btn.addEventListener('click', () => openLightbox(btn));
});
lightboxClose.addEventListener('click', closeLightbox);
lightbox.addEventListener('click', (e) => { if (e.target === lightbox) closeLightbox(); });
document.addEventListener('keydown', (e) => { if (e.key === 'Escape') closeLightbox(); });
</script>
'''

    html_doc += '''</body>
</html>
'''

    with open(out_path, "w") as f:
        f.write(html_doc)

    print(f"{args.pngdir}: {len(items)} plots -> {out_path}")


if __name__ == "__main__":
    main()
