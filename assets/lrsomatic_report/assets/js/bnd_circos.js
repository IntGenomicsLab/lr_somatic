// Breakend circos drawn in the browser from window.BND_DATA (see R/circos_bnd.R) so sectors re-lay-out per filter; 1000x1000 user space scaled by viewBox
(function () {
  "use strict";

  var SVG_NS = "http://www.w3.org/2000/svg";

  // Geometry in the 1000x1000 user space, ordered outward; the ring stays inside the box to leave room for labels
  var CX = 500, CY = 500;
  var R_LINK  = 314;   // arcs terminate here, just inside the ideogram
  var R_IDEO  = 322;   // ideogram ring, inner edge
  var R_IDEO2 = 344;   // ideogram ring, outer edge
  var R_TICK  = 352;   // gene connector, inner end
  var R_LABEL = 378;   // gene labels sit on this radius
  var R_CHROM = 300;   // chromosome names, inside the ring

  var GAP_DEG      = 2;    // between adjacent sectors
  var GAP_LAST_DEG = 6;    // after the last, so the ring has a visible seam
  var START_DEG    = 90;   // 12 o'clock, matching the genome-wide circos
  var MIN_LABEL_SEP_DEG = 3.4;  // greedy de-overlap target for gene labels

  // Giemsa stains, as circlize draws them.
  var STAIN = {
    gneg: "#f5f2ec", gpos25: "#d5cfc4", gpos50: "#b3aa9a", gpos75: "#8d8271",
    gpos100: "#6b6153", gpos: "#6b6153", acen: "#b8593f", gvar: "#9d94c4",
    stalk: "#8fa2b8"
  };

  function el(name, attrs) {
    var n = document.createElementNS(SVG_NS, name);
    for (var k in attrs) if (attrs[k] !== null && attrs[k] !== undefined) {
      n.setAttribute(k, attrs[k]);
    }
    return n;
  }

  // --- Layout ------------------------------------------------------------

  // Angular extent per visible chromosome, proportional to length after the inter-sector gaps
  function layout(chroms, lenOf) {
    var total = 0, i;
    for (i = 0; i < chroms.length; i++) total += lenOf(chroms[i]);
    var gaps = GAP_DEG * Math.max(chroms.length - 1, 0) + GAP_LAST_DEG;
    var usable = 360 - gaps;
    if (total <= 0 || usable <= 0) return {};

    var out = {}, at = START_DEG;
    for (i = 0; i < chroms.length; i++) {
      var span = usable * (lenOf(chroms[i]) / total);
      // Angles run clockwise from 12 o'clock, which is how circlize lays these out.
      out[chroms[i]] = { start: at, span: span, len: lenOf(chroms[i]) };
      at -= span + GAP_DEG;
    }
    return out;
  }

  function angleOf(sectors, chrom, pos) {
    var s = sectors[chrom];
    if (!s) return null;
    var f = s.len > 0 ? Math.min(Math.max(pos / s.len, 0), 1) : 0;
    return s.start - f * s.span;
  }

  function pt(deg, r) {
    var rad = deg * Math.PI / 180;
    return [CX + r * Math.cos(rad), CY - r * Math.sin(rad)];
  }

  // Annular sector between two angles, as a filled path.
  function ringPath(a0, a1, r0, r1) {
    var p0 = pt(a0, r1), p1 = pt(a1, r1), p2 = pt(a1, r0), p3 = pt(a0, r0);
    var large = Math.abs(a1 - a0) > 180 ? 1 : 0;
    // sweep 1 = clockwise on the outer edge, because angles decrease as we go round.
    return "M" + p0[0] + "," + p0[1] +
           "A" + r1 + "," + r1 + " 0 " + large + " 1 " + p1[0] + "," + p1[1] +
           "L" + p2[0] + "," + p2[1] +
           "A" + r0 + "," + r0 + " 0 " + large + " 0 " + p3[0] + "," + p3[1] + "Z";
  }

  // Push labels apart until none is closer than MIN_LABEL_SEP_DEG (forward then backward pass); connectors show the displacement
  function deoverlap(items) {
    if (items.length < 2) return items;
    items.sort(function (a, b) { return b.angle - a.angle; });
    var i;
    for (i = 1; i < items.length; i++) {
      var minA = items[i - 1].angle - MIN_LABEL_SEP_DEG;
      if (items[i].angle > minA) items[i].angle = minA;
    }
    for (i = items.length - 2; i >= 0; i--) {
      var maxA = items[i + 1].angle + MIN_LABEL_SEP_DEG;
      if (items[i].angle < maxA) items[i].angle = maxA;
    }
    return items;
  }

  // --- Drawing -----------------------------------------------------------

  function draw(host, D, visibleIds, visibleGenes, selected) {
    while (host.firstChild) host.removeChild(host.firstChild);

    var lenOf = {}, i;
    for (i = 0; i < D.chromosomes.length; i++) lenOf[D.chromosomes[i]] = D.lengths[i];

    // Only the links this filter leaves visible; null means "no filter yet".
    var links = D.links.filter(function (l) {
      return visibleIds === null || visibleIds.has(l.id);
    });

    var touched = {};
    links.forEach(function (l) { touched[l.chromA] = 1; touched[l.chromB] = 1; });
    var chroms = D.chromosomes.filter(function (c) { return touched[c]; });

    if (!chroms.length) {
      var note = document.createElement("p");
      note.className = "bnd-circos-empty";
      note.textContent = "No breakends match the current filter.";
      host.appendChild(note);
      return { arcs: 0, genes: 0, chroms: 0 };
    }

    var sectors = layout(chroms, function (c) { return lenOf[c] || 0; });

    var svg = el("svg", {
      viewBox: "0 0 1000 1000",
      preserveAspectRatio: "xMidYMid meet",
      role: "img",
      "aria-label": "Breakend circos over " + chroms.length + " chromosomes"
    });

    var gBands  = el("g", { class: "bnd-bands" });
    var gBodies = el("g", { class: "bnd-bodies" });
    var gLines  = el("g", { class: "bnd-lines" });
    var gLabels = el("g", { class: "bnd-labels" });
    var gChrom  = el("g", { class: "bnd-chroms" });
    var gLinks  = el("g", { class: "bnd-links" });
    // Arcs last so they sit above the ring; labels above those again.
    [gLinks, gBands, gBodies, gLines, gChrom, gLabels].forEach(function (g) {
      svg.appendChild(g);
    });

    // Ideogram bands.
    D.cytobands.forEach(function (b) {
      if (!sectors[b.chrom]) return;
      var a0 = angleOf(sectors, b.chrom, b.start);
      var a1 = angleOf(sectors, b.chrom, b.end);
      if (a0 === null || a1 === null || a0 === a1) return;
      gBands.appendChild(el("path", {
        d: ringPath(a0, a1, R_IDEO, R_IDEO2),
        class: "bnd-band",
        fill: STAIN[b.stain] || STAIN.gneg,
        "data-chrom": b.chrom
      }));
    });

    // Sector outline, so a chromosome with sparse banding still reads as one block.
    chroms.forEach(function (c) {
      var s = sectors[c];
      gBands.appendChild(el("path", {
        d: ringPath(s.start, s.start - s.span, R_IDEO, R_IDEO2),
        class: "bnd-sector",
        "data-chrom": c
      }));
    });

    // Chromosome names, inside the ring, upright.
    chroms.forEach(function (c) {
      var s = sectors[c];
      var p = pt(s.start - s.span / 2, R_CHROM);
      var t = el("text", {
        x: p[0], y: p[1], class: "bnd-chrom-label", "data-chrom": c,
        "text-anchor": "middle", "dominant-baseline": "middle"
      });
      t.textContent = c.replace(/^chr/, "");
      gChrom.appendChild(t);
    });

    // Gene track: only the genes the visible rows actually name.
    var genes = D.genes.filter(function (g) {
      return sectors[g.chrom] && visibleGenes.has(g.gene);
    });

    // Gene bodies get a floor of ~half a degree: they are markers, not spans to scale
    var placed = deoverlap(genes.map(function (g) {
      var mid = (g.start + g.end) / 2;
      return {
        gene: g.gene, panels: g.panels, chrom: g.chrom,
        at: angleOf(sectors, g.chrom, mid),
        angle: angleOf(sectors, g.chrom, mid),
        a0: angleOf(sectors, g.chrom, g.start),
        a1: angleOf(sectors, g.chrom, g.end)
      };
    }).filter(function (g) { return g.at !== null; }));

    placed.forEach(function (g) {
      var half = Math.max(Math.abs(g.a0 - g.a1) / 2, 0.25);
      gBodies.appendChild(el("path", {
        d: ringPath(g.at + half, g.at - half, R_IDEO, R_IDEO2),
        class: "bnd-gene-body",
        "data-gene": g.gene, "data-panels": g.panels
      }));

      // Connector from the locus out to wherever de-overlapping moved the label.
      var p0 = pt(g.at, R_TICK), p1 = pt(g.angle, R_LABEL - 6);
      gLines.appendChild(el("polyline", {
        points: p0[0] + "," + p0[1] + " " + p1[0] + "," + p1[1],
        class: "bnd-gene-line", "data-gene": g.gene
      }));

      // Labels read outward on the right half and inward on the left, so each needs its own rotate()
      var flip = Math.cos(g.angle * Math.PI / 180) < 0;
      var lp = pt(g.angle, R_LABEL);
      var rot = flip ? (180 - g.angle) : -g.angle;
      var t = el("text", {
        x: lp[0], y: lp[1],
        class: "bnd-gene-label",
        "data-gene": g.gene, "data-panels": g.panels,
        "text-anchor": flip ? "end" : "start",
        "dominant-baseline": "middle",
        transform: "rotate(" + rot + " " + lp[0] + " " + lp[1] + ")"
      });
      t.textContent = g.gene;
      gLabels.appendChild(t);
    });

    // Arcs: quadratic Bezier with the control point at the centre (circlize-style bundling)
    links.forEach(function (l) {
      var aA = angleOf(sectors, l.chromA, l.posA);
      var aB = angleOf(sectors, l.chromB, l.posB);
      if (aA === null || aB === null) return;
      var pA = pt(aA, R_LINK), pB = pt(aB, R_LINK);
      var picked = selected.has(l.id);
      var arc = el("path", {
        d: "M" + pA[0] + "," + pA[1] + "Q" + CX + "," + CY + " " + pB[0] + "," + pB[1],
        class: "bnd-link" + (picked ? " is-selected" : ""),
        "data-svid": l.id,
        "data-svclass": l.svclass,
        "data-chrom-a": l.chromA,
        "data-chrom-b": l.chromB
      });
      // Selected arcs go last within the group so they draw over their neighbours.
      if (picked) gLinks.appendChild(arc); else gLinks.insertBefore(arc, gLinks.firstChild);
    });

    host.appendChild(svg);
    return { arcs: links.length, genes: placed.length, chroms: chroms.length };
  }

  // --- Public surface ----------------------------------------------------

  function normalise(raw) {
    return {
      chromosomes: raw.chromosomes || [],
      lengths: raw.lengths || [],
      cytobands: (raw.cytobands || []).map(function (r) {
        return { chrom: r[0], start: r[1], end: r[2], stain: r[3] };
      }),
      links: (raw.links || []).map(function (r) {
        return { id: r[0], svclass: r[1], chromA: r[2], posA: r[3],
                 chromB: r[4], posB: r[5] };
      }),
      genes: (raw.genes || []).map(function (r) {
        return { chrom: r[0], start: r[1], end: r[2], gene: r[3], panels: r[4] };
      })
    };
  }

  var DATA = null;

  window.bndCircos = {
    // Draw into `host` for a filter state (visibleIds/visibleGenes/selected Sets, null = everything); returns {arcs, genes, chroms}
    render: function (host, visibleIds, visibleGenes, selected) {
      if (!host) return { arcs: 0, genes: 0, chroms: 0 };
      if (!DATA) {
        if (!window.BND_DATA) return { arcs: 0, genes: 0, chroms: 0 };
        DATA = normalise(window.BND_DATA);
      }
      return draw(host, DATA, visibleIds, visibleGenes || new Set(),
                  selected || new Set());
    }
  };
})();
