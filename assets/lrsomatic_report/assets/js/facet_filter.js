// Checkbox-dropdown ("tickbox") column filters for the SNV and SV tables. R publishes values and column maps (window.*_FACETS / window.*_COLS); the menu lives in <body> as position:fixed because DT's scroll containers clip it; no ids or document-delegated handlers (DT clones the header every draw); ticked state lives in a closure; this predicate must stay LAST in ext.search so the live counts measure the rows every other filter leaves
(function () {
  "use strict";

  // consequence is ~30 terms; impact is 4. Below this a find box is just clutter.
  var SEARCH_THRESHOLD = 12;
  // One draw over 167k rows runs two ext.search predicates per row, so coalesce rapid ticks.
  var DRAW_DEBOUNCE_MS = 120;
  var LABEL_MAX = 16;

  // Tables that get facets, keyed on the globals they publish (same test as the gene-panel predicate)
  var TABLES = [
    { key: "snv",
      elem:      function () { return window.snvTableElem; },
      cols:      function () { return window.SNV_COLS; },
      facets:    function () { return window.SNV_FACETS; },
      requested: function () { return window.SNV_FACET_COLS; } },
    { key: "sv",
      elem:      function () { return window.svTableElem; },
      cols:      function () { return window.SV_COLS; },
      facets:    function () { return window.SV_FACETS; },
      requested: function () { return window.SV_FACET_COLS; } }
  ];

  // key -> { colName -> Set of ticked tokens }; null is the "(none)" bucket, an absent column is unconstrained
  var STATE    = { snv: {}, sv: {} };
  var ATTACHED = {};
  var WIDGETS  = [];      // every built widget, for the "clear all" control
  var PENDING  = {};      // per-table draw debounce timers
  var OPEN     = null;    // the one open menu

  // Counting machinery keyed by table; FCOLS = facet columns resolved in *_COLS, everything else parallel to it
  var FCOLS   = {};       // key -> [colName]
  var VOCAB   = {};       // key -> colName -> { strs: [token|null], map: Map, noneId: int }
  var CIDX    = {};       // key -> colName -> { starts: Int32Array, ids: Int32Array }
  var COUNTS  = {};       // key -> colName -> [rows per vocabulary id]
  var TOUCHED = {};       // key -> has the predicate run since the last flush?

  // Per-pass scratch indexed by FCOLS position, reused across rows
  var SC_SRC = [], SC_LO = [], SC_HI = [], SC_TMP = [];

  function tableFor(node) {
    for (var i = 0; i < TABLES.length; i++) if (TABLES[i].elem() === node) return TABLES[i];
    return null;
  }

  function esc(s) {
    return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;")
                    .replace(/>/g, "&gt;").replace(/"/g, "&quot;");
  }

  // ---- matching: separators come from R so tokens mirror how the cell was built; NA and "" fall in the null bucket ----
  function tokens(v, sep) {
    if (v === null || v === undefined) return [null];
    var s = String(v).trim();
    if (s === "") return [null];
    if (!sep) return [s];
    var out = s.split(sep).map(function (x) { return x.trim(); })
                          .filter(function (x) { return x.length > 0; });
    return out.length ? out : [null];
  }

  // ---- vocabulary: tokens interned to integer ids, seeded in payload order, auto-extending for the fallback split path ----
  function addTok(v, tok) {
    var id = v.strs.length;
    v.strs.push(tok);
    if (tok === null) v.noneId = id; else v.map.set(tok, id);
    return id;
  }

  function idOf(v, tok) {
    if (tok === null || tok === undefined) {
      return v.noneId >= 0 ? v.noneId : addTok(v, null);
    }
    var id = v.map.get(tok);
    return id === undefined ? addTok(v, tok) : id;
  }

  function newVocab(def) {
    var v = { strs: [], map: new Map(), noneId: -1 };
    var vals = (def && def.values) || [];
    for (var i = 0; i < vals.length; i++) addTok(v, vals[i][0]);
    return v;
  }

  // Distinct ids of one cell into `out`, deduplicated per row like js_facet_defs()
  function idsInto(voc, toks, out) {
    var n = 0;
    for (var i = 0; i < toks.length; i++) {
      var id = idOf(voc, toks[i]), dup = false;
      for (var k = 0; k < n; k++) if (out[k] === id) { dup = true; break; }
      if (!dup) out[n++] = id;
    }
    return n;
  }

  function resetCounts(key) {
    var cnt = COUNTS[key], voc = VOCAB[key];
    if (!cnt) return;
    Object.keys(cnt).forEach(function (nm) {
      var a = cnt[nm], n = voc[nm].strs.length;
      a.length = n;
      for (var i = 0; i < n; i++) a[i] = 0;
    });
  }

  // Set up counting structures for a table on first sight (from the predicate or attach())
  function ensureTable(t) {
    if (FCOLS[t.key]) return FCOLS[t.key].length > 0;
    var defs = t.facets(), cols = t.cols();
    if (!defs || !cols) return false;             // not published yet — retry on the next row
    var names = [], voc = {}, cnt = {};
    Object.keys(defs).forEach(function (nm) {
      if (cols[nm] === undefined) return;         // renamed column: attach() warns loudly
      names.push(nm);
      voc[nm] = newVocab(defs[nm]);
      cnt[nm] = [];
    });
    FCOLS[t.key]  = names;
    VOCAB[t.key]  = voc;
    COUNTS[t.key] = cnt;
    resetCounts(t.key);
    return names.length > 0;
  }

  function bump(arr, src, lo, hi) {
    for (var k = lo; k < hi; k++) {
      var id = src[k];
      arr[id] = (arr[id] || 0) + 1;
    }
  }

  // Tokenise every row of every facet column once into a CSR index keyed by DataTables row index (~6 MB for 167k rows)
  function buildIndex(t, api) {
    var names = FCOLS[t.key];
    if (!names || !names.length) return;
    var defs = t.facets(), cols = t.cols();
    var idxs = api.rows().indexes().toArray();
    var data = api.rows().data();
    var maxRow = -1, k;
    for (k = 0; k < idxs.length; k++) if (idxs[k] > maxRow) maxRow = idxs[k];
    if (maxRow < 0) return;

    var byIdx = new Array(maxRow + 1);
    for (k = 0; k < idxs.length; k++) byIdx[idxs[k]] = data[k];

    var idx = {}, scratch = [];
    for (var c = 0; c < names.length; c++) {
      var nm = names[c], voc = VOCAB[t.key][nm], sep = defs[nm].sep, col = cols[nm];
      var starts = new Int32Array(maxRow + 2);
      var buf = new Int32Array(Math.max(16, (maxRow + 1) * 2));
      var total = 0;
      for (var r = 0; r <= maxRow; r++) {
        var row = byIdx[r];
        var n = (row === undefined) ? 0 : idsInto(voc, tokens(row[col], sep), scratch);
        if (total + n > buf.length) {
          var grown = new Int32Array(Math.max(buf.length * 2, total + n));
          grown.set(buf);
          buf = grown;
        }
        for (var q = 0; q < n; q++) buf[total + q] = scratch[q];
        total += n;
        starts[r + 1] = total;
      }
      idx[nm] = { starts: starts, ids: buf.subarray(0, total) };
    }
    CIDX[t.key] = idx;
  }

  // One ext.search predicate for both tables: OR within a column, AND across (DataTables ANDs it with the panel and text filters); a token rule rather than a regex avoids escaping bugs. Also accumulates the exclude-own-column counts; `counter` restarts at 0 per pass
  function facetPredicate(settings, searchData, index, rowData, counter) {
    var t = tableFor(settings.nTable);
    if (!t) return true;                                  // another DT table in the report
    if (!ensureTable(t)) return true;                     // no facet column on this table

    var key = t.key;
    if (counter === 0 || !TOUCHED[key]) resetCounts(key);
    TOUCHED[key] = true;

    var names = FCOLS[key], sel = STATE[key], voc = VOCAB[key], cnt = COUNTS[key];
    var cols = t.cols(), defs = t.facets(), idx = CIDX[key];
    // Use rowData (typed), not searchData (rendered): the SV table renders locus/size differently
    var row = rowData || searchData;
    var fails = 0, failed = -1, i;

    for (i = 0; i < names.length; i++) {
      var nm = names[i], chosen = sel[nm];
      var hasSel = !!(chosen && chosen.size);
      // After one failed column, only that column can still be counted for this row
      if (!hasSel && fails > 0) { SC_LO[i] = SC_HI[i] = 0; continue; }

      var src, lo, hi;
      var ci = idx && idx[nm];
      if (ci && index >= 0 && index + 1 < ci.starts.length) {
        src = ci.ids; lo = ci.starts[index]; hi = ci.starts[index + 1];
      } else {
        var out = SC_TMP[i] || (SC_TMP[i] = []);
        hi  = idsInto(voc[nm], tokens(row[cols[nm]], defs[nm].sep), out);
        src = out; lo = 0;
      }
      SC_SRC[i] = src; SC_LO[i] = lo; SC_HI[i] = hi;

      if (hasSel) {
        var hit = false, strs = voc[nm].strs;
        for (var k = lo; k < hi; k++) {
          if (chosen.has(strs[src[k]])) { hit = true; break; }
        }
        if (!hit) {
          fails++; failed = i;
          if (fails > 1) break;                           // counts for nobody; stop reading
        }
      }
    }

    if (fails > 1) return false;
    if (fails === 1) {                                    // AND across columns
      bump(cnt[names[failed]], SC_SRC[failed], SC_LO[failed], SC_HI[failed]);
      return false;
    }
    for (i = 0; i < names.length; i++) {
      bump(cnt[names[i]], SC_SRC[i], SC_LO[i], SC_HI[i]);
    }
    return true;
  }

  $.fn.dataTable.ext.search.push(facetPredicate);

  function scheduleDraw(w) {
    clearTimeout(PENDING[w.key]);
    PENDING[w.key] = setTimeout(function () {
      // resetPaging default: a filter change returns to the top
      w.api.draw();
    }, DRAW_DEBOUNCE_MS);
  }

  // ---- the open menu -------------------------------------------------------
  var EDGE = 8;    // keep this much clear of every viewport edge

  function closeOpen() {
    if (!OPEN) return;
    var w = OPEN;
    w.$menu.prop("hidden", true);
    w.$btn.attr("aria-expanded", "false");
    OPEN = null;
  }

  // Fixed menu placed from the button's viewport rect: left-aligned, flipped right/up and clamped when there is no room
  function placeMenu(w) {
    var b    = w.$btn[0].getBoundingClientRect();
    var menu = w.$menu[0];
    var mw   = menu.offsetWidth, mh = menu.offsetHeight;

    var left = b.left;
    if (left + mw > window.innerWidth - EDGE) left = b.right - mw;
    left = Math.max(EDGE, Math.min(left, window.innerWidth - EDGE - mw));

    var top = b.bottom + 4;
    if (top + mh > window.innerHeight - EDGE) {
      var above = b.top - 4 - mh;
      top = above >= EDGE ? above
                          : Math.max(EDGE, window.innerHeight - EDGE - mh);
    }

    menu.style.left = Math.round(left) + "px";
    menu.style.top  = Math.round(top) + "px";
  }

  function openMenu(w) {
    if (OPEN === w) { closeOpen(); return; }
    closeOpen();
    w.$menu.prop("hidden", false);           // measurable only once it is not display:none
    w.$btn.attr("aria-expanded", "true");
    OPEN = w;
    placeMenu(w);
    // preventScroll: letting the browser reveal the focused field would jump the table
    if (w.$find.length) w.$find[0].focus({ preventScroll: true });
  }

  // A fixed menu must follow its button: one capturing document scroll listener covers the page and the scroll containers
  function followButton() {
    if (!OPEN) return;
    var w = OPEN;
    var b = w.$btn[0].getBoundingClientRect();
    var clip = w.clip ? w.clip.getBoundingClientRect() : null;
    // Column scrolled out from under the menu: close it
    if ((b.width === 0 && b.height === 0) ||
        b.bottom < 0 || b.top > window.innerHeight ||
        (clip && (b.right <= clip.left || b.left >= clip.right))) {
      closeOpen();
      return;
    }
    placeMenu(w);
  }

  // ---- widget: button in the header cell, menu appended to <body> (see the header note) ----
  function buttonMarkup(name) {
    return '<div class="facet">' +
             '<button type="button" class="facet__btn" aria-expanded="false"' +
             ' aria-haspopup="true" title="Filter ' + esc(name) + ' by value">' +
               '<span class="facet__label">All</span>' +
               '<span class="facet__caret" aria-hidden="true">&#9662;</span>' +
             '</button>' +
           '</div>';
  }

  function menuMarkup(name, def) {
    var withFind = def.values.length > SEARCH_THRESHOLD;
    var opts = def.values.map(function (v) {
      var tok = v[0], n = v[1];
      var isNone = (tok === null);
      var label = isNone ? "(none)" : tok;
      return '<label class="facet__opt' + (isNone ? " facet__opt--none" : "") + '"' +
             ' title="' + esc(label) + '">' +
             '<input type="checkbox"' +
             (isNone ? ' data-none="1"' : ' data-tok="' + esc(tok) + '"') + '>' +
             '<span class="facet__v">' + esc(label) + '</span>' +
             '<span class="facet__n">' + Number(n).toLocaleString() + '</span>' +
             '</label>';
    }).join("");

    return '<div class="facet__menu" role="group" aria-label="' + esc(name) +
           ' values" hidden>' +
             '<div class="facet__tools">' +
               (withFind ? '<input type="text" class="facet__find"' +
                           ' placeholder="Find value…" autocomplete="off">' : "") +
               '<button type="button" class="facet__all" title="Tick all shown">All' +
               '</button>' +
               '<button type="button" class="facet__none" title="Untick all shown">None' +
               '</button>' +
             '</div>' +
             '<div class="facet__list">' + opts + '</div>' +
             '<p class="facet__foot">counts follow the other active filters</p>' +
           '</div>';
  }

  function tokenOf($cb) {
    return $cb.data("none") ? null : String($cb.attr("data-tok"));
  }

  // A column with nothing ticked is absent from STATE, keeping the per-row fast path cheap
  var NO_SELECTION = new Set();

  function readSel(w) {
    return STATE[w.key][w.name] || NO_SELECTION;
  }

  function ensureSel(w) {
    var sel = STATE[w.key][w.name];
    if (!sel) { sel = new Set(); STATE[w.key][w.name] = sel; }
    return sel;
  }

  function pruneSel(w) {
    var sel = STATE[w.key][w.name];
    if (sel && sel.size === 0) delete STATE[w.key][w.name];
  }

  function trunc(s) {
    return s.length > LABEL_MAX ? s.slice(0, LABEL_MAX - 1) + "…" : s;
  }

  function refreshLabel(w) {
    var sel = readSel(w), n = sel.size, text, title = "";
    if (n === 0) {
      text = "All";
    } else if (n === 1) {
      var only = sel.values().next().value;
      var lab  = (only === null) ? "(none)" : only;
      text  = trunc(lab);
      title = lab;
    } else {
      text  = n + " selected";
      title = Array.from(sel).map(function (v) { return v === null ? "(none)" : v; })
                             .join(", ");
    }
    w.$label.text(text);
    w.$btn.attr("title", title || ("Filter " + w.name + " by value"));
    w.$container.toggleClass("facet--active", n > 0);
  }

  // Paint the counts from the last filtering pass, touching only changed numbers
  function flushCounts(key) {
    var cnt = COUNTS[key];
    if (!cnt) return;
    for (var i = 0; i < WIDGETS.length; i++) {
      var w = WIDGETS[i];
      if (w.key !== key) continue;
      var arr = cnt[w.name];
      if (!arr) continue;
      for (var j = 0; j < w.optIds.length; j++) {
        var n = arr[w.optIds[j]] || 0;
        if (w.optLast[j] === n) continue;
        w.optLast[j] = n;
        w.numNodes[j].textContent = n.toLocaleString();
        // The baseline count moves into the tooltip
        w.optNodes[j].title = w.optLabel[j] + " — " + n.toLocaleString() + " of " +
                              w.optBase[j].toLocaleString() + " rows";
        w.optNodes[j].classList.toggle("facet__opt--empty", n === 0);
      }
    }
  }

  function anyActive() {
    return TABLES.some(function (t) {
      return Object.keys(STATE[t.key]).some(function (c) {
        return STATE[t.key][c] && STATE[t.key][c].size > 0;
      });
    });
  }

  // Reset lives in the report's control bar (a ticked filter can be scrolled out of sight); guarded for older templates
  function syncClearButton() {
    var el = document.getElementById("facet-clear");
    if (el) el.hidden = !anyActive();
  }

  function clearAll() {
    WIDGETS.forEach(function (w) {
      delete STATE[w.key][w.name];
      w.$menu.find('input[type="checkbox"]').prop("checked", false);
      refreshLabel(w);
    });
    syncClearButton();
    TABLES.forEach(function (t) {
      var node = t.elem();
      if (node && $.fn.dataTable.isDataTable(node)) $(node).DataTable().draw();
    });
  }

  function build(t, api, $td, name, def) {
    var $wrap  = $td.children("div").first();      // DT's div.form-group.has-feedback
    var $input = $wrap.children("input");

    // Hidden, not removed: DT holds a reference to this input; we never write to it
    $input.hide().attr("tabindex", "-1").attr("aria-hidden", "true");
    $wrap.children("span.glyphicon").hide();

    var $container = $(buttonMarkup(name));
    $wrap.append($container);                      // inside $wrap, see note 1 at the top
    var $menu = $(menuMarkup(name, def)).appendTo(document.body);

    // What clips the button, so an open menu can tell its column scrolled away
    var $cont = $(api.table().container());
    var $clip = $cont.find(".dataTables_scrollHead").first();

    var w = {
      key: t.key, name: name, api: api,
      clip: ($clip.length ? $clip : $cont)[0],
      $container: $container,
      $btn:   $container.find(".facet__btn"),
      $label: $container.find(".facet__label"),
      $menu:  $menu,
      $find:  $menu.find(".facet__find"),
      $list:  $menu.find(".facet__list")
    };

    // Option rows in payload order, paired with vocabulary ids; cached raw nodes for flushCounts()
    var voc = VOCAB[t.key][name];
    w.optNodes = w.$list.children(".facet__opt").toArray();
    w.numNodes = w.optNodes.map(function (el) { return el.querySelector(".facet__n"); });
    w.optIds   = def.values.map(function (v) { return idOf(voc, v[0]); });
    w.optBase  = def.values.map(function (v) { return Number(v[1]) || 0; });
    w.optLabel = def.values.map(function (v) { return v[0] === null ? "(none)" : v[0]; });
    w.optLast  = def.values.map(function () { return -1; });   // -1: nothing painted yet

    WIDGETS.push(w);

    w.$btn.on("click", function (e) {
      e.stopPropagation();
      openMenu(w);
    });

    // Delegated within our own node, never from `document` (the header clone would match too)
    w.$menu.on("change", 'input[type="checkbox"]', function () {
      var sel = ensureSel(w), tok = tokenOf($(this));
      if (this.checked) sel.add(tok); else sel.delete(tok);
      pruneSel(w);
      refreshLabel(w);
      syncClearButton();
      scheduleDraw(w);
    });

    w.$menu.find(".facet__all, .facet__none").on("click", function () {
      var on  = $(this).hasClass("facet__all");
      var sel = ensureSel(w);
      // Act only on the values currently shown by the find box
      w.$list.children(".facet__opt").not(".facet__opt--hidden")
        .find('input[type="checkbox"]').each(function () {
          this.checked = on;
          var tok = tokenOf($(this));
          if (on) sel.add(tok); else sel.delete(tok);
        });
      pruneSel(w);
      refreshLabel(w);
      syncClearButton();
      scheduleDraw(w);
    });

    if (w.$find.length) {
      w.$find.on("input", function () {
        var q = this.value.trim().toLowerCase();
        w.$list.children(".facet__opt").each(function () {
          var txt = $(this).find(".facet__v").text().toLowerCase();
          $(this).toggleClass("facet__opt--hidden", q.length > 0 && txt.indexOf(q) < 0);
        });
      });
    }

    // Bound to both nodes: Escape must work from the find box and the button
    w.$container.add(w.$menu).on("keydown", function (e) {
      if (e.key === "Escape" || e.keyCode === 27) {
        closeOpen();
        w.$btn.trigger("focus");
      }
    });

    refreshLabel(w);
    return w;
  }

  // ---- attach: re-push the predicate so it filters last ----
  function moveToEnd() {
    var list = $.fn.dataTable.ext.search;
    var at = list.indexOf(facetPredicate);
    if (at >= 0 && at !== list.length - 1) {
      list.splice(at, 1);
      list.push(facetPredicate);
    }
  }

  function attach(t) {
    if (ATTACHED[t.key]) return;
    var node = t.elem();
    if (!node || !$.fn.dataTable.isDataTable(node)) return;
    var defs = t.facets(), cols = t.cols();
    if (!defs || !cols) return;

    var api    = $(node).DataTable();
    var $head  = $(api.table().header());          // the live thead, wherever DT moved it
    var $rows  = $head.children("tr");
    if ($rows.length < 2) return;                  // no filter row: nothing to replace
    ATTACHED[t.key] = true;
    ensureTable(t);
    moveToEnd();

    var $label  = $rows.first().children("th,td");
    var $filter = $rows.last().children("td");
    var visIdx  = api.columns(":visible").indexes().toArray();

    if ($filter.length !== visIdx.length) {
      console.warn("facet_filter: " + t.key + " filter row has " + $filter.length +
                   " cells for " + visIdx.length +
                   " visible columns — leaving the text filters alone.");
      return;
    }

    // R drops columns it cannot offer a dropdown for and knitr swallows the message, so log it here (info: expected for `callers`/`caller`)
    var requested = t.requested() || [];
    var missing = requested.filter(function (n) { return !(n in defs); });
    if (missing.length) {
      console.info("facet_filter: no value list for " + t.key + " column(s) " +
                   missing.join(", ") + " (fewer than " +
                   "2 distinct values, or not in the table) — they keep their text filter.");
    }

    Object.keys(defs).forEach(function (name) {
      var idx = cols[name];
      if (idx === undefined) {
        console.warn("facet_filter: '" + name + "' is not in " + t.key.toUpperCase() +
                     "_COLS — was the column renamed?");
        return;
      }
      var pos = visIdx.indexOf(idx);
      if (pos < 0) {
        console.warn("facet_filter: '" + name + "' is a hidden column.");
        return;
      }
      // Cross-check: the header label above the chosen cell must be the expected column
      var got = $.trim($label.eq(pos).text());
      if (got !== name) {
        console.warn("facet_filter: header mismatch for '" + name + "' at visible column " +
                     pos + " (found '" + got + "') — skipping.");
        return;
      }
      var $td = $filter.eq(pos);
      if ($td.attr("data-type") !== "character") {
        // Factor/logical columns already carry DT's selectize control
        console.warn("facet_filter: '" + name + "' has data-type '" +
                     $td.attr("data-type") + "' — leaving DT's own filter in place.");
        return;
      }
      build(t, api, $td, name, defs[name]);
    });

    syncClearButton();

    // Speed only: the predicate splits live until the index exists; a failure here must not break the filter
    try {
      buildIndex(t, api);
    } catch (err) {
      console.warn("facet_filter: could not index " + t.key +
                   " tokens, falling back to splitting per draw.", err);
    }

    // DataTables fires `search` right after the pass that filled the counts
    $(node).on("search.dt", function (e) {
      if (e.target !== node) return;
      if (!TOUCHED[t.key]) resetCounts(t.key);
      TOUCHED[t.key] = false;
      flushCounts(t.key);
    });

    // R's counts are the unfiltered baseline and --gene-panel can pre-filter, so force one pass now
    api.draw(false);
  }

  // Defer a tick so DT's own filter-row handlers exist; the load sweep covers tables initialised before this script
  $(document).on("init.dt", function (e, settings) {
    var t = tableFor(settings.nTable);
    if (t) setTimeout(function () { attach(t); }, 0);
  });

  $(window).on("load", function () {
    TABLES.forEach(function (t) { if (t.elem()) attach(t); });
    var clear = document.getElementById("facet-clear");
    if (clear) clear.addEventListener("click", clearAll);
    syncClearButton();
  });

  // Single outside-click handler on document, testing containment against the widget's two nodes
  $(document).on("mousedown", function (e) {
    if (!OPEN) return;
    if (OPEN.$container[0].contains(e.target) || OPEN.$menu[0].contains(e.target)) return;
    closeOpen();
  });

  // Re-place the fixed menu on any scroll (capturing)
  document.addEventListener("scroll", followButton, true);
  window.addEventListener("resize", followButton);

  // Ticked state, exposed for the console and tests/js/test-facet-predicate.js
  window.facetFilterState = STATE;

  // Live counts of the last pass as { table: { column: [[token, n], ...] } }
  window.facetFilterCounts = function () {
    var out = {};
    Object.keys(COUNTS).forEach(function (key) {
      var cnt = COUNTS[key], voc = VOCAB[key];
      out[key] = {};
      Object.keys(cnt).forEach(function (nm) {
        // Iterate the vocabulary, not the count array: unmatched values are 0, not gaps
        var strs = voc[nm].strs, pairs = [];
        for (var id = 0; id < strs.length; id++) pairs.push([strs[id], cnt[nm][id] || 0]);
        out[key][nm] = pairs;
      });
    });
    return out;
  };
})();
