/* Full Circle: client for the JSON files written by python/scripts/export_website.py */
"use strict";

const NS = "http://www.w3.org/2000/svg";
const $ = (s, root = document) => root.querySelector(s);
const STATES = { AL: "Alabama", AK: "Alaska", AZ: "Arizona", AR: "Arkansas", CA: "California", CO: "Colorado", CT: "Connecticut", DE: "Delaware", FL: "Florida", GA: "Georgia", HI: "Hawaii", ID: "Idaho", IL: "Illinois", IN: "Indiana", IA: "Iowa", KS: "Kansas", KY: "Kentucky", LA: "Louisiana", ME: "Maine", MD: "Maryland", MA: "Massachusetts", MI: "Michigan", MN: "Minnesota", MS: "Mississippi", MO: "Missouri", MT: "Montana", NE: "Nebraska", NV: "Nevada", NH: "New Hampshire", NJ: "New Jersey", NM: "New Mexico", NY: "New York", NC: "North Carolina", ND: "North Dakota", OH: "Ohio", OK: "Oklahoma", OR: "Oregon", PA: "Pennsylvania", RI: "Rhode Island", SC: "South Carolina", SD: "South Dakota", TN: "Tennessee", TX: "Texas", UT: "Utah", VT: "Vermont", VA: "Virginia", WA: "Washington", WV: "West Virginia", WI: "Wisconsin", WY: "Wyoming" };
const PARTY = { D: "Democrat", R: "Republican", I: "Independent" };
const MONTHS = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"];
const LATEST = 119;
const IDX = (v) => Math.round(100 * v);   // polarization index: 100 x circular variance
// Context events for the trend chart (timing, not cause)
const EVENTS = [
  { c: 104, text: "1994: Republicans win the House for the first time in 40 years.", url: "https://en.wikipedia.org/wiki/1994_United_States_House_of_Representatives_elections" },
  { c: 112, text: "2010: Tea Party wave election.", url: "https://en.wikipedia.org/wiki/2010_United_States_House_of_Representatives_elections" },
  { c: 114, text: "2015: House Freedom Caucus forms; Speaker Boehner resigns.", url: "https://en.wikipedia.org/wiki/Freedom_Caucus" },
  { c: 118, text: "2023: Speaker McCarthy removed, with 8 Republicans joining all Democrats.", url: "https://en.wikipedia.org/wiki/Kevin_McCarthy" },
  { c: 119, text: "2025: Massie (R) and Khanna (D) force the Epstein files vote past leadership; it passes 427 to 1.", url: "https://www.news9.com/national-politics-oklahoma/epstein-discharge-petition-secures-final-signature-needed-to-force-house-vote-on-releasing-files" },
  { c: 119, text: "2025: Record 43-day government shutdown ends.", url: "https://www.nhpr.org/2025-11-12/house-votes-to-end-the-longest-ever-government-shutdown" },
  { c: 119, text: "2026: War with Iran begins; House war powers votes split Republicans.", url: "https://www.npr.org/2026/06/03/nx-s1-5845102/house-iran-war-powers-vote" },
  { c: 119, text: "2026: Massie loses his primary to a Trump-endorsed challenger.", url: "https://www.pbs.org/newshour/politics/trump-backed-gallrein-defeats-rep-thomas-massie-in-gop-primary" },
];

const ord = (n) => n + (n % 100 >= 11 && n % 100 <= 13 ? "th" : ["th", "st", "nd", "rd"][n % 10] || "th");
const years = (c) => `${1787 + 2 * c}–${1789 + 2 * c}`;
const fmtDate = (d) => { const [y, m, dd] = d.split("-").map(Number); return `${MONTHS[m - 1]} ${dd}, ${y}`; };
const esc = (s) => String(s).replace(/[&<>"]/g, (ch) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[ch]));
const fold = (s) => s.normalize("NFD").replace(/[̀-ͯ]/g, "").toLowerCase();
const dist = (st, d) => (d === 0 ? `${st} at-large` : `${st}-${d}`);
const bill = (b) => b.replace(/^(HCONRES|HJRES|HRES|HR|SCONRES|SJRES|SRES|S)(\d+)$/, (_, p, n) => ({ HR: "H.R.", HRES: "H.Res.", HJRES: "H.J.Res.", HCONRES: "H.Con.Res.", S: "S.", SRES: "S.Res.", SJRES: "S.J.Res.", SCONRES: "S.Con.Res." }[p] + " " + n));
function el(tag, attrs, parent) { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
function txt(parent, x, y, s, attrs = {}) { const t = el("text", { x, y, ...attrs }, parent); t.textContent = s; return t; }
function clear(n) { while (n.firstChild) n.firstChild.remove(); }

/* ---------- data ---------- */
const cache = {};
function getJSON(name) {
  if (!cache[name]) cache[name] = fetch(`data/${name}.json`).then((r) => { if (!r.ok) throw new Error(`${name}: ${r.status}`); return r.json(); });
  return cache[name];
}
let INDEX, PEOPLE, BYID;
const houseInfo = (c) => INDEX.houses.find((h) => h.congress === c);
async function house(c) {
  const H = await getJSON(`H${c}`);
  if (!H.byId) {
    H.list = H.members.map((m, k) => ({ k, id: m[0], party: m[1], state: m[2], d: m[3], angle: m[4], rank: m[5], lo: m[6], hi: m[7], votes: m[8], placed: m[9], name: BYID[m[0]].name }));
    H.byId = Object.fromEntries(H.list.map((m) => [m.id, m]));
    H.rollList = H.rolls.map((r, j) => ({ j, roll: r[0], date: r[1], bill: r[2], question: r[3], desc: r[4], result: r[5], yea: r[6], nay: r[7], score: r[8], key: fold(`${r[2]} ${bill(r[2])} ${r[3]} ${r[4]} ${r[1]} roll ${r[0]}`).replace(/\./g, "") }));
  }
  return H;
}
async function votesOf(c) {
  const V = await getJSON(`V${c}`);
  if (!V.byRoll) V.byRoll = Object.fromEntries(V.rolls.map((r, j) => [r, V.votes[j]]));
  return V;
}

/* ---------- tooltip ---------- */
const tip = $("#tip");
function showTip(e, html) { tip.innerHTML = html; tip.hidden = false; const x = Math.min(Math.max(e.clientX, 130), innerWidth - 130); tip.style.left = x + "px"; tip.style.top = e.clientY + "px"; }
function hideTip() { tip.hidden = true; }

/* ---------- the circle ---------- */
// Members at their true angles (0 = top, between the two parties); dots stack in layers on both sides of the ring.
// Dense Houses get smaller dots so the stack stays readable.
function layout(list, R) {
  // members at their true angles (0 = top, between the two parties); where they crowd, dots stack outward
  // like a histogram bent around the circle. Dense Houses get smaller dots so stacks stay short.
  const sorted = list.slice().sort((a, b) => a.angle - b.angle);
  for (const r of [2.7, 2.5, 2.3, 2.1, 1.9, 1.75, 1.6]) {
    const step = 2 * r + 0.45, rad = (L) => R + L * step;
    const last = []; let maxL = 0;
    const pos = sorted.map((m) => {
      let L = 0;
      while (last[L] !== undefined && (m.angle - last[L]) * rad(L) < 2 * r + 0.5) L++;
      last[L] = m.angle; maxL = Math.max(maxL, L);
      return { m, L };
    });
    if (maxL * step <= 92 || r === 1.6) return { r, pts: pos.map(({ m, L }) => ({ m, rr: rad(L), x: 300 + rad(L) * Math.sin(m.angle), y: 0 - rad(L) * Math.cos(m.angle) })) };
  }
}
let svgUid = 0;
function arrowDefs(svg) {
  const id = `arw${++svgUid}`;
  const defs = el("defs", {}, svg);
  const mk = el("marker", { id, viewBox: "0 0 10 10", refX: 8, refY: 5, markerWidth: 7, markerHeight: 7, orient: "auto-start-reverse" }, defs);
  el("path", { d: "M0,0 L10,5 L0,10 z", fill: "var(--ink2)" }, mk);
  return id;
}
// curved arrow + label along the circle, from angle a0 to a1 (radians, 0 = top)
function curvedLabel(svg, cy, rad, a0, a1, label, mid) {
  const p = (a) => [300 + rad * Math.sin(a), cy - rad * Math.cos(a)];
  const [x0, y0] = p(a0), [x1, y1] = p(a1);
  const sweep = a1 > a0 ? 1 : 0;
  el("path", { d: `M ${x0} ${y0} A ${rad} ${rad} 0 0 ${sweep} ${x1} ${y1}`, fill: "none", stroke: "var(--ink2)", "stroke-width": 1.6, "marker-end": `url(#${mid})` }, svg);
  const id = `lp${++svgUid}`, lr = rad + 10;
  const q = (a) => [300 + lr * Math.sin(a), cy - lr * Math.cos(a)];
  const [u0, v0] = q(Math.min(a0, a1)), [u1, v1] = q(Math.max(a0, a1));
  el("path", { id, d: `M ${u0} ${v0} A ${lr} ${lr} 0 0 1 ${u1} ${v1}`, fill: "none" }, svg);
  const tx = el("text", { style: "font-weight:600" }, svg);
  const tp = el("textPath", { href: `#${id}`, startOffset: "50%", "text-anchor": "middle" }, tx);
  tp.textContent = label;
}
function drawCircle(svg, H, opts) {
  clear(svg);
  const cy = opts.cy, R = 168;
  if (!H._lay) H._lay = layout(H.list, R);
  const { r, pts } = H._lay;
  const info = houseInfo(H.congress);
  el("circle", { cx: 300, cy, r: R, fill: "none", stroke: "var(--line)", "stroke-width": 1, "stroke-dasharray": "2 5" }, svg);
  const mid = arrowDefs(svg);
  curvedLabel(svg, cy, 272, -0.25, -0.95, "More liberal", mid);
  curvedLabel(svg, cy, 272, 0.25, 0.95, "More conservative", mid);
  const minA = Math.min(...H.list.map((m) => m.angle)), maxA = Math.max(...H.list.map((m) => m.angle));
  // the back of the circle, where the two ends would meet
  const a0 = maxA + 0.12, a1 = 2 * Math.PI + minA - 0.12;
  if (a1 > a0) {
    const p = (a) => `${300 + R * Math.sin(a)} ${cy - R * Math.cos(a)}`;
    el("path", { d: `M ${p(a0)} A ${R} ${R} 0 ${a1 - a0 > Math.PI ? 1 : 0} 1 ${p(a1)}`, fill: "none", stroke: "var(--ink3)", "stroke-width": 1.6, "stroke-dasharray": "5 6" }, svg);
  }
  if (info.chi_0[0] >= 0.2) {
    txt(svg, 300, cy + R + 30, "Opposite ends, same vote", { "text-anchor": "middle", style: "font-weight:700; fill: var(--ink)" });
    txt(svg, 300, cy + R + 50, "Here the two ends join against the middle", { "text-anchor": "middle" });
  } else {
    txt(svg, 300, cy + R + 30, "The two ends rarely meet here", { "text-anchor": "middle", style: "font-weight:600" });
    txt(svg, 300, cy + R + 50, "this House still behaves like a line", { "text-anchor": "middle" });
  }

  const vs = opts.votes;
  const g = el("g", {}, svg);
  for (const p of pts) {
    const a = { cx: p.x.toFixed(1), cy: (cy + p.y).toFixed(1), r };
    if (!vs) a.fill = `var(--${p.m.party})`;
    else {
      const ch = vs[p.m.k];
      if (ch === "Y") a.fill = "var(--yes)";
      else if (ch === "N") a.fill = "var(--no)";
      else { a.fill = "none"; a.stroke = "var(--ink3)"; a["stroke-width"] = 0.9; a.r = Math.max(r - 0.45, 1); }
    }
    el("circle", a, g);
  }
  const ring = el("circle", { r: 6.5, fill: "none", stroke: "var(--ink)", "stroke-width": 2, visibility: "hidden" }, svg);
  const mark = (p) => { if (!p) { ring.setAttribute("visibility", "hidden"); return; } ring.setAttribute("cx", p.x); ring.setAttribute("cy", cy + p.y); ring.setAttribute("visibility", "visible"); };
  const nearest = (e) => {
    const pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY;
    const q = pt.matrixTransform(svg.getScreenCTM().inverse());
    let best = null, bd = 196;
    for (const p of pts) { const d = (p.x - q.x) ** 2 + (cy + p.y - q.y) ** 2; if (d < bd) { bd = d; best = p; } }
    return best;
  };
  svg.onpointermove = (e) => { const p = nearest(e); mark(p); if (p) showTip(e, opts.tip(p.m)); else hideTip(); };
  svg.onpointerleave = () => { if (!opts.keep || !opts.keep()) mark(null); hideTip(); };
  svg.onclick = (e) => { const p = nearest(e); if (p && opts.click) opts.click(p.m); };
  if (opts.selected) mark(pts.find((p) => p.m.id === opts.selected));
  return { mark: (id) => mark(pts.find((p) => p.m.id === id)) };
}
function memberTip(H, m, vs) {
  const v = vs ? vs[m.k] : null;
  return `<b>${esc(m.name)}</b><br>${PARTY[m.party]}, ${dist(m.state, m.d)}<br>Rank ${m.rank} of ${H.nRank} from the liberal end` +
    (vs ? `<br>Voted: ${v === "Y" ? "Yes" : v === "N" ? "No" : "did not vote (absent or present)"}` : "");
}

/* ---------- vote tallies ---------- */
function groupShares(H, vs, which) {
  const ms = H.list.slice().sort((a, b) => a.rank - b.rank || a.k - b.k), n = ms.length, q = Math.round(n / 10);
  const groups = [["Most liberal 10%", ms.slice(0, q)], ["Democrats in between", ms.slice(q, n - q).filter((m) => m.party === "D")],
    ["Republicans in between", ms.slice(q, n - q).filter((m) => m.party === "R")], ["Most conservative 10%", ms.slice(n - q)]];
  return groups.map(([lab, g]) => {
    const cast = g.filter((m) => vs[m.k] === "Y" || vs[m.k] === "N");
    return [lab, cast.length ? cast.filter((m) => vs[m.k] === which).length / cast.length : 0];
  });
}
function tallyHTML(H, vs, which) {
  return groupShares(H, vs, which).map(([lab, s]) => `<span class="lab">${lab}</span><div class="track"><div class="fill" style="width:${(100 * s).toFixed(0)}%"></div></div><span class="mono">${(100 * s).toFixed(0)}%</span>`).join("");
}
// yes / no split per group, in the vote colors
function splitHTML(H, vs) {
  return groupShares(H, vs, "Y").map(([lab, y]) => `<span class="lab">${lab}</span><div class="track split"><span class="y" style="width:${(100 * y).toFixed(0)}%"></span><span class="n" style="width:${(100 - 100 * y).toFixed(0)}%"></span></div><span class="mono">${(100 * y).toFixed(0)}–${(100 - 100 * y).toFixed(0)}</span>`).join("");
}

/* ---------- trend charts ---------- */
function trendFrame(svg, W, H, yMax, yStep) {
  clear(svg);
  const L = 50, Rm = 14, T = 22, B = 34;
  const x = (c) => L + (c - 99.5) / 20 * (W - L - Rm), y = (v) => T + (1 - v / yMax) * (H - T - B);
  for (const h of INDEX.houses) el("rect", { x: x(h.congress - 0.5), y: T, width: x(1) - x(0), height: H - T - B, fill: h.majority === "D" ? "var(--majD)" : "var(--majR)" }, svg);
  for (let v = 0; v <= yMax + 1e-9; v += yStep) {
    el("line", { x1: L, x2: W - Rm, y1: y(v), y2: y(v), stroke: "var(--line)", "stroke-width": 1 }, svg);
    txt(svg, L - 8, y(v) + 4, String(Math.round(v)), { "text-anchor": "end", class: "tick" });
  }
  for (const c of [100, 104, 108, 112, 116, 119]) txt(svg, x(c), H - B + 22, String(1787 + 2 * c), { "text-anchor": "middle", class: "tick" });
  return { x, y };
}
function series(svg, x, y, key, color, label) {
  const d = INDEX.houses, v = (h, i) => 100 * h[key][i];
  const band = d.map((h) => `${x(h.congress)},${y(v(h, 2))}`).join(" ") + " " + d.slice().reverse().map((h) => `${x(h.congress)},${y(v(h, 1))}`).join(" ");
  el("polygon", { points: band, fill: color, opacity: .16 }, svg);
  el("polyline", { points: d.map((h) => `${x(h.congress)},${y(v(h, 0))}`).join(" "), fill: "none", stroke: color, "stroke-width": 2.2, "stroke-linejoin": "round" }, svg);
  for (const h of d) {
    el("circle", { cx: x(h.congress), cy: y(v(h, 0)), r: 3.6, fill: color, stroke: "var(--surface)", "stroke-width": 1.5 }, svg);
    const hit = el("rect", { x: x(h.congress - 0.5), y: y(v(h, 0)) - 14, width: x(1) - x(0), height: 28, fill: "transparent" }, svg);
    hit.addEventListener("pointermove", (e) => showTip(e, `<b>${ord(h.congress)} Congress</b> (${years(h.congress)})<br>${label}: ${v(h, 0).toFixed(1)} (range ${v(h, 1).toFixed(1)}–${v(h, 2).toFixed(1)})<br>χ = ${h[key][0].toFixed(3)} · ${PARTY[h.majority]} majority${h.congress === LATEST ? "<br>In progress" : ""}`));
    hit.addEventListener("pointerleave", hideTip);
  }
}
function drawTrends() {
  let svg = $("#trend-all");
  let f = trendFrame(svg, 600, 340, 45, 10);
  series(svg, f.x, f.y, "chi_0", "var(--ink)", "Polarization index");
  const top = INDEX.houses.reduce((a, b) => (b.chi_0[0] > a.chi_0[0] ? b : a));
  txt(svg, f.x(top.congress) - 10, f.y(100 * top.chi_0[2]) - 6, "Record high", { class: "ann", "text-anchor": "end" });
  // event markers along the top: one badge per Congress, numbered like the list below
  for (const c of [...new Set(EVENTS.map((e) => e.c))]) {
    const nums = EVENTS.map((e, j) => [e.c, j + 1]).filter(([cc]) => cc === c).map(([, j]) => j);
    const label = nums.length > 1 ? `${nums[0]}–${nums[nums.length - 1]}` : String(nums[0]);
    const w = 10 + 7 * label.length, xx = Math.min(f.x(c), 600 - 14 - w / 2);
    el("line", { x1: f.x(c), x2: f.x(c), y1: 20, y2: f.y(0), stroke: "var(--ink3)", "stroke-dasharray": "2 3", opacity: .6 }, svg);
    el("rect", { x: xx - w / 2, y: 2, width: w, height: 18, rx: 9, fill: "var(--ink)" }, svg);
    txt(svg, xx, 15, label, { "text-anchor": "middle", style: "font: 600 11px var(--f-body); fill: var(--bg)" });
  }
  $("#events").innerHTML = EVENTS.map((ev) => `<li>${esc(ev.text)} <a href="${ev.url}" rel="noopener" target="_blank">source</a></li>`).join("");
  $("#trend-legend").innerHTML = `<span><i class="sq" style="background:var(--majD);box-shadow:inset 0 0 0 1px var(--D)"></i>Democratic majority</span><span><i class="sq" style="background:var(--majR);box-shadow:inset 0 0 0 1px var(--R)"></i>Republican majority</span>`;
  svg = $("#trend-party");
  const pmax = Math.max(...INDEX.houses.map((h) => Math.max(h.chi_D[2], h.chi_R[2]))) * 100;
  f = trendFrame(svg, 600, 300, Math.ceil(pmax / 2) * 2, 2);
  series(svg, f.x, f.y, "chi_D", "var(--D)", "Democrats");
  series(svg, f.x, f.y, "chi_R", "var(--R)", "Republicans");
  // headline ratio, first vs latest Congress
  const h0 = INDEX.houses[0], h1 = INDEX.houses[INDEX.houses.length - 1];
  $("#party-ratio").textContent = `In ${h0.years[0]}, Democrats were about ${Math.round(h0.chi_D[0] / h0.chi_R[0])} times as spread out as Republicans. Today it is the reverse: Republicans are about ${Math.round(h1.chi_R[0] / h1.chi_D[0])} times as spread out as Democrats.`;
}
// what an index value looks like: two equal sides whose spread gives that circular variance
function drawExplain() {
  const svg = $("#explain-svg"); clear(svg);
  [[10, 100], [25, 300], [40, 500]].forEach(([v, cx]) => {
    const half = Math.acos(1 - v / 100);   // two point blocs at +-half have circular variance v/100
    el("circle", { cx, cy: 86, r: 62, fill: "none", stroke: "var(--line)", "stroke-dasharray": "2 4" }, svg);
    for (const [side, col] of [[-1, "var(--D)"], [1, "var(--R)"]]) for (let i = -3; i <= 3; i++) {
      const a = side * half + i * 0.07;
      el("circle", { cx: cx + 62 * Math.sin(a), cy: 86 - 62 * Math.cos(a), r: 4, fill: col }, svg);
    }
    txt(svg, cx, 178, `Index ${v}`, { "text-anchor": "middle", style: "font-weight:700; fill: var(--ink)" });
  });
  txt(svg, 100, 196, "close to a line", { "text-anchor": "middle" });
  txt(svg, 500, 196, "like the House today", { "text-anchor": "middle" });
}
// one ring per Congress, 1987 inside to today outside; optional highlighted member path
async function drawRings(svg, person) {
  const R = await getJSON("rings");
  clear(svg);
  const r0 = 70, dr = 11, rad = (c) => r0 + (c - 100) * dr, cy = 300;
  const g = el("g", {}, svg);
  for (const h of R) {
    const rr = rad(h.congress);
    el("circle", { cx: 300, cy, r: rr, fill: "none", stroke: "var(--line)", "stroke-width": .6 }, g);
    let d = { D: "", R: "", I: "" };
    h.angle.forEach((a, i) => {
      const s = Math.sin(a), c = Math.cos(a);
      d[h.party[i]] += `M${(300 + (rr - 3.5) * s).toFixed(1)} ${(cy - (rr - 3.5) * c).toFixed(1)}L${(300 + (rr + 3.5) * s).toFixed(1)} ${(cy - (rr + 3.5) * c).toFixed(1)}`;
    });
    for (const k of ["D", "R", "I"]) if (d[k]) el("path", { d: d[k], stroke: `var(--${k})`, "stroke-width": 1.2, opacity: person ? .18 : .4, fill: "none" }, g);
  }
  for (const c of [100, 108, 116, 119]) txt(svg, 300, cy + rad(c) + 4, String(1787 + 2 * c), { "text-anchor": "middle", class: "tick", style: "paint-order: stroke; stroke: var(--surface); stroke-width: 4px" });
  if (person) {
    const pts = person.h.filter((h) => h[11] !== null && h[11] !== undefined).map((h) => [h, 300 + rad(h[0]) * Math.sin(h[11]), cy - rad(h[0]) * Math.cos(h[11])]);
    if (pts.length > 1) el("polyline", { points: pts.map(([, x, y]) => `${x.toFixed(1)},${y.toFixed(1)}`).join(" "), fill: "none", stroke: "var(--ink)", "stroke-width": 1.6 }, svg);
    for (const [h, x, y] of pts) {
      el("circle", { cx: x, cy: y, r: 4.5, fill: `var(--${h[1]})`, stroke: "var(--ink)", "stroke-width": 1.4 }, svg);
      const hit = el("circle", { cx: x, cy: y, r: 10, fill: "transparent" }, svg);
      hit.addEventListener("pointermove", (e) => showTip(e, `<b>${ord(h[0])} Congress</b> (${years(h[0])})<br>Rank ${h[4]} of ${h[7]} (range ${h[5]}–${h[6]})`));
      hit.addEventListener("pointerleave", hideTip);
    }
  }
}
// line vs circle: the same "middle against both ends" vote
function drawWhy() {
  const svg = $("#why-svg"); if (!svg) return; clear(svg);
  const cols = ["var(--D)", "var(--D)", "var(--D)", "var(--D)", "var(--R)", "var(--R)", "var(--R)", "var(--R)"];
  // line: one cut, left vs right
  el("line", { x1: 20, x2: 290, y1: 100, y2: 100, stroke: "var(--ink3)", "stroke-width": 2 }, svg);
  cols.forEach((c, i) => el("circle", { cx: 35 + i * 34, cy: 100, r: 8, fill: c }, svg));
  el("line", { x1: 153, x2: 153, y1: 60, y2: 140, stroke: "var(--ink)", "stroke-width": 2, "stroke-dasharray": "5 4" }, svg);
  txt(svg, 86, 52, "one side", { "text-anchor": "middle" }); txt(svg, 222, 52, "other side", { "text-anchor": "middle" });
  txt(svg, 155, 222, "A line: left against right", { "text-anchor": "middle", style: "font-weight:700; fill: var(--ink)" });
  txt(svg, 155, 240, "both ends can never vote together", { "text-anchor": "middle" });
  // circle: the no arc wraps around the back
  const cx = 480, cy = 100, R = 58, p = (a, r = R) => `${cx + r * Math.sin(a)} ${cy - r * Math.cos(a)}`;
  el("path", { d: `M ${p(-0.75, R + 14)} A ${R + 14} ${R + 14} 0 0 1 ${p(0.75, R + 14)}`, fill: "none", stroke: "var(--yes)", "stroke-width": 5, "stroke-linecap": "round" }, svg);
  el("path", { d: `M ${p(0.95, R + 14)} A ${R + 14} ${R + 14} 0 1 1 ${p(-0.95, R + 14)}`, fill: "none", stroke: "var(--no)", "stroke-width": 5, "stroke-linecap": "round" }, svg);
  el("circle", { cx, cy, r: R, fill: "none", stroke: "var(--ink3)", "stroke-width": 2 }, svg);
  [-1.65, -1.25, -0.6, -0.25, 0.25, 0.6, 1.25, 1.65].forEach((a, i) => el("circle", { cx: cx + R * Math.sin(a), cy: cy - R * Math.cos(a), r: 8, fill: cols[i] }, svg));
  txt(svg, cx, cy - R - 22, "yes: the middle", { "text-anchor": "middle", style: "fill: var(--yes); font-weight:600" });
  txt(svg, cx, cy + R + 34, "no: both ends", { "text-anchor": "middle", style: "fill: var(--no); font-weight:600" });
  txt(svg, cx, 222, "A circle: the middle against both ends", { "text-anchor": "middle", style: "font-weight:700; fill: var(--ink)" });
}

/* ---------- views and routing ---------- */
const VIEWS = ["story", "explore", "members", "member", "methods"];
function show(view) {
  for (const v of VIEWS) $(`#view-${v}`).hidden = v !== view;
  const tab = view === "member" ? "members" : view;
  for (const a of document.querySelectorAll("nav.tabs a")) { if (a.dataset.view === tab) a.setAttribute("aria-current", "page"); else a.removeAttribute("aria-current"); }
}
async function route() {
  hideTip();
  const h = (location.hash || "#story").slice(1);
  let m;
  if ((m = h.match(/^explore(?:-(\d{3}))?$/))) { show("explore"); await renderExplore(m[1] ? +m[1] : state.c); }
  else if ((m = h.match(/^member-([A-Z]\d{6})$/))) { show("member"); renderMember(m[1]); }
  else if (h === "members") { show("members"); }
  else if (h === "methods") { show("methods"); }
  else { show("story"); }
  if (!state.firstRoute) window.scrollTo(0, 0);
  state.firstRoute = false;
}
const state = { c: LATEST, mode: "party", vote: {}, pendingVote: null, selected: null, firstRoute: true };

/* ---------- story ---------- */
async function renderStory() {
  drawTrends(); drawExplain(); drawWhy();
  drawRings($("#rings"));
  const H = await house(LATEST);
  drawCircle($("#hero-circle"), H, { cy: 300, tip: (m) => memberTip(H, m), click: (m) => { location.hash = `member-${m.id}`; } });
  const [H118, V118] = await Promise.all([house(118), votesOf(118)]);
  $("#story-tally").innerHTML = tallyHTML(H118, V118.byRoll[837], "N");
  $("#story-vote-link").addEventListener("click", (e) => { e.preventDefault(); state.pendingVote = { c: 118, roll: 837 }; location.hash = "explore-118"; });
}

/* ---------- explore ---------- */
let exHandle = null;
async function renderExplore(c) {
  if (!houseInfo(c)) c = LATEST;
  state.c = c;
  const info = houseInfo(c);
  if (state.pendingVote && state.pendingVote.c === c) { state.vote[c] = state.pendingVote.roll; state.mode = "vote"; state.pendingVote = null; }
  $("#ex-select").value = c;
  $("#ex-prev").disabled = c <= 100; $("#ex-next").disabled = c >= LATEST;
  $("#ex-title").textContent = `${ord(c)} Congress`;
  $("#ex-sub").textContent = `${years(c)}${c === LATEST ? ", in progress" : ""} · ${PARTY[info.majority]} majority`;
  const nEst = info.nFitted + info.nPlaced;
  $("#ex-tiles").innerHTML = `
    <div class="tile"><div class="k">Polarization index</div><div class="v">${IDX(info.chi_0[0])} <small class="muted" style="font-size:15px;font-weight:500">/ 100</small></div><div class="s">range ${IDX(info.chi_0[1])}–${IDX(info.chi_0[2])}</div></div>
    <div class="tile"><div class="k">Members estimated</div><div class="v">${nEst}</div><div class="s">${info.nD} D · ${info.nR} R${info.nI ? ` · ${info.nI} I` : ""}${info.nSkipped ? ` · ${info.nSkipped} with too few votes` : ""}</div></div>
    <div class="tile"><div class="k">Within-party index</div><div class="v"><span style="color:var(--D)">${(100 * info.chi_D[0]).toFixed(1)}</span> · <span style="color:var(--R)">${(100 * info.chi_R[0]).toFixed(1)}</span></div><div class="s">Democrats · Republicans</div></div>
    <div class="tile"><div class="k">Roll-call votes</div><div class="v">${info.nRolls.toLocaleString()}</div><div class="s">${fmtDate(info.first)} to ${fmtDate(info.last)}</div></div>`;
  const notes = [];
  if (c === LATEST) notes.push(`In session. Votes through ${fmtDate(info.last)}.`);
  if (c === 103) notes.push("Pat Schroeder's position is uncertain (see Methods).");
  $("#ex-note").hidden = !notes.length; $("#ex-note").textContent = notes.join(" ");
  $("#ex-filter").value = ""; $("#ex-vsearch").value = ""; lineupAll = false; $("#ex-vresults").innerHTML = "";
  if (state.selected && state.selected.c !== c) state.selected = null;
  const H = await house(c);
  if (state.c !== c) return;
  renderFeatured(H, info);
  renderLineup(H, "");
  await drawExplore(H);
}
async function drawExplore(H) {
  const c = H.congress, roll = state.vote[c];
  let vs = null;
  if (state.mode === "vote" && roll) { const V = await votesOf(c); if (state.c !== c) return; vs = V.byRoll[roll]; }
  for (const b of document.querySelectorAll(".toggle button")) b.setAttribute("aria-pressed", b.dataset.mode === state.mode);
  exHandle = drawCircle($("#ex-circle"), H, {
    cy: 300, votes: vs, tip: (m) => memberTip(H, m, vs),
    click: (m) => { state.selected = { c, id: m.id }; renderSelected(H, vs); exHandle.mark(m.id); },
    keep: () => !!state.selected, selected: state.selected && state.selected.id,
  });
  $("#ex-legend").innerHTML = vs
    ? `<span><i style="background:var(--yes)"></i>Voted yes</span><span><i style="background:var(--no)"></i>Voted no</span><span><i style="box-shadow:inset 0 0 0 1.5px var(--ink3)"></i>Didn't vote (absent or present)</span>`
    : `<span><i style="background:var(--D)"></i>Democrat</span><span><i style="background:var(--R)"></i>Republican</span>${H.list.some((m) => m.party === "I") ? `<span><i style="background:var(--I)"></i>Independent</span>` : ""}`;
  renderVoteCard(H, vs);
  renderSelected(H, vs);
  for (const b of document.querySelectorAll("#ex-featured button, #ex-vresults button")) b.setAttribute("aria-pressed", state.mode === "vote" && +b.dataset.roll === roll);
}
function rollTitle(r, f) { return (f && f.title) || r.desc || r.question || `Roll call ${r.roll}`; }
function renderVoteCard(H, vs) {
  const box = $("#ex-vote"), roll = state.vote[H.congress];
  if (state.mode !== "vote") { box.hidden = true; return; }
  box.hidden = false;
  if (!roll) { box.innerHTML = `<p class="vd">Pick a vote below.</p>`; return; }
  const r = H.rollList.find((x) => x.roll === roll);
  const f = (houseInfo(H.congress).featured || []).find((x) => x.roll === roll);
  box.innerHTML = `<div class="vt">${esc(rollTitle(r, f))}</div>
    <div class="vm">${fmtDate(r.date)}${r.bill ? ` · ${bill(r.bill)}` : ""} · ${esc(r.question)} · ${esc(r.result || "")} ${r.yea}–${r.nay}</div>
    ${f && f.text ? `<p class="vd">${esc(f.text)}</p>` : ""}
    <div class="small muted">Yes–no split, in %</div><div class="tally">${splitHTML(H, vs)}</div>
    <button class="linkbtn" id="ex-clearvote">Color by party instead</button>`;
  $("#ex-clearvote").onclick = () => { state.mode = "party"; drawExplore(H); };
}
function renderSelected(H, vs) {
  const box = $("#ex-selected");
  const s = state.selected && state.selected.c === H.congress ? H.byId[state.selected.id] : null;
  if (!s) { box.hidden = true; return; }
  const v = vs ? vs[s.k] : null;
  box.hidden = false;
  box.innerHTML = `<b>${esc(s.name)}</b><p>${PARTY[s.party]}, ${dist(s.state, s.d)} · rank ${s.rank} of ${H.nRank} (range ${s.lo}–${s.hi})${vs ? ` · voted ${v === "Y" ? "yes" : v === "N" ? "no" : "did not vote"}` : ""}</p><p><a href="#member-${s.id}">Open member page</a></p>`;
}
function rollButton(r, f) {
  return `<li><button data-roll="${r.roll}" aria-pressed="false"><span class="t">${esc(rollTitle(r, f)).slice(0, 140)}</span><span class="m">${fmtDate(r.date)}${r.bill ? ` · ${bill(r.bill)}` : ""}${f ? "" : ` · ${esc(r.question)}`} · ${r.yea}–${r.nay}</span></button></li>`;
}
function renderFeatured(H, info) {
  const fs = info.featured || [];
  $("#ex-featured").innerHTML = fs.length ? fs.map((f) => rollButton(H.rollList.find((r) => r.roll === f.roll), f)).join("") : `<li class="small muted">No strong examples in this Congress.</li>`;
}
function pickVote(roll) { state.vote[state.c] = roll; state.mode = "vote"; house(state.c).then(drawExplore); }
function searchVotes(H, q) {
  const box = $("#ex-vresults");
  const terms = fold(q).replace(/\./g, "").split(/\s+/).filter(Boolean);
  if (!terms.length) { box.innerHTML = ""; return; }
  const joined = terms.join("");
  const hits = H.rollList.filter((r) => terms.every((t) => r.key.includes(t)) || r.key.includes(joined));
  box.innerHTML = hits.length ? hits.slice(0, 60).map((r) => rollButton(r)).join("") + (hits.length > 60 ? `<li class="none">${hits.length - 60} more. Add a word to narrow the list.</li>` : "")
    : `<li class="none">No vote matches "${esc(q)}". Try a bill number such as H.R. 6172, or a word like budget.</li>`;
}
let lineupAll = false;
function renderLineup(H, q) {
  const f = fold(q.trim());
  const rows = H.list.slice().sort((a, b) => a.rank - b.rank || a.name.localeCompare(b.name))
    .filter((m) => !f || fold(m.name).includes(f) || m.state.toLowerCase() === f || fold(STATES[m.state] || "").includes(f));
  const more = !f && !lineupAll && rows.length > 40;
  $("#ex-more").hidden = !more; $("#ex-more").textContent = `Show all ${rows.length} members`;
  const glyph = (m) => `<svg class="glyph" viewBox="0 0 22 22" aria-hidden="true"><circle cx="11" cy="11" r="8" fill="none" stroke="var(--line)" stroke-width="1.5"/><circle cx="${(11 + 8 * Math.sin(m.angle)).toFixed(1)}" cy="${(11 - 8 * Math.cos(m.angle)).toFixed(1)}" r="3" fill="var(--${m.party})"/></svg>`;
  $("#ex-lineup").innerHTML = rows.length ? (more ? rows.slice(0, 40) : rows).map((m) => `<a class="lrow" role="listitem" href="#member-${m.id}"><span class="r">${m.rank}</span>${glyph(m)}<span class="n">${esc(m.name)} <small>${dist(m.state, m.d)}</small></span><span class="rg">${m.lo}–${m.hi}</span></a>`).join("")
    : `<p class="muted">No member matches "${esc(q)}".</p>`;
}

/* ---------- members ---------- */
function latestRanked(p) { for (let i = p.h.length - 1; i >= 0; i--) if (p.h[i][4] !== null) return p.h[i]; return null; }
function sentence(h) {
  const [, , , , rank, , , n] = h;
  const p = (rank - 1) / (n - 1);
  if (rank === 1) return "The most liberal-voting member of the House";
  if (rank === n) return "The most conservative-voting member of the House";
  if (p > 0.45 && p < 0.55) return "Votes near the middle of the House";
  if (p <= 0.45) return `Votes more liberal than ${Math.round(100 * (n - rank) / (n - 1))}% of other members`;
  return `Votes more conservative than ${Math.round(100 * (rank - 1) / (n - 1))}% of other members`;
}
function personMeta(p) {
  const last = p.h[p.h.length - 1], first = p.h[0];
  return `${last[1]} · ${dist(last[2], last[3])} · ${1787 + 2 * first[0]}–${last[0] === LATEST ? "now" : 1789 + 2 * last[0]}`;
}
function bindMemberSearch() {
  const q = $("#m-q"), res = $("#m-res");
  let hits = [], on = -1;
  const draw = () => {
    res.innerHTML = hits.length ? hits.map((p, i) => `<li><button data-id="${p.id}" class="${i === on ? "on" : ""}"><span>${esc(p.name)}</span><span class="m">${personMeta(p)}</span></button></li>`).join("")
      : `<li class="none">No member matches "${esc(q.value)}". Try a last name.</li>`;
    res.hidden = false;
  };
  q.addEventListener("input", () => {
    const s = fold(q.value.trim());
    if (s.length < 2) { res.hidden = true; return; }
    const cur = (p) => p.h[p.h.length - 1][0] === LATEST ? 0 : 1;
    hits = PEOPLE.filter((p) => p.key.includes(s)).sort((a, b) => cur(a) - cur(b) || a.name.localeCompare(b.name)).slice(0, 15);
    on = -1; draw();
  });
  q.addEventListener("keydown", (e) => {
    if (res.hidden || !hits.length) return;
    if (e.key === "ArrowDown") { on = Math.min(on + 1, hits.length - 1); draw(); e.preventDefault(); }
    else if (e.key === "ArrowUp") { on = Math.max(on - 1, 0); draw(); e.preventDefault(); }
    else if (e.key === "Enter") { const p = hits[Math.max(on, 0)]; location.hash = `member-${p.id}`; res.hidden = true; q.value = ""; }
    else if (e.key === "Escape") res.hidden = true;
  });
  res.addEventListener("click", (e) => { const b = e.target.closest("button"); if (!b) return; res.hidden = true; q.value = ""; location.hash = `member-${b.dataset.id}`; });
  document.addEventListener("click", (e) => { if (!e.target.closest(".search")) res.hidden = true; });
  const sel = $("#m-state");
  sel.innerHTML += Object.entries(STATES).map(([k, v]) => `<option value="${k}">${v}</option>`).join("");
  sel.addEventListener("change", () => {
    const st = sel.value;
    if (!st) { $("#m-statelist").innerHTML = ""; return; }
    const rows = [];
    for (const p of PEOPLE) for (const h of p.h) if (h[0] === LATEST && h[2] === st) rows.push([p, h]);
    rows.sort((a, b) => a[1][3] - b[1][3] || a[0].name.localeCompare(b[0].name));
    $("#m-statelist").innerHTML = rows.map(([p, h]) => `<a class="scard" href="#member-${p.id}"><b>${esc(p.name)}</b><span>${PARTY[h[1]]} · ${h[3] === 0 ? "At-large" : "District " + h[3]}</span><span>${h[4] === null ? "Too few votes to estimate" : sentence(h)}</span></a>`).join("")
      + `<p class="small muted" style="grid-column:1/-1">119th Congress districts. Your 2026 district may differ.</p>`;
  });
}
function renderMember(id) {
  const p = BYID[id], box = $("#member");
  if (!p) { box.innerHTML = `<p>No member with this id. <a href="#members">Search for a member</a>.</p>`; return; }
  document.title = `${p.name} · Full Circle`;
  const last = p.h[p.h.length - 1], lr = latestRanked(p);
  const parties = [...new Set(p.h.map((h) => h[1]))];
  const pct = (r, n) => 100 * (r - 1) / (n - 1);
  const notes = [];
  if (last[0] === LATEST) notes.push(`Votes through ${fmtDate(houseInfo(LATEST).last)}. Congress still in session.`);
  if (lr && (lr[9] || lr[6] - lr[5] > lr[7] * 0.2)) notes.push(`Wide range${lr[9] ? " (missed many votes)" : ""}: treat as approximate.`);
  for (const h of p.h) if (h[10]) notes.push(`${ord(h[0])} Congress: ${h[10]}.`);
  if (lr && lr !== last) notes.push(`Too few votes in the ${ord(last[0])} Congress; showing the ${ord(lr[0])}.`);
  const card = lr ? `
      <p class="big">${sentence(lr)}${lr === last && last[0] === LATEST ? "" : ` (${ord(lr[0])} Congress)`}.</p>
      <div class="scale" aria-hidden="true"><div class="ln"></div><div class="ci" style="left:${pct(lr[5], lr[7])}%;width:${Math.max(pct(lr[6], lr[7]) - pct(lr[5], lr[7]), 1.5)}%"></div><div class="pt" style="left:${pct(lr[4], lr[7])}%"></div></div>
      <div class="ends"><span>Most liberal</span><span>Most conservative</span></div>
      <dl class="facts">
        <dt>Rank</dt><dd><span class="mono">${lr[4]}</span> of <span class="mono">${lr[7]}</span> (1 = most liberal)</dd>
        <dt>Likely range</dt><dd><span class="mono">${lr[5]}</span>–<span class="mono">${lr[6]}</span></dd>
        <dt>Votes counted</dt><dd class="mono">${lr[8].toLocaleString()}</dd>
      </dl>` : `<p class="big">Too few votes to estimate a position.</p>`;
  box.innerHTML = `
    <div class="mcard">
      <h2 class="who">${esc(p.name)}</h2>
      <div class="sub"><span class="pill" style="background:var(--${last[1]})">${PARTY[last[1]]}</span><span>${STATES[last[2]] || last[2]}, ${last[3] === 0 ? "at-large" : "district " + last[3]}</span>${parties.length > 1 ? `<span class="small muted">Also served as ${parties.filter((x) => x !== last[1]).map((x) => PARTY[x]).join(", ")}</span>` : ""}</div>
      ${card}
      ${notes.map((n) => `<p class="note flag">${esc(n)}</p>`).join("")}
      <p class="note">Based on floor votes only. <a href="https://bioguide.congress.gov/search/bio/${p.id}" rel="noopener" target="_blank">Official biography</a></p>
    </div>
    <div>
      <figure class="panel">
        <figcaption><b>Seat on the circle, ${ord(lr ? lr[0] : last[0])} Congress</b><span>Ringed dot: ${esc(p.name)}</span></figcaption>
        <svg id="m-circle" class="circ" viewBox="0 0 600 560" role="img" aria-label="Member highlighted on the circle"></svg>
      </figure>
      <figure class="panel" style="margin-top:16px">
        <figcaption><b>Over the years</b><span>One ring per Congress, 1987 inside. Dots: this member.</span></figcaption>
        <svg id="m-rings" class="circ" viewBox="0 0 600 600" role="img" aria-label="Member's position on each Congress ring"></svg>
      </figure>
      <div class="tablewrap"><table class="htable">
        <thead><tr><th>Congress</th><th>Party</th><th>District</th><th>Rank</th><th>Range</th><th>Votes</th></tr></thead>
        <tbody>${p.h.slice().reverse().map((h) => `<tr><td><a href="#explore-${h[0]}">${ord(h[0])}</a> <span class="muted small">${1787 + 2 * h[0]}</span></td><td>${h[1]}</td><td>${dist(h[2], h[3])}</td><td>${h[4] === null ? "too few votes" : `${h[4]} of ${h[7]}`}</td><td>${h[4] === null ? "" : `${h[5]}–${h[6]}`}</td><td>${h[8].toLocaleString()}</td></tr>`).join("")}</tbody>
      </table></div>
    </div>`;
  drawRings($("#m-rings"), p);
  const hc = lr ? lr[0] : last[0];
  house(hc).then((H) => { if (!$("#m-circle")) return; drawCircle($("#m-circle"), H, { cy: 300, tip: (m) => memberTip(H, m), click: (m) => { location.hash = `member-${m.id}`; }, selected: id, keep: () => true }); });
}
/* ---------- start ---------- */
async function start() {
  try {
    [INDEX, PEOPLE] = await Promise.all([getJSON("index"), getJSON("members")]);
  } catch (e) {
    $("#main").insertAdjacentHTML("afterbegin", `<p class="callout wrap">The data did not load (${esc(e.message)}). Reload the page to try again.</p>`);
    return;
  }
  BYID = Object.fromEntries(PEOPLE.map((p) => [p.id, p]));
  for (const p of PEOPLE) p.key = fold(p.name);
  for (const s of document.querySelectorAll(".lastvote")) s.textContent = fmtDate(INDEX.lastVote);
  const sel = $("#ex-select");
  sel.innerHTML = INDEX.houses.slice().reverse().map((h) => `<option value="${h.congress}">${ord(h.congress)} Congress (${years(h.congress)})</option>`).join("");
  sel.addEventListener("change", () => { location.hash = `explore-${sel.value}`; });
  $("#ex-prev").addEventListener("click", () => { location.hash = `explore-${state.c - 1}`; });
  $("#ex-next").addEventListener("click", () => { location.hash = `explore-${state.c + 1}`; });
  for (const b of document.querySelectorAll(".toggle button")) b.addEventListener("click", () => { state.mode = b.dataset.mode; house(state.c).then(drawExplore); });
  const pick = (e) => { const b = e.target.closest("button[data-roll]"); if (b) pickVote(+b.dataset.roll); };
  $("#ex-featured").addEventListener("click", pick);
  $("#ex-vresults").addEventListener("click", pick);
  $("#ex-vsearch").addEventListener("input", (e) => house(state.c).then((H) => searchVotes(H, e.target.value)));
  $("#ex-filter").addEventListener("input", (e) => house(state.c).then((H) => renderLineup(H, e.target.value)));
  $("#ex-more").addEventListener("click", () => { lineupAll = true; house(state.c).then((H) => renderLineup(H, "")); });
  bindMemberSearch();
  window.addEventListener("hashchange", () => { document.title = "Full Circle"; route(); });
  route();
  renderStory();
}
start();
