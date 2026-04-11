import { useState, useEffect, useRef, useCallback } from "react";
import * as d3 from "d3";

// ─── BLOSUM62 Matrix ─────────────────────────────────────────────────────────
const BLOSUM62 = {
  A:{A:4,R:-1,N:-2,D:-2,C:0,Q:-1,E:-1,G:0,H:-2,I:-1,L:-1,K:-1,M:-1,F:-2,P:-1,S:1,T:0,W:-3,Y:-2,V:0},
  R:{A:-1,R:5,N:0,D:-2,C:-3,Q:1,E:0,G:-2,H:0,I:-3,L:-2,K:2,M:-1,F:-3,P:-2,S:-1,T:-1,W:-3,Y:-2,V:-3},
  N:{A:-2,R:0,N:6,D:1,C:-3,Q:0,E:0,G:0,H:1,I:-3,L:-3,K:0,M:-2,F:-3,P:-2,S:1,T:0,W:-4,Y:-2,V:-3},
  D:{A:-2,R:-2,N:1,D:6,C:-3,Q:0,E:2,G:-1,H:-1,I:-3,L:-4,K:-1,M:-3,F:-3,P:-1,S:0,T:-1,W:-4,Y:-3,V:-3},
  C:{A:0,R:-3,N:-3,D:-3,C:9,Q:-3,E:-4,G:-3,H:-3,I:-1,L:-1,K:-3,M:-1,F:-2,P:-3,S:-1,T:-1,W:-2,Y:-2,V:-1},
  Q:{A:-1,R:1,N:0,D:0,C:-3,Q:5,E:2,G:-2,H:0,I:-3,L:-2,K:1,M:0,F:-3,P:-1,S:0,T:-1,W:-2,Y:-1,V:-2},
  E:{A:-1,R:0,N:0,D:2,C:-4,Q:2,E:5,G:-2,H:0,I:-3,L:-3,K:1,M:-2,F:-3,P:-1,S:0,T:-1,W:-3,Y:-2,V:-2},
  G:{A:0,R:-2,N:0,D:-1,C:-3,Q:-2,E:-2,G:6,H:-2,I:-4,L:-4,K:-2,M:-3,F:-3,P:-2,S:0,T:-2,W:-2,Y:-3,V:-3},
  H:{A:-2,R:0,N:1,D:-1,C:-3,Q:0,E:0,G:-2,H:8,I:-3,L:-3,K:-1,M:-2,F:-1,P:-2,S:-1,T:-2,W:-2,Y:2,V:-3},
  I:{A:-1,R:-3,N:-3,D:-3,C:-1,Q:-3,E:-3,G:-4,H:-3,I:4,L:2,K:-3,M:1,F:0,P:-3,S:-2,T:-1,W:-3,Y:-1,V:3},
  L:{A:-1,R:-2,N:-3,D:-4,C:-1,Q:-2,E:-3,G:-4,H:-3,I:2,L:4,K:-2,M:2,F:0,P:-3,S:-2,T:-1,W:-2,Y:-1,V:1},
  K:{A:-1,R:2,N:0,D:-1,C:-3,Q:1,E:1,G:-2,H:-1,I:-3,L:-2,K:5,M:-1,F:-3,P:-1,S:0,T:-1,W:-3,Y:-2,V:-2},
  M:{A:-1,R:-1,N:-2,D:-3,C:-1,Q:0,E:-2,G:-3,H:-2,I:1,L:2,K:-1,M:5,F:0,P:-2,S:-1,T:-1,W:-1,Y:-1,V:1},
  F:{A:-2,R:-3,N:-3,D:-3,C:-2,Q:-3,E:-3,G:-3,H:-1,I:0,L:0,K:-3,M:0,F:6,P:-4,S:-2,T:-2,W:1,Y:3,V:-1},
  P:{A:-1,R:-2,N:-2,D:-1,C:-3,Q:-1,E:-1,G:-2,H:-2,I:-3,L:-3,K:-1,M:-2,F:-4,P:7,S:-1,T:-1,W:-4,Y:-3,V:-2},
  S:{A:1,R:-1,N:1,D:0,C:-1,Q:0,E:0,G:0,H:-1,I:-2,L:-2,K:0,M:-1,F:-2,P:-1,S:4,T:1,W:-3,Y:-2,V:-2},
  T:{A:0,R:-1,N:0,D:-1,C:-1,Q:-1,E:-1,G:-2,H:-2,I:-1,L:-1,K:-1,M:-1,F:-2,P:-1,S:1,T:5,W:-2,Y:-2,V:0},
  W:{A:-3,R:-3,N:-4,D:-4,C:-2,Q:-2,E:-3,G:-2,H:-2,I:-3,L:-2,K:-3,M:-1,F:1,P:-4,S:-3,T:-2,W:11,Y:2,V:-3},
  Y:{A:-2,R:-2,N:-2,D:-3,C:-2,Q:-1,E:-2,G:-3,H:2,I:-1,L:-1,K:-2,M:-1,F:3,P:-3,S:-2,T:-2,W:2,Y:7,V:-1},
  V:{A:0,R:-3,N:-3,D:-3,C:-1,Q:-2,E:-2,G:-3,H:-3,I:3,L:1,K:-2,M:1,F:-1,P:-2,S:-2,T:0,W:-3,Y:-1,V:4},
};

// ─── Bioinformatics Engine ───────────────────────────────────────────────────
function parseFASTA(text) {
  const sequences = [];
  let currentHeader = null, currentSeq = [];
  for (const line of text.split(/\r?\n/)) {
    const t = line.trim();
    if (!t) continue;
    if (t.startsWith(">")) {
      if (currentHeader !== null) sequences.push({ header: currentHeader, sequence: currentSeq.join("") });
      currentHeader = t.slice(1).trim();
      currentSeq = [];
    } else {
      currentSeq.push(t.toUpperCase());
    }
  }
  if (currentHeader !== null) sequences.push({ header: currentHeader, sequence: currentSeq.join("") });
  return sequences;
}

function needlemanWunsch(seqA, seqB, gapPenalty = -2) {
  const n = seqA.length, m = seqB.length;
  const scoreAA = (a, b) => {
    if (a === "-" || b === "-") return gapPenalty;
    if (BLOSUM62[a] && BLOSUM62[a][b] !== undefined) return BLOSUM62[a][b];
    return a === b ? 2 : -1;
  };
  const dp = Array.from({ length: n + 1 }, (_, i) => Array.from({ length: m + 1 }, (_, j) => 0));
  for (let i = 0; i <= n; i++) dp[i][0] = i * gapPenalty;
  for (let j = 0; j <= m; j++) dp[0][j] = j * gapPenalty;
  for (let i = 1; i <= n; i++)
    for (let j = 1; j <= m; j++)
      dp[i][j] = Math.max(
        dp[i-1][j-1] + scoreAA(seqA[i-1], seqB[j-1]),
        dp[i-1][j] + gapPenalty,
        dp[i][j-1] + gapPenalty
      );
  let i = n, j = m;
  const aln1 = [], aln2 = [];
  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && dp[i][j] === dp[i-1][j-1] + scoreAA(seqA[i-1], seqB[j-1])) {
      aln1.unshift(seqA[i-1]); aln2.unshift(seqB[j-1]); i--; j--;
    } else if (i > 0 && dp[i][j] === dp[i-1][j] + gapPenalty) {
      aln1.unshift(seqA[i-1]); aln2.unshift("-"); i--;
    } else {
      aln1.unshift("-"); aln2.unshift(seqB[j-1]); j--;
    }
  }
  return [aln1.join(""), aln2.join("")];
}

function jaccardKmer(a, b, k = 4) {
  const kmers = s => { const set = new Set(); for (let i = 0; i <= s.length - k; i++) set.add(s.slice(i, i + k)); return set; };
  const sa = kmers(a), sb = kmers(b);
  let inter = 0;
  for (const x of sa) if (sb.has(x)) inter++;
  const union = sa.size + sb.size - inter;
  return union === 0 ? 0 : 1 - inter / union;
}

function starMSA(sequences) {
  const n = sequences.length;
  if (n === 1) return [sequences[0].sequence];
  const dists = Array.from({ length: n }, () => Array(n).fill(0));
  for (let i = 0; i < n; i++)
    for (let j = i + 1; j < n; j++) {
      const d = jaccardKmer(sequences[i].sequence, sequences[j].sequence);
      dists[i][j] = dists[j][i] = d;
    }
  const sums = dists.map(row => row.reduce((a, b) => a + b, 0));
  const center = sums.indexOf(Math.min(...sums));
  const pairAlignments = {};
  let centerRef = sequences[center].sequence;
  for (let i = 0; i < n; i++) {
    if (i === center) continue;
    const [ac, ai] = needlemanWunsch(sequences[center].sequence, sequences[i].sequence);
    pairAlignments[i] = { centerAligned: ac, otherAligned: ai };
    if (ac.length > centerRef.length) centerRef = ac;
  }
  const rows = Array.from({ length: n }, () => "");
  const anyPair = Object.values(pairAlignments)[0];
  rows[center] = anyPair ? anyPair.centerAligned : sequences[center].sequence;
  for (let i = 0; i < n; i++) {
    if (i === center) continue;
    rows[i] = pairAlignments[i].otherAligned;
  }
  const maxLen = Math.max(...rows.map(r => r.length));
  return rows.map(r => r.padEnd(maxLen, "-"));
}

function computeDistances(msa, model) {
  const n = msa.length;
  const mat = Array.from({ length: n }, () => Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      let matches = 0, total = 0;
      for (let k = 0; k < msa[i].length; k++) {
        if (msa[i][k] !== "-" && msa[j][k] !== "-") {
          total++;
          if (msa[i][k] === msa[j][k]) matches++;
        }
      }
      const p = total === 0 ? 1 : 1 - matches / total;
      let d;
      if (model === "jaccard") {
        d = jaccardKmer(msa[i].replace(/-/g, ""), msa[j].replace(/-/g, ""));
      } else if (model === "poisson") {
        d = p >= 1 ? 2 : -Math.log(1 - p);
      } else {
        d = p;
      }
      mat[i][j] = mat[j][i] = Math.max(0, d);
    }
  }
  return mat;
}

function neighborJoining(distMatrix, labels) {
  let n = distMatrix.length;
  let D = distMatrix.map(row => [...row]);
  let nodes = labels.map((l) => ({ name: l, children: [], branchLength: 0 }));
  while (nodes.length > 2) {
    const m = nodes.length;
    const r = Array(m).fill(0);
    for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) r[i] += D[i][j];
    let minQ = Infinity, pi = 0, pj = 1;
    for (let i = 0; i < m; i++)
      for (let j = i + 1; j < m; j++) {
        const q = (m - 2) * D[i][j] - r[i] - r[j];
        if (q < minQ) { minQ = q; pi = i; pj = j; }
      }
    const diu = (D[pi][pj] + (r[pi] - r[pj]) / (m - 2)) / 2;
    const dju = D[pi][pj] - diu;
    const newNode = { name: null, children: [nodes[pi], nodes[pj]], branchLength: 0 };
    nodes[pi].branchLength = Math.max(0, diu);
    nodes[pj].branchLength = Math.max(0, dju);
    const remaining = nodes.map((_, idx) => idx).filter(idx => idx !== pi && idx !== pj);
    const newD = [];
    for (const k of remaining) {
      newD.push(Math.max(0, (D[pi][k] + D[pj][k] - D[pi][pj]) / 2));
    }
    const newNodes = remaining.map(idx => nodes[idx]);
    newNodes.push(newNode);
    const sz = newNodes.length;
    const newDM = Array.from({ length: sz }, () => Array(sz).fill(0));
    for (let a = 0; a < remaining.length; a++)
      for (let b = 0; b < remaining.length; b++)
        newDM[a][b] = D[remaining[a]][remaining[b]];
    for (let a = 0; a < remaining.length; a++) {
      newDM[a][sz - 1] = newD[a];
      newDM[sz - 1][a] = newD[a];
    }
    D = newDM;
    nodes = newNodes;
  }
  if (nodes.length === 2) {
    nodes[0].branchLength = D[0][1] / 2;
    nodes[1].branchLength = D[0][1] / 2;
    return { name: null, children: nodes, branchLength: 0 };
  }
  return nodes[0];
}

function upgma(distMatrix, labels) {
  let D = distMatrix.map(row => [...row]);
  let clusters = labels.map((l) => ({ node: { name: l, children: [], branchLength: 0 }, size: 1, height: 0 }));
  while (clusters.length > 1) {
    let minD = Infinity, pi = 0, pj = 1;
    for (let i = 0; i < clusters.length; i++)
      for (let j = i + 1; j < clusters.length; j++)
        if (D[i][j] < minD) { minD = D[i][j]; pi = i; pj = j; }
    const newHeight = minD / 2;
    clusters[pi].node.branchLength = Math.max(0, newHeight - clusters[pi].height);
    clusters[pj].node.branchLength = Math.max(0, newHeight - clusters[pj].height);
    const newCluster = {
      node: { name: null, children: [clusters[pi].node, clusters[pj].node], branchLength: 0 },
      size: clusters[pi].size + clusters[pj].size,
      height: newHeight,
    };
    const remaining = clusters.map((_, idx) => idx).filter(idx => idx !== pi && idx !== pj);
    const sz = remaining.length + 1;
    const newD = Array.from({ length: sz }, () => Array(sz).fill(0));
    for (let a = 0; a < remaining.length; a++)
      for (let b = 0; b < remaining.length; b++)
        newD[a][b] = D[remaining[a]][remaining[b]];
    for (let a = 0; a < remaining.length; a++) {
      const si = clusters[pi].size, sj = clusters[pj].size;
      const d = (si * D[remaining[a]][pi] + sj * D[remaining[a]][pj]) / (si + sj);
      newD[a][sz - 1] = d;
      newD[sz - 1][a] = d;
    }
    D = newD;
    clusters = [...remaining.map(idx => clusters[idx]), newCluster];
  }
  return clusters[0].node;
}

function toNewick(node) {
  if (!node.children || node.children.length === 0) return `${node.name}:${(node.branchLength || 0).toFixed(5)}`;
  return `(${node.children.map(toNewick).join(",")})${node.name || ""}:${(node.branchLength || 0).toFixed(5)}`;
}

function bootstrapSupport(msa, labels, distModel, treeMethod, reps = 100) {
  const len = msa[0].length;
  const getLeafSets = (node) => {
    const collect = (n) => {
      if (!n.children || !n.children.length) return [n.name];
      return n.children.flatMap(collect);
    };
    const result = [];
    const traverse = (n) => {
      if (n.children && n.children.length) {
        result.push(collect(n).sort().join("|"));
        n.children.forEach(traverse);
      }
    };
    traverse(node);
    return result;
  };
  const refClades = new Set(getLeafSets(
    treeMethod === "upgma"
      ? upgma(computeDistances(msa, distModel), labels)
      : neighborJoining(computeDistances(msa, distModel), labels)
  ));
  const counts = {};
  refClades.forEach(k => { counts[k] = 0; });
  for (let r = 0; r < reps; r++) {
    const cols = Array.from({ length: len }, () => Math.floor(Math.random() * len));
    const rMSA = msa.map(seq => cols.map(c => seq[c]).join(""));
    const rDist = computeDistances(rMSA, distModel);
    const rTree = treeMethod === "upgma" ? upgma(rDist, labels) : neighborJoining(rDist, labels);
    const rClades = new Set(getLeafSets(rTree));
    for (const k of refClades) if (rClades.has(k)) counts[k]++;
  }
  const result = {};
  for (const k of refClades) result[k] = Math.round((counts[k] / reps) * 100);
  return result;
}

// ─── Demo FASTA ──────────────────────────────────────────────────────────────
const DEMO_FASTA = `>Human
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
>Chimpanzee
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
>Gorilla
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGQG
>Mouse
MVHLTDAEKAAVSCLWGKVNSDEVGGEALGRLLVVYPWTQRFFDSFGNLSSASAIMGNPK
>Rat
MVHLTDAEKAAVSCLWGKVNSDEVGGEALGRLLVVYPWTQRFFDSFGNLSSASAIMGNPQ`;

// ─── AI Analysis ─────────────────────────────────────────────────────────────
async function analyzeWithAI(newick, labels) {
  const prompt = `You are an evolutionary biologist. Analyze this phylogenetic tree.
Newick: ${newick}
Taxa: ${labels.join(", ")}

Respond ONLY with valid JSON (no markdown, no backticks, no explanation outside JSON):
{"summary":"2-3 sentence evolutionary summary","clades":[{"name":"clade name","members":["sp1","sp2"],"interpretation":"biological meaning"}],"keyFindings":["finding1","finding2","finding3"]}`;

  const response = await fetch("https://api.anthropic.com/v1/messages", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      model: "claude-sonnet-4-20250514",
      max_tokens: 1000,
      messages: [{ role: "user", content: prompt }],
    }),
  });
  const data = await response.json();
  const text = (data.content || []).map(b => b.text || "").join("").replace(/```json|```/g, "").trim();
  try { return JSON.parse(text); } catch { return { summary: text, clades: [], keyFindings: [] }; }
}

// ─── Tree Visualizer ─────────────────────────────────────────────────────────
function TreeVisualizer({ treeData, bootstrapValues, circular, phylogram, onNodeClick, annotations }) {
  const svgRef = useRef(null);

  useEffect(() => {
    if (!treeData || !svgRef.current) return;
    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();

    const margin = { top: 30, right: 180, bottom: 30, left: 50 };
    const leaves = d3.hierarchy(treeData, d => d.children).leaves().length;
    const height = Math.max(360, leaves * 44);
    const width = 780;
    const iW = width - margin.left - margin.right;
    const iH = height - margin.top - margin.bottom;

    svg.attr("viewBox", `0 0 ${width} ${height}`).attr("width", "100%").attr("height", height);

    const defs = svg.append("defs");
    const glowFilter = defs.append("filter").attr("id", "glow");
    glowFilter.append("feGaussianBlur").attr("stdDeviation", "3").attr("result", "coloredBlur");
    const merge = glowFilter.append("feMerge");
    merge.append("feMergeNode").attr("in", "coloredBlur");
    merge.append("feMergeNode").attr("in", "SourceGraphic");

    const purpleGlow = defs.append("filter").attr("id", "purpleGlow");
    purpleGlow.append("feGaussianBlur").attr("stdDeviation", "4").attr("result", "coloredBlur");
    const merge2 = purpleGlow.append("feMerge");
    merge2.append("feMergeNode").attr("in", "coloredBlur");
    merge2.append("feMergeNode").attr("in", "SourceGraphic");

    const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

    if (circular) {
      const radius = Math.min(iW, iH) / 2 - 20;
      const cx = iW / 2, cy = iH / 2;
      const tree = d3.cluster().size([2 * Math.PI, radius]);
      const root = d3.hierarchy(treeData, d => d.children);
      tree(root);

      g.append("g").selectAll("path")
        .data(root.links())
        .join("path")
        .attr("fill", "none")
        .attr("stroke", "#1e3a8a")
        .attr("stroke-width", 1.5)
        .attr("d", d3.linkRadial().angle(d => d.x).radius(d => d.y))
        .attr("transform", `translate(${cx},${cy})`);

      const node = g.append("g").selectAll("g")
        .data(root.descendants())
        .join("g")
        .attr("transform", d => {
          const x = cx + d.y * Math.cos(d.x - Math.PI / 2);
          const y = cy + d.y * Math.sin(d.x - Math.PI / 2);
          return `translate(${x},${y})`;
        })
        .style("cursor", "pointer")
        .on("click", (_, d) => onNodeClick(d.data));

      node.append("circle")
        .attr("r", d => d.children ? 3 : 5)
        .attr("fill", d => annotations[d.data.name] ? "#a855f7" : (d.children ? "#1d4ed8" : "#06b6d4"))
        .attr("stroke", d => d.children ? "#3b82f6" : "#22d3ee")
        .attr("stroke-width", 1)
        .attr("filter", d => annotations[d.data.name] ? "url(#purpleGlow)" : "url(#glow)");

      node.filter(d => !d.children).append("text")
        .attr("dy", "0.31em")
        .attr("x", d => d.x < Math.PI ? 10 : -10)
        .attr("text-anchor", d => d.x < Math.PI ? "start" : "end")
        .attr("transform", d => d.x >= Math.PI ? "rotate(180)" : null)
        .text(d => d.data.name || "")
        .attr("fill", d => annotations[d.data.name] ? "#c084fc" : "#94a3b8")
        .attr("font-size", "11px")
        .attr("font-family", "monospace");
    } else {
      const tree = d3.cluster().size([iH, iW]);
      const root = d3.hierarchy(treeData, d => d.children);
      tree(root);

      if (phylogram) {
        const getDepth = (node) => {
          let d = 0, n = node;
          while (n.parent) { d += n.data.branchLength || 0; n = n.parent; }
          return d;
        };
        const maxD = d3.max(root.descendants(), getDepth) || 1;
        root.each(d => { d.y = (getDepth(d) / maxD) * iW; });
      }

      // Links (elbow style)
      g.append("g").selectAll("path")
        .data(root.links())
        .join("path")
        .attr("fill", "none")
        .attr("stroke", "#1e3a8a")
        .attr("stroke-width", 1.8)
        .attr("stroke-opacity", 0.9)
        .attr("d", d => `M${d.target.y},${d.target.x}H${d.source.y}V${d.source.x}`);

      // Bootstrap labels
      if (bootstrapValues && Object.keys(bootstrapValues).length > 0) {
        g.append("g").selectAll("text.bs")
          .data(root.descendants().filter(d => d.children && d.depth > 0))
          .join("text")
          .attr("class", "bs")
          .attr("x", d => d.y - 4)
          .attr("y", d => d.x - 5)
          .attr("font-size", "9px")
          .attr("text-anchor", "end")
          .attr("fill", d => {
            const key = d.leaves().map(l => l.data.name).sort().join("|");
            const v = bootstrapValues[key] || 0;
            return v >= 90 ? "#34d399" : v >= 70 ? "#fbbf24" : "#f87171";
          })
          .text(d => {
            const key = d.leaves().map(l => l.data.name).sort().join("|");
            return bootstrapValues[key] ? `${bootstrapValues[key]}` : "";
          });
      }

      // Nodes
      const node = g.append("g").selectAll("g")
        .data(root.descendants())
        .join("g")
        .attr("transform", d => `translate(${d.y},${d.x})`)
        .style("cursor", "pointer")
        .on("click", (_, d) => onNodeClick(d.data));

      node.append("circle")
        .attr("r", d => d.children ? 4 : 5)
        .attr("fill", d => annotations[d.data.name] ? "#a855f7" : (d.children ? "#1d4ed8" : "#06b6d4"))
        .attr("stroke", d => annotations[d.data.name] ? "#c084fc" : (d.children ? "#3b82f6" : "#22d3ee"))
        .attr("stroke-width", 1.5)
        .attr("filter", d => annotations[d.data.name] ? "url(#purpleGlow)" : "url(#glow)");

      node.filter(d => !d.children).append("text")
        .attr("x", 10)
        .attr("dy", "0.32em")
        .attr("fill", d => annotations[d.data.name] ? "#c084fc" : "#94a3b8")
        .attr("font-size", "12px")
        .attr("font-family", "monospace")
        .text(d => d.data.name || "");

      if (phylogram) {
        node.filter(d => d.parent && (d.data.branchLength || 0) > 0.001).append("text")
          .attr("x", d => -(d.y - (d.parent ? d.parent.y : 0)) / 2)
          .attr("dy", "-5px")
          .attr("text-anchor", "middle")
          .attr("fill", "#334155")
          .attr("font-size", "9px")
          .text(d => d.data.branchLength ? d.data.branchLength.toFixed(3) : "");
      }
    }
  }, [treeData, bootstrapValues, circular, phylogram, annotations]);

  const downloadSVG = () => {
    const el = svgRef.current;
    const clone = el.cloneNode(true);
    clone.style.background = "#0f172a";
    const str = new XMLSerializer().serializeToString(clone);
    const blob = new Blob([str], { type: "image/svg+xml" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = "phylotree.svg"; a.click();
    URL.revokeObjectURL(url);
  };

  return (
    <div style={{ position: "relative" }}>
      <svg ref={svgRef} style={{ background: "transparent", display: "block" }} />
      <button onClick={downloadSVG} style={{ position: "absolute", top: "8px", right: "8px", padding: "4px 12px", fontSize: "11px", fontFamily: "monospace", background: "rgba(6,182,212,0.15)", border: "1px solid rgba(6,182,212,0.4)", borderRadius: "6px", color: "#67e8f9", cursor: "pointer" }}>
        ↓ SVG
      </button>
    </div>
  );
}

// ─── Distance Matrix ─────────────────────────────────────────────────────────
function DistanceMatrix({ matrix, labels }) {
  if (!matrix || !labels) return null;
  const max = Math.max(...matrix.flat().filter(v => v > 0), 0.001);
  return (
    <div style={{ overflowX: "auto" }}>
      <table style={{ fontSize: "10px", borderCollapse: "separate", borderSpacing: "2px" }}>
        <thead>
          <tr>
            <th style={{ padding: "2px 6px", color: "#334155" }}></th>
            {labels.map(l => (
              <th key={l} style={{ padding: "2px 6px", color: "#64748b", fontFamily: "monospace", writingMode: "vertical-lr", maxHeight: "70px", overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}>
                {l.slice(0, 8)}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {matrix.map((row, i) => (
            <tr key={i}>
              <td style={{ padding: "3px 8px", color: "#94a3b8", fontFamily: "monospace", whiteSpace: "nowrap" }}>{labels[i].slice(0, 10)}</td>
              {row.map((val, j) => {
                const intensity = val / max;
                const bg = i === j ? "rgba(6,182,212,0.15)" : `rgba(30,64,175,${0.08 + intensity * 0.65})`;
                return (
                  <td key={j} style={{ padding: "3px 8px", textAlign: "center", background: bg, borderRadius: "3px", color: i === j ? "#06b6d4" : "#94a3b8", minWidth: "48px" }}>
                    {val.toFixed(3)}
                  </td>
                );
              })}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

// ─── Main App ─────────────────────────────────────────────────────────────────
export default function PhyloGenAI() {
  const [sequences, setSequences] = useState(null);
  const [loading, setLoading] = useState(false);
  const [treeRoot, setTreeRoot] = useState(null);
  const [distMatrix, setDistMatrix] = useState(null);
  const [newick, setNewick] = useState("");
  const [bootstrapVals, setBootstrapVals] = useState({});
  const [aiResult, setAiResult] = useState(null);
  const [aiLoading, setAiLoading] = useState(false);
  const [error, setError] = useState("");
  const [activeTab, setActiveTab] = useState("ai");
  const [annotations, setAnnotations] = useState({});
  const [selectedNode, setSelectedNode] = useState(null);
  const [annotLabel, setAnnotLabel] = useState("");
  const [annotNote, setAnnotNote] = useState("");
  const [circular, setCircular] = useState(false);
  const [phylogram, setPhylogram] = useState(true);
  const [distModel, setDistModel] = useState("percent");
  const [treeMethod, setTreeMethod] = useState("nj");
  const [runBootstrap, setRunBootstrap] = useState(false);
  const [dragOver, setDragOver] = useState(false);
  const [showSettings, setShowSettings] = useState(false);
  const fileRef = useRef();

  const runAnalysis = useCallback(async (text) => {
    setLoading(true); setError(""); setTreeRoot(null); setDistMatrix(null); setAiResult(null); setBootstrapVals({});
    await new Promise(r => setTimeout(r, 30));
    try {
      const seqs = parseFASTA(text);
      if (seqs.length < 3) throw new Error("Need at least 3 sequences in FASTA format.");
      const labels = seqs.map(s => s.header);
      const msa = starMSA(seqs);
      const mat = computeDistances(msa, distModel);
      const root = treeMethod === "upgma" ? upgma(mat, labels) : neighborJoining(mat, labels);
      const nwk = toNewick(root);
      setDistMatrix(mat);
      setTreeRoot(root);
      setNewick(nwk);
      setSequences(seqs);
      if (runBootstrap) {
        const bs = bootstrapSupport(msa, labels, distModel, treeMethod, 100);
        setBootstrapVals(bs);
      }
    } catch (e) {
      setError(e.message);
    }
    setLoading(false);
  }, [distModel, treeMethod, runBootstrap]);

  const handleFile = (file) => {
    const reader = new FileReader();
    reader.onload = e => runAnalysis(e.target.result);
    reader.readAsText(file);
  };

  const runAI = async () => {
    if (!newick || !sequences) return;
    setAiLoading(true);
    const result = await analyzeWithAI(newick, sequences.map(s => s.header));
    setAiResult(result);
    setAiLoading(false);
  };

  const handleAnnotate = () => {
    if (!selectedNode || !selectedNode.name) return;
    setAnnotations(prev => ({ ...prev, [selectedNode.name]: { label: annotLabel, note: annotNote } }));
    setAnnotLabel(""); setAnnotNote("");
  };

  const S = {
    card: { background: "rgba(15,23,42,0.75)", border: "1px solid rgba(30,64,175,0.4)", borderRadius: "12px", padding: "18px", backdropFilter: "blur(10px)" },
    label: { fontSize: "10px", color: "#64748b", letterSpacing: "0.12em", marginBottom: "5px" },
    input: { background: "rgba(2,6,23,0.8)", border: "1px solid rgba(30,64,175,0.5)", borderRadius: "6px", color: "#cbd5e1", padding: "6px 10px", fontSize: "12px", width: "100%", fontFamily: "monospace", boxSizing: "border-box" },
    btn: (active) => ({ background: active ? "rgba(6,182,212,0.2)" : "rgba(15,23,42,0.6)", border: `1px solid ${active ? "rgba(6,182,212,0.6)" : "rgba(30,64,175,0.4)"}`, borderRadius: "6px", color: active ? "#67e8f9" : "#94a3b8", padding: "5px 14px", fontSize: "11px", cursor: "pointer", fontFamily: "monospace" }),
    sectionTitle: { fontSize: "12px", color: "#67e8f9", letterSpacing: "0.12em", marginBottom: "12px" },
  };

  return (
    <div style={{ minHeight: "100vh", background: "#020617", color: "#e2e8f0", fontFamily: "'Courier New', monospace", overflowX: "hidden" }}>
      {/* Ambient */}
      <div style={{ position: "fixed", top: -120, left: -120, width: 520, height: 520, borderRadius: "50%", background: "radial-gradient(circle, rgba(6,182,212,0.07) 0%, transparent 70%)", pointerEvents: "none", zIndex: 0 }} />
      <div style={{ position: "fixed", bottom: -100, right: -100, width: 560, height: 560, borderRadius: "50%", background: "radial-gradient(circle, rgba(29,78,216,0.07) 0%, transparent 70%)", pointerEvents: "none", zIndex: 0 }} />

      {/* Header */}
      <header style={{ position: "sticky", top: 0, zIndex: 100, background: "rgba(2,6,23,0.88)", backdropFilter: "blur(12px)", borderBottom: "1px solid rgba(30,64,175,0.3)", padding: "0 28px" }}>
        <div style={{ maxWidth: 1380, margin: "0 auto", display: "flex", alignItems: "center", justifyContent: "space-between", height: 56 }}>
          <div style={{ display: "flex", alignItems: "center", gap: 14 }}>
            <span style={{ fontSize: 22, filter: "drop-shadow(0 0 10px #06b6d4)" }}>🧬</span>
            <div>
              <div style={{ fontSize: 16, fontWeight: "bold", letterSpacing: "0.12em", color: "#67e8f9" }}>PhyloGenAI</div>
              <div style={{ fontSize: 10, color: "#334155", letterSpacing: "0.22em" }}>GENOMIC INTELLIGENCE UNIT</div>
            </div>
          </div>
          <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
            <span style={{ width: 7, height: 7, borderRadius: "50%", background: "#34d399", display: "inline-block", boxShadow: "0 0 7px #34d399", animation: "pulse 2s infinite" }} />
            <span style={{ fontSize: 11, color: "#34d399", letterSpacing: "0.12em" }}>SYSTEM ONLINE</span>
          </div>
        </div>
      </header>

      <div style={{ maxWidth: 1380, margin: "0 auto", padding: "24px 20px", display: "grid", gridTemplateColumns: treeRoot ? "minmax(0,1fr) 340px" : "minmax(0,1fr)", gap: 18, position: "relative", zIndex: 1 }}>

        {/* LEFT */}
        <div style={{ display: "flex", flexDirection: "column", gap: 16 }}>
          {/* Upload Card */}
          <div style={S.card}>
            <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 14 }}>
              <div style={S.sectionTitle}>⬆ SEQUENCE INPUT</div>
              <div style={{ display: "flex", gap: 6 }}>
                <button onClick={() => setShowSettings(s => !s)} style={S.btn(showSettings)}>⚙ Settings</button>
                <button onClick={() => runAnalysis(DEMO_FASTA)} style={S.btn(false)}>🔬 Demo</button>
              </div>
            </div>

            {showSettings && (
              <div style={{ marginBottom: 14, padding: 12, background: "rgba(2,6,23,0.6)", borderRadius: 8, border: "1px solid rgba(30,64,175,0.25)", display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10 }}>
                <div>
                  <div style={S.label}>DISTANCE MODEL</div>
                  <select value={distModel} onChange={e => setDistModel(e.target.value)} style={S.input}>
                    <option value="jaccard">Jaccard k-mer</option>
                    <option value="percent">Percent Identity</option>
                    <option value="poisson">Poisson Corrected</option>
                  </select>
                </div>
                <div>
                  <div style={S.label}>TREE METHOD</div>
                  <select value={treeMethod} onChange={e => setTreeMethod(e.target.value)} style={S.input}>
                    <option value="nj">Neighbor-Joining</option>
                    <option value="upgma">UPGMA</option>
                  </select>
                </div>
                <div style={{ gridColumn: "span 2", display: "flex", alignItems: "center", gap: 8 }}>
                  <input type="checkbox" id="bs" checked={runBootstrap} onChange={e => setRunBootstrap(e.target.checked)} style={{ accentColor: "#06b6d4", cursor: "pointer" }} />
                  <label htmlFor="bs" style={{ fontSize: 11, color: "#94a3b8", cursor: "pointer" }}>Bootstrap Validation (100 replicates) — adds ~5–10 s</label>
                </div>
              </div>
            )}

            <div
              onDragOver={e => { e.preventDefault(); setDragOver(true); }}
              onDragLeave={() => setDragOver(false)}
              onDrop={e => { e.preventDefault(); setDragOver(false); const f = e.dataTransfer.files[0]; if (f) handleFile(f); }}
              onClick={() => fileRef.current.click()}
              style={{ border: `2px dashed ${dragOver ? "#06b6d4" : "rgba(30,64,175,0.45)"}`, borderRadius: 8, padding: "28px 20px", textAlign: "center", cursor: "pointer", transition: "all 0.25s", background: dragOver ? "rgba(6,182,212,0.04)" : "transparent", boxShadow: dragOver ? "0 0 22px rgba(6,182,212,0.12)" : "none" }}>
              <div style={{ fontSize: 30, marginBottom: 8 }}>🧬</div>
              <div style={{ fontSize: 13, color: "#64748b" }}>Drop <span style={{ color: "#67e8f9" }}>.fasta</span> / <span style={{ color: "#67e8f9" }}>.fa</span> here, or click to upload</div>
              <div style={{ fontSize: 11, color: "#1e293b", marginTop: 4 }}>Minimum 3 sequences required</div>
            </div>
            <input ref={fileRef} type="file" accept=".fasta,.fa,.txt" style={{ display: "none" }} onChange={e => { if (e.target.files[0]) handleFile(e.target.files[0]); }} />

            {error && (
              <div style={{ marginTop: 10, padding: "8px 12px", background: "rgba(239,68,68,0.1)", border: "1px solid rgba(239,68,68,0.3)", borderRadius: 6, color: "#f87171", fontSize: 12 }}>
                ⚠ {error}
              </div>
            )}
          </div>

          {/* Tree Card */}
          {(treeRoot || loading) && (
            <div style={S.card}>
              <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 14, flexWrap: "wrap", gap: 8 }}>
                <div style={S.sectionTitle}>🌿 PHYLOGENETIC TREE</div>
                <div style={{ display: "flex", gap: 6 }}>
                  <button onClick={() => setCircular(c => !c)} style={S.btn(circular)}>◎ Radial</button>
                  <button onClick={() => setPhylogram(p => !p)} style={S.btn(phylogram)}>⟿ Phylogram</button>
                </div>
              </div>
              {loading ? (
                <div style={{ display: "flex", alignItems: "center", justifyContent: "center", height: 200, gap: 12 }}>
                  <div style={{ width: 22, height: 22, border: "2px solid rgba(6,182,212,0.25)", borderTop: "2px solid #06b6d4", borderRadius: "50%", animation: "spin 0.9s linear infinite" }} />
                  <span style={{ color: "#475569", fontSize: 12 }}>Computing phylogeny…</span>
                </div>
              ) : (
                <TreeVisualizer
                  treeData={treeRoot}
                  bootstrapValues={bootstrapVals}
                  circular={circular}
                  phylogram={phylogram}
                  onNodeClick={node => { setSelectedNode(node); setActiveTab("annotate"); }}
                  annotations={annotations}
                />
              )}
              {runBootstrap && Object.keys(bootstrapVals).length > 0 && (
                <div style={{ display: "flex", gap: 14, fontSize: 10, marginTop: 8 }}>
                  <span style={{ color: "#34d399" }}>■ ≥90%</span>
                  <span style={{ color: "#fbbf24" }}>■ ≥70%</span>
                  <span style={{ color: "#f87171" }}>■ &lt;70%</span>
                  <span style={{ color: "#475569" }}>bootstrap support</span>
                </div>
              )}
            </div>
          )}
        </div>

        {/* RIGHT SIDEBAR */}
        {treeRoot && (
          <div style={{ display: "flex", flexDirection: "column", gap: 14 }}>
            {/* Tab Bar */}
            <div style={{ display: "flex", gap: 4, background: "rgba(15,23,42,0.7)", border: "1px solid rgba(30,64,175,0.35)", borderRadius: 10, padding: 5 }}>
              {[["ai", "🤖 AI"], ["annotate", "📌 Annotate"], ["matrix", "📊 Matrix"]].map(([k, lbl]) => (
                <button key={k} onClick={() => setActiveTab(k)} style={{ ...S.btn(activeTab === k), flex: 1, fontSize: 10, padding: "6px 2px" }}>{lbl}</button>
              ))}
            </div>

            {/* AI Tab */}
            {activeTab === "ai" && (
              <div style={S.card}>
                <div style={S.sectionTitle}>🤖 EVOLUTIONARY INSIGHTS</div>
                {!aiResult && !aiLoading && (
                  <button onClick={runAI} style={{ ...S.btn(false), width: "100%", padding: "10px 0", fontSize: 12 }}>
                    ✨ Analyze with AI
                  </button>
                )}
                {aiLoading && (
                  <div style={{ display: "flex", alignItems: "center", gap: 10, color: "#475569", fontSize: 12 }}>
                    <div style={{ width: 15, height: 15, border: "2px solid rgba(6,182,212,0.3)", borderTop: "2px solid #06b6d4", borderRadius: "50%", animation: "spin 0.9s linear infinite" }} />
                    Consulting evolutionary database…
                  </div>
                )}
                {aiResult && (
                  <div style={{ display: "flex", flexDirection: "column", gap: 12 }}>
                    <div style={{ padding: 12, background: "rgba(6,182,212,0.05)", borderRadius: 8, border: "1px solid rgba(6,182,212,0.18)" }}>
                      <div style={{ ...S.label, marginBottom: 6 }}>SUMMARY</div>
                      <p style={{ fontSize: 12, color: "#94a3b8", lineHeight: 1.65, margin: 0 }}>{aiResult.summary}</p>
                    </div>
                    {aiResult.clades?.length > 0 && (
                      <div>
                        <div style={S.label}>DETECTED CLADES</div>
                        {aiResult.clades.map((c, i) => (
                          <div key={i} style={{ marginBottom: 7, padding: "9px 11px", background: "rgba(29,78,216,0.09)", borderRadius: 7, borderLeft: "2px solid #1d4ed8" }}>
                            <div style={{ fontSize: 11, color: "#60a5fa", marginBottom: 3 }}>{c.name}</div>
                            <div style={{ fontSize: 10, color: "#64748b", marginBottom: 3 }}>{c.members?.join(" · ")}</div>
                            <div style={{ fontSize: 11, color: "#94a3b8" }}>{c.interpretation}</div>
                          </div>
                        ))}
                      </div>
                    )}
                    {aiResult.keyFindings?.length > 0 && (
                      <div>
                        <div style={S.label}>KEY FINDINGS</div>
                        {aiResult.keyFindings.map((f, i) => (
                          <div key={i} style={{ display: "flex", gap: 8, marginBottom: 5, fontSize: 11, color: "#94a3b8" }}>
                            <span style={{ color: "#06b6d4", flexShrink: 0 }}>›</span>
                            <span>{f}</span>
                          </div>
                        ))}
                      </div>
                    )}
                    <button onClick={() => setAiResult(null)} style={{ ...S.btn(false), fontSize: 11 }}>↻ Re-analyze</button>
                  </div>
                )}
              </div>
            )}

            {/* Annotate Tab */}
            {activeTab === "annotate" && (
              <div style={S.card}>
                <div style={S.sectionTitle}>📌 NODE ANNOTATIONS</div>
                {selectedNode ? (
                  <div style={{ display: "flex", flexDirection: "column", gap: 10 }}>
                    <div style={{ padding: "7px 11px", background: "rgba(168,85,247,0.1)", border: "1px solid rgba(168,85,247,0.3)", borderRadius: 6, fontSize: 12, color: "#c084fc" }}>
                      Selected: <strong>{selectedNode.name || "(internal node)"}</strong>
                    </div>
                    <div>
                      <div style={S.label}>LABEL</div>
                      <input value={annotLabel} onChange={e => setAnnotLabel(e.target.value)} placeholder="e.g. Outgroup" style={S.input} />
                    </div>
                    <div>
                      <div style={S.label}>NOTES</div>
                      <textarea value={annotNote} onChange={e => setAnnotNote(e.target.value)} placeholder="Biological notes..." rows={3} style={{ ...S.input, resize: "vertical" }} />
                    </div>
                    <button onClick={handleAnnotate} style={{ ...S.btn(true), width: "100%", padding: "8px 0" }}>+ Save Annotation</button>
                    {annotations[selectedNode.name] && (
                      <div style={{ padding: "7px 10px", background: "rgba(168,85,247,0.06)", borderRadius: 6, border: "1px solid rgba(168,85,247,0.2)", fontSize: 11, color: "#a78bfa" }}>
                        ✓ Node glows purple on tree
                      </div>
                    )}
                  </div>
                ) : (
                  <div style={{ padding: "20px 0", textAlign: "center", color: "#334155", fontSize: 12 }}>
                    Click any node on the tree to annotate it
                  </div>
                )}
                {Object.keys(annotations).length > 0 && (
                  <div style={{ marginTop: 14 }}>
                    <div style={S.label}>SAVED ANNOTATIONS</div>
                    {Object.entries(annotations).map(([name, ann]) => (
                      <div key={name} style={{ marginBottom: 6, padding: "8px 10px", background: "rgba(2,6,23,0.5)", borderRadius: 6, border: "1px solid rgba(168,85,247,0.2)" }}>
                        <div style={{ fontSize: 11, color: "#c084fc" }}>{name}</div>
                        {ann.label && <div style={{ fontSize: 10, color: "#64748b", marginTop: 2 }}>{ann.label}</div>}
                        {ann.note && <div style={{ fontSize: 11, color: "#94a3b8", marginTop: 3 }}>{ann.note}</div>}
                      </div>
                    ))}
                  </div>
                )}
              </div>
            )}

            {/* Matrix Tab */}
            {activeTab === "matrix" && (
              <div style={S.card}>
                <div style={S.sectionTitle}>📊 DISTANCE MATRIX</div>
                {distMatrix && sequences ? (
                  <DistanceMatrix matrix={distMatrix} labels={sequences.map(s => s.header)} />
                ) : (
                  <div style={{ color: "#334155", fontSize: 12, textAlign: "center", padding: "20px 0" }}>No matrix yet</div>
                )}
              </div>
            )}

            {/* Newick */}
            <div style={{ ...S.card, padding: 14 }}>
              <div style={{ ...S.label, marginBottom: 7 }}>NEWICK STRING</div>
              <div style={{ background: "rgba(2,6,23,0.9)", borderRadius: 6, padding: "8px 10px", fontSize: 10, color: "#334155", fontFamily: "monospace", wordBreak: "break-all", maxHeight: 80, overflowY: "auto", lineHeight: 1.5 }}>
                {newick || "—"}
              </div>
            </div>
          </div>
        )}
      </div>

      <style>{`
        @keyframes spin { to { transform: rotate(360deg); } }
        @keyframes pulse { 0%,100%{opacity:1} 50%{opacity:0.4} }
        * { box-sizing: border-box; }
        ::-webkit-scrollbar { width: 4px; height: 4px; }
        ::-webkit-scrollbar-track { background: rgba(15,23,42,0.4); }
        ::-webkit-scrollbar-thumb { background: rgba(30,64,175,0.5); border-radius: 2px; }
        select option { background: #0f172a; color: #cbd5e1; }
      `}</style>
    </div>
  );
}
