# 🧬 PhyloGenAI — Genomic Intelligence Unit

> A fully client-side, AI-powered phylogenetic tree generator and evolutionary analysis tool. Upload raw biological sequences (FASTA), align them, compute genetic distances, build phylogenetic trees, and get deep evolutionary insights powered by AI.

---

## ✨ Features

### 🔬 Bioinformatics Engine (100% Client-Side)
- **FASTA Parser** — robust multi-line sequence parsing with edge-case handling
- **Needleman-Wunsch** — global pairwise alignment with BLOSUM62 scoring matrix
- **Star MSA** — multiple sequence alignment via Jaccard k-mer center selection
- **3 Distance Models** — Jaccard k-mer, Percent Identity, Poisson-corrected
- **2 Tree Methods** — Neighbor-Joining (NJ) and UPGMA
- **Bootstrap Resampling** — 100-replicate confidence scoring for internal clades

### 🌿 Interactive Tree Visualization (D3.js)
- Linear and Radial (circular) layouts
- Cladogram vs. Phylogram (branch-length scaled) modes
- Hover glow effects and click-to-select nodes
- Bootstrap support values color-coded (≥90% green, ≥70% amber, <70% red)
- SVG export with dark-mode background

### 🤖 AI-Powered Analysis
- Sends Newick tree string + taxa to Claude AI
- Returns structured evolutionary summary, detected clades, and key findings
- Biological interpretation by an AI acting as an evolutionary biologist

### 🖥️ UI Features
- Drag & drop `.fasta` / `.fa` file upload
- Click any node to annotate it with labels and notes (turns purple with glow)
- Color-coded distance matrix panel
- Demo data: Human, Chimp, Gorilla, Mouse, Rat hemoglobin sequences
- Deep-tech dark UI with glassmorphism and cyan glow effects

---

## 🚀 Getting Started

### Prerequisites
- Node.js 18+ 
- npm

### Installation

```bash
# 1. Clone the repository
git clone https://github.com/your-username/phylogenai.git
cd phylogenai

# 2. Install dependencies
npm install

# 3. Start development server
npm run dev
```

Open [http://localhost:5173](http://localhost:5173) in your browser.

### Build for Production

```bash
npm run build
npm run preview
```

---

## 📁 File Structure

```
phylogenai/
├── index.html          ← App entry point
├── main.jsx            ← React root mount
├── App.jsx             ← Root component
├── PhyloGenAI.jsx      ← Main application (all logic + UI)
├── index.css           ← Global styles + Tailwind
├── package.json        ← Dependencies
├── vite.config.js      ← Vite configuration
├── tailwind.config.js  ← Tailwind configuration
├── postcss.config.js   ← PostCSS configuration
└── .gitignore
```

---

## 🔬 How to Use

1. **Upload a FASTA file** — drag & drop or click to browse (min. 3 sequences)
2. **Or click Demo** — loads 5 mammalian hemoglobin sequences instantly
3. **Configure settings** — choose distance model and tree building method
4. **View the tree** — toggle between linear/radial and cladogram/phylogram
5. **Click a node** — annotate it with custom labels and notes
6. **Run AI Analysis** — get evolutionary insights from Claude AI
7. **Check Matrix tab** — view the color-coded pairwise distance matrix
8. **Export** — download the tree as an SVG file

---

## 🧮 Data Flow

```
FASTA Input
    ↓
parseFASTA()
    ↓
starMSA()          ← Needleman-Wunsch + BLOSUM62
    ↓
computeDistances() ← Jaccard / Percent Identity / Poisson
    ↓
neighborJoining()  ← or upgma()
    ↓
bootstrapSupport() ← optional 100 replicates
    ↓
toNewick() ──────→ Claude AI Analysis
    ↓
TreeVisualizer()   ← D3.js SVG rendering
```

---

## 🛠️ Tech Stack

| Technology | Purpose |
|---|---|
| React 18 | UI framework |
| D3.js v7 | Tree visualization |
| Vite 5 | Build tool |
| Tailwind CSS | Styling |
| Claude API | AI evolutionary analysis |

---

## ☁️ Deployment

This project is deployed on **Vercel**. To deploy your own:

1. Push all files to a GitHub repository
2. Go to [vercel.com](https://vercel.com) and import the repo
3. Vercel auto-detects Vite — no extra configuration needed
4. Click **Deploy**

Live demo: [phylogen-ai.vercel.app](https://phylogen-ai.vercel.app)

---

## 📜 License

MIT © 2026 PhyloGenAI
