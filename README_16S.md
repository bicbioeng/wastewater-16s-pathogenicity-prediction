<img alt="Genome-Resolved 16S Pathogenicity Prediction" src="assets/pipeline_overview.png" />

# Wastewater 16S Pathogenicity Prediction

Genome-resolved prediction of bacterial pathogenicity from **16S rRNA genes**, using a
**two-model ensemble** (identify + classify) with a precision-biased fusion layer. This is
the 16S counterpart of the whole-genome MAG pipeline
[`wastewater-mag-pathogenicity-prediction`](https://github.com/bicbioeng/wastewater-mag-pathogenicity-prediction),
built within the GOMICS framework. Every call is tied back to a reconstructed genome (MAG),
so a flagged organism can be followed into its source genome for downstream analysis.

| | Question | Method |
|---|---|---|
| **Model A** — identify | *Which organism is this?* | BLAST nearest-neighbour vs a **26,877**-sequence 16S RefSeq reference → taxonomy + rank |
| **Model B** — classify | *Is it pathogenic?* | canonical **6-mer (2,080)** vector → class-balanced **Random Forest** → P(pathogen) |
| **Fusion** — decide | *One honest call* | precision-biased rule layer with a reject option → PATHOGEN / OPPORTUNIST / NON-PATHOGEN / REVIEW / UNKNOWN |

**Model performance (ground truth, 1,432 labeled 16S genes across 127 genera):**
leave-one-genus-out AUROC **0.81** (accuracy 0.71) · in-distribution AUROC **0.99** (accuracy 0.99).

---

## 📦 Getting Started

This section explains how to clone the repository, fetch the required datasets and models,
and run the analysis notebook.

---

## 🔁 Clone the Repository

```bash
git clone https://github.com/bicbioeng/wastewater-16s-pathogenicity-prediction.git
cd wastewater-16s-pathogenicity-prediction
```

---

## 📂 Dataset & Model Setup

Large inputs (the Model A reference, the trained Model B, and the extracted 16S) are hosted
externally to keep the repository lightweight. Place them under `datasets/` and
`saved_model_16s/` as shown below.

### 1️⃣ Create the required directories

```bash
mkdir -p datasets/05_16s_from_mags datasets/06_16s_pathogen_calls saved_model_16s/ref_16s
```

### 2️⃣ Download the trained model + labels (shared Drive)

Download into `saved_model_16s/`:

```bash
# from the shared Drive folder (or gdown / rclone), place these files in saved_model_16s/:
#   kmer16s_model.pkl              trained Model B (Random Forest)
#   pathogen_16s_catalog.json      fusion catalog
#   genome_16s_pathogen.fasta      826 pathogen 16S genes (Model B labels)
#   genome_16s_nonpathogen.fasta   606 non-pathogen 16S genes
# example with gdown:
# gdown --folder <shared-drive-folder-url> -O saved_model_16s/
```

### 3️⃣ Download the Model A reference (NCBI 16S RefSeq, ~26,877)

```bash
wget -O saved_model_16s/ref_16s/16S_ribosomal_RNA.tar.gz \
  https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
tar -xzf saved_model_16s/ref_16s/16S_ribosomal_RNA.tar.gz -C saved_model_16s/ref_16s
```

### 4️⃣ Download the genome-resolved 16S (this study)

```bash
# place all_16s.fna (≈1,321 representative 16S) into:
#   datasets/05_16s_from_mags/all_16s.fna
```

### ✅ Expected directory structure

```
wastewater-16s-pathogenicity-prediction/
├── assets/
│   └── pipeline_overview.png
├── wastewater_16s_pathogenicity_prediction.ipynb   Objective-per-cell runnable notebook
├── extract_16s.py                                  barrnap + GFF → genome-resolved 16S
├── pathogen_predict.py                             unified inference tool (`16s` mode)
├── train_16s_model.py                              k-mer featurization + Random Forest
├── validate_16s_model.py                           ground-truth validation harness
├── kmer_feature_importance_16s.py                  top canonical k-mers
├── feature_importance_genes_go.py                  (bridge) pangenome feature → gene → GO
├── go_pathway_prediction.py                        (bridge) GO enrichment
├── saved_model_16s/
│   ├── kmer16s_model.pkl
│   ├── pathogen_16s_catalog.json
│   ├── genome_16s_pathogen.fasta
│   ├── genome_16s_nonpathogen.fasta
│   └── ref_16s/                                    Model A BLAST DB (16S RefSeq)
└── datasets/
    ├── 05_16s_from_mags/
    │   └── all_16s.fna
    └── 06_16s_pathogen_calls/                      inference output
```

---

## 🧠 Running the Analysis

Ensure Jupyter Notebook or JupyterLab is installed. The full pipeline runs objective by
objective in:

```
wastewater_16s_pathogenicity_prediction.ipynb
```

Launch it with:

```bash
jupyter notebook wastewater_16s_pathogenicity_prediction.ipynb
```

The notebook covers, one cell per objective:

1. Build the Model A reference (16S RefSeq → BLAST DB)
2. Build the Model B labeled set (826 pathogen + 606 non-pathogen)
3. Train Model B (canonical 6-mer Random Forest + fragment augmentation)
4. Extract genome-resolved 16S from the wastewater MAGs
5. Inference — identify → classify → fuse
6. **Validate the model against ground truth** (stratified + leave-one-genus-out CV)
7. k-mer feature importance (what drives Model B)
8. Community composition, diversity & phylogeny
9. **Genome-resolved bridge** — pangenome + GO enrichment on flagged pathogens

### 🔬 Model validation (accuracy)

```bash
python validate_16s_model.py \
  --labeled-fasta saved_model_16s/labeled_16s_full.fasta \
  --kmer-k 6 --folds 5 --augment 2 --min-genus-size 2 \
  --outdir model_b_validation
```

Reports **stratified k-fold** (in-distribution) and **leave-one-genus-out**
(generalization) accuracy, plus ROC / confusion / calibration figures. Report the
leave-one-genus-out numbers as the model's accuracy.

### 🧬 Inference on new 16S

```bash
python pathogen_predict.py 16s --input datasets/05_16s_from_mags/all_16s.fna \
  --output datasets/06_16s_pathogen_calls
```

---

## ⚙️ Requirements

**Python (≥ 3.10):** `pandas numpy scikit-learn matplotlib seaborn scipy statsmodels goatools requests tqdm biopython nbformat`

**Bioinformatics tools (conda / bioconda):**

```bash
conda create -n gomics_16s python=3.11 -y && conda activate gomics_16s
pip install pandas numpy scikit-learn matplotlib seaborn scipy statsmodels goatools requests tqdm biopython nbformat
conda install -c bioconda blast barrnap ppanggolin mafft fasttree -y
```

`blast` + `barrnap` cover the core 16S workflow; `ppanggolin`, `mafft`, `fasttree` are used
only by the phylogeny (Obj. 8) and the genome-resolved bridge (Obj. 9).

---

## 📝 Notes

- All notebooks assume the directory structure shown above; paths are relative to the project root.
- Large models and datasets are hosted externally (shared Drive / NCBI) to keep the repository lightweight.
- **Model accuracy** should be reported from the **leave-one-genus-out** cross-validation (leakage-controlled); the in-distribution figure is an optimistic bound.
- The 16S method ships as a subcommand of the same CLI as the whole-genome modes — identical UX.

---

## 📧 Contact

BicBioEng / GEAR Center, University of South Dakota. For questions or issues, please open a
GitHub issue or contact the research team.
