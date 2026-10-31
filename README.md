# Evaluating the robustness of k-mer-based ANI estimation to regional variation in mutation rates

Code, data lists and results for the CMPS 355 (Bioinformatics Algorithms) research paper by
Faysal Elestwani and Nanar Aintablian, American University of Beirut, Fall 2026–2027.

All notebooks were run on Google Colab (CPU runtime). To rerun them, put this folder in Google Drive as
`MyDrive/CMPS355_ANI_Project` and open the notebooks in numbered order. Each notebook installs the tools it needs.

## Notebooks

| Notebook | What it does | Paper section |
|---|---|---|
| `1_Pilot_and_RealGenomes.ipynb` | Pilot simulations on a random 1 Mb sequence; first comparison on real genomes | Methods 2.1, 2.2 (pilot) |
| `2_RealGenome_Estimators.ipynb` | Real-genome estimators, development and held-out pairs, the four aligners, frozen predictions for G | Methods 2.1, 2.3, 2.4 |
| `3_SyntheticPattern_Simulations.ipynb` | Synthetic-pattern simulations (1,260 pairs; "Experiment 6" in the code) | Methods 2.2 |
| `4_Final_Tables.ipynb` | Two-class estimator, AICc model choice, Slope-SpaM, final tables | Methods 2.3, 2.5 |
| `5_RealProfiles_WindowEstimator.ipynb` | WT estimator, real rate profiles, real-profile simulations (360 pairs; "Experiment 7" in the code), sketches | Methods 2.2, 2.3 |

## Folders

- `data/` – genome lists (`v2_genomes.csv`), genome pairs (`v2_pairs.csv`) and reference ANI (`v2_reference_ani.csv`).
  The genome FASTA files are **not** included because of their size. They are public NCBI RefSeq assemblies and can be
  downloaded by accession, for example with the NCBI Datasets tool:
  `datasets download genome accession GCF_000597845.1 --include genome`
- `results/realsim/` – synthetic-pattern simulations: raw results, per-pair k-mer counts (`curves/`), the simulation script `realsim.py`
- `results/exp7/` – real-profile simulations: raw results, the 24 rate profiles, the script `realsim7.py`
- `results/Part2/`, `results/multitool/` – real-genome results (development and held-out pairs, all aligners)
- `results/final/` – final tables and figures used in the paper, and `key_numbers.txt` / `key_numbers_v2.txt`

## Reproducibility

Simulated pair *j* uses random seed 10000 + *j* (synthetic-pattern simulations) or 20000 + *j* (real-profile simulations).
Held-out real pairs were drawn with seed 2026. Tool versions are listed in Section 2.5 of the paper.
