Repository for analysis scripts and processed data for the paper:
"Regional endemism instead of unlimited gene flow: Phylogeography of nivicolous myxomycetes in the mountains of Eurasia" (Ecosphere; manuscript ID ECS25-0940, in revision)

Raw and processed data

- Raw data (sequences, sample metadata, environmental layers):
  - Published on Zenodo: https://doi.org/10.5281/zenodo.17447098
  - Large files (e.g. rasters) are hosted on Zenodo and not included in this repository.

- Processed input files included in this repository:
  - _data/processed_data/
    - singleton176_mafft.fas
    - unique473_mafft.fas
  - These are the files used for analyses with the included Python scripts.

Reproduce the analyses
1. Download the raw data archive from Zenodo (https://doi.org/10.5281/zenodo.17447098) and extract it. The archive contains `01_Raw_Data/`, `02_Processed_Data/`, and `03_Environmental_Data/`.
2. Open R in the project root (where `renv.lock` is) and restore environment:
   - `install.packages("renv")` (if needed)
   - `renv::restore()`
3. Run scripts in order (exact filenames):
   - `_script/01_data_cleaning_processing.R`
   - `_script/02_basic_statistics.R`
   - `_script/03_geographic_differentiation.R`
   - `_script/04_haplotype_beta_diversity_analysis.R`
   - `_script/05_species_ecological_analysis.R`
   - `_script/SA_sensitivity_analysis.R` (sensitivity analyses; can be run after script 02)
   - Utility functions: `_script/utils_ecology_myxo.R`
4. Run these scripts at any time to analyse the sequence data with Python:
   - `_script/Double peak detection.py`
   - `_script/Single_substitution_detect_among_singleton.py`

Outputs from the sensitivity analyses are saved in `_results/`:
- `SA_GAP_subsampling_summary.csv`
- `SA_Japan_TNS_endemicity.csv`

How to cite
- Paper: to be added upon publication.
- Repository: https://doi.org/10.5281/zenodo.17593131 (this links to all versions; for the specific version, see the Zenodo record)

License
- Code: MIT (see LICENSE)
- Data: CC-BY-4.0 (see DATA_LICENSE)

Zenodo / release notes
- Create a GitHub release (tag) and Zenodo will mint a new version DOI. Example:
  ```
  git tag -a v1.1.0 -m "Revision update for Ecosphere submission"
  git push origin v1.1.0
  gh release create v1.1.0 --title "v1.1.0" --notes "Revision update: sensitivity analyses added"
  ```
