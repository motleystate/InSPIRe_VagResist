# InSPIRe Vaginal Resistome Paper 

### Content of this subfolder 
```
─ InSPIRe_VaginalResistome_paper/   : Folder with main article downstream analyses (notebooks + data).  
   └─ data/             :   Includes all raw and intermediate data required to replicate the main results presented in this study.
   │  ├─ supp_materials :   Dictionaries, trees and gene lengths for normalization.
   │  ├─ Norm/ :            Normalized count tables for species, ARGs and MGEs.
   │  ├─ Output/ :          Intermediate analysis outputs (e.g. SparCC results) used in notebooks.
   │  ├─ Raw/:              Raw count tables for ARGs and MGEs. 
   │  inspire_VR_metadata/: InSPIRe cohort metadata with variables used in this study.
   │  ExternalDataset_CST/: CST categorization for the external dataset.
   └─ Notebooks/:           Notebooks presenting the main analysis of the manuscript (figures and statistical analysis).
``` 

| Notebook | Information |
|----------|-------------|
| `InSPIRe_Analysis1_ProcessingData` | Data processing and normalization. |
| `InSPIRe_Analysis2_DiversityAnalysis` | Global distributions of ARGs and MGEs and diversity analyses across CSTs and clinical outcomes. |
| `InSPIRe_Analysis3_StatisticalTesting` | Statistical analyses of the main results presented in the manuscript. |
| `InSPIRe_Analysis_Supplementary` | Supplementary analyses, including PRDI conceptual validation. |


### Dependancies 
To reproduce notebooks analyses, you must have Python and R installed with the following libraries:

**for python:**
* `Pandas` 
* `numpy` 
* `matplotlib` and `seaborn` 
* `scikit-bio `
* `json` 

**for R:**  
* `tidyverse`
* `tibble`
* `dplyr`
* `reporttools`
* `dunn.test`
* `glmmTMB`
* `emmeans`
* `DHARMa`
* `performance`
* `see`
* `vegan`
* `ggplot2`
* `Maaslin2`

### Processing data before running notebooks 
Data need to be extracted before running analysis of the notebooks   
```bash
unzip /data.zip 
```