Figure_1 — Scripts for initial data processing
Figure_2 — Cross‑tissue DC subtype conservation analyses
Figure_3 — DC developmental gene‑regulatory programs
Figure_4 — scVelo‑inferred pseudotime trajectories
Figure_5 — DC–disease associations   
Figure_6 — SLE data analyses

**Operating System:** Tested on standard macOS and Linux computing environments.

**Software:**
* R (4.3.0)
* Python (3.11.0)
* Cell Ranger (2.0.2)
* PLINK (1.9)
* FlowJo (v10.2) - *Used for Flow Cytometry data analysis*

**Python Packages:**
* scVelo (0.3)
* SCENIC+ (1.0)
* velocyto (0.17.17)
* CellHint (1.0)
* CellTypist (1.7.1)

**R Packages:**
* Seurat (5.4)
* Signac (1.16)
* MiloR (1.12)
* sccomp (3.22)
* scE2G (1.0)
* Pando (1.0)
* CytoTRACE2 (1.1)
* gchromVAR (1.0)
* SCAVENGE (1.0)
* regioneR (1.44)
* Rsubread (2.26)
* qs2

# Demo
library(qs2)
library(Seurat)
DC_demo <- qs2::qs_read("DC_demo.qs2")
