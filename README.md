[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23289539-blue)](https://doi.org/10.5281/zenodo.23289539) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.23289650-blue)](https://doi.org/10.5281/zenodo.23289650) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)


# UPDSF_v6_0_1

Environment: 🐍 Python 3.8+

## Author: Seyed Mohammad Reza Hashemi (Reza Hashemi)

#### ORCID: 0009-0002-0645-5180


## UNIFIED PREBIOTIC DNA SELECTION FRAMEWORK (UPDSF) v6.0.1 Chemical Filtering and T-over-U Enrichment Model

### The UPDSF v6.0.1 model demonstrates chemical filtering of T over U under modeled prebiotic conditions, peaking at T = 90.0 °C and pH 11.0. The total CSI of 2.624× is the result of balanced contributions from both direct environmental discrimination and indirect pool drift.
The model is robust across multiple orders of magnitude of UV flux, with hydrolysis remaining the dominant selective force under alkaline hydrothermal conditions. It is critical to note that this model demonstrates chemical filtering, not Darwinian evolution. However, these results suggest that in alkaline hydrothermal scenarios, differential nucleobase stability provided a significant selective advantage for thymine incorporation into early genetic polymers, potentially paving the way for the transition from RNA to DNA.


---
## 🛠️ Parameter Provenance & Model Transparency

To ensure the scientific integrity and reproducibility of the *UPDSF v6.0.1 model, we employ a rigorous Provenance-Tracking System for all model parameters. Instead of using "hard-coded" constants, every numerical value is encapsulated in a `Source` object.

<br>
<br>

📂 **Complete Provenance-Tracked Parameter List for UPDSF v6.0.1 Model**


The Complete Provenance-Tracked Parameter List for UPDSF v6.0.1 Model can be found in the 📄 [Detailed Parameter List](data/Parameter_List_for_UPDSF_v6.0.1_Model.docx)

<br>
<br>

### Why a Detailed Parameter List?

The use of a structured Parameter List (available in `Parameter List for UPDSF v6.0.1 Model.docx`) serves four critical purposes:

1. Scientific Reproducibility: By explicitly listing all 95 core parameters, any researcher can replicate the exact simulation environment and verify the results without ambiguity.
2. Source Accountability (Provenance): Each parameter is tagged with a specific Kind, allowing for a clear distinction between:
     📚 *Literature: Values derived from peer-reviewed experimental data (e.g., Lindahl 1993, Ferris 1996).
     💡 *Assumed: Heuristics based on chemical intuition where direct data is unavailable.
     ⚙️ *Imposed: Modeling choices necessary for framework stability or specific scenario testing.
     📈 *Derived: Values computed mathematically from other reference constants.
3. Sensitivity Analysis Baseline:* The detailed list provides a baseline for "What-if" scenarios. By knowing exactly which parameters are assumed versus literature-based, we can pinpoint the model's vulnerabilities and perform targeted sensitivity sweeps.
4. *Transparency in "Chemical Filtering":* It allows us to isolate the exact physical drivers (e.g., the $\Delta E_a$ between Thymine and Uracil) that lead to the observed Chemical Selection Index (CSI), proving that the results are a consequence of fundamental kinetics rather than arbitrary tuning.

---
This commitment to transparency transforms the model from a "black box" into an open scientific framework for prebiotic chemistry research.
---

---
## 🧪 Note on Parameter Distribution & Model Robustness

A transparent review of our provenance tracking reveals that approximately 46% of the model parameters are based on direct literature values, while the remainder are assumed, imposed, or derived. In many traditional models, this might be seen as a limitation; however, in the context of Prebiotic Chemistry, this distribution is a deliberate and necessary architectural choice for the following reasons:

### 1. The "Data Gap" in Prebiotic Conditions
Experimental data for nucleobase stability at extreme conditions (e.g., $T = 90^\circ\text{C}$ and $\text{pH} = 11.0$) is sparse or non-existent in current literature. To bridge this gap, we employ physicochemical heuristics—educated assumptions based on known chemical trends—to allow the model to explore plausible prebiotic scenarios that have not yet been tested in a wet lab.

### 2. Heuristics as Hypotheses
The `Assumed` and `Imposed` parameters are not "guesses"; they are formalized hypotheses*. By explicitly tagging them as such, we distinguish between known chemical facts and model-driven predictions. This transparency allows other researchers to easily adjust these specific values to test the boundaries of our conclusions.

### 3. Robustness through Sensitivity Analysis
To compensate for the reliance on assumed values, we conducted extensive *Global Sensitivity Sweeps. Our results demonstrate that the enrichment of Thymine over Uracil is not dependent on the fine-tuning of a single assumed parameter, but emerges from the fundamental kinetic gap (e.g., $\Delta E_a$) provided by the literature-based constants. 

### 4. A Framework for Future Experimentation
We view the `Assumed` parameters as targets for future experimental validation*. The UPDSF v6.0.1 model identifies exactly which "unknowns" (e.g., specific deamination rates of free monomers) most significantly impact the results, thereby providing a roadmap for experimentalists to prioritize their lab work.

---
By disclosing the exact ratio of literature-to-assumed parameters, we move away from "black-box" simulations toward a verifiable, iterative scientific process.
---

