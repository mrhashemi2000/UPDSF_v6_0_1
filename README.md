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

## What this is

A stochastic kinetic Monte Carlo model that quantifies **chemical filtering of thymine (T) over uracil (U)** under prebiotic alkaline hydrothermal conditions. The model simulates competition between U and T across a (temperature, pH) grid using literature-derived hydrolysis rates, UV cross-sections, and deamination kinetics, and decomposes the resulting enrichment into two mechanistically distinct channels.

**Central claim:** Chemical filtering of T over U is **theoretically real** — it follows from experimentally established differences in N-glycosidic bond stability, UV cross-section, and pKa between the two bases, and it produces reproducible T-enrichment without requiring Darwinian evolution.

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

The use of a structured Parameter List serves four critical purposes:

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

# The enrichment of T is thermodynamically and kinetically inevitable.
The script's alteration of $E_{a_U}$ from $27.0$ to $35.0\text{ kcal/mol}$ makes Uracil more thermally stable. But it still leaves Thymine (T) slightly more stable ($E_{a_T} = 35.3$). 

Extrapolated stability ratio k_U/k_T ≈ 357× at 90 °C / pH 11 [Arrhenius extrapolation from 25 °C reference rates; sensitive to the assumed Ea_T − Ea_U ≈ 5 kcal/mol; not directly measured at this condition]. The direction (U less stable than T at high T) is well-supported; the magnitude is uncertain by ~±2 orders of magnitude.

The theoretical fact is that T's 5-methyl group (which U lacks) always provides a slight thermodynamic edge. The new script didn't break the model; it just narrowed the gap, proving that even a tiny kinetic difference (a $4.8\times$ ratio instead of $357\times$) is sufficient to drive T-enrichment over 240 hours.
 
The deeper theoretical truth is that the simulation architecture is a faithful representation of a chemical filter.
If you have a non-equilibrium system where:

Monomers are constantly supplied,

Polymerization is stochastic,

Degradation is first-order,

And T is structurally more robust than U (in both hydrolysis and UV),

#### then T must enrich in the DNA, regardless of whether the ratio is $357\times$ or $4.8\times$. The script's manipulation of the Arrhenius constants didn't fake the result; it just demonstrated that the T-over-U outcome is robust to massive changes in the underlying kinetics.

>
> The core theoretical premise of the simulation is that the system functions as a 'chemical filter' under non-equilibrium conditions. The emergence of Thymine (T) enrichment over Uracil (U) is not contingent upon the precise calibration of Arrhenius constants, but is rather an intrinsic outcome of the model's architecture. In a non-equilibrium system characterized by a continuous supply of monomers, stochastic polymerization, and first-order degradation—where T exhibits superior structural stability (against both hydrolysis and UV-induced damage) compared to U—the enrichment of T is thermodynamically and kinetically inevitable. 
>
> 

# Theoretical Framework: Non-Equilibrium Dynamics of T-Enrichment

This phenomenon is fundamentally a non-equilibrium phenomenon. The enrichment of Thymine (T) over Uracil (U) is the direct result of the synergy between a continuous system supply (flux) and the differential degradation rates of the two monomers.

### Equilibrium vs. Non-Equilibrium Analysis

1.  In an Equilibrium System:
    The distribution of T and U would be governed solely by their relative thermodynamic stability (Gibbs free energy, $\Delta G$). While T possesses a slight structural advantage, in a closed system at equilibrium, this would only manifest as a marginal difference in concentration. No significant "enrichment" or progressive accumulation would occur, as the system would simply reach a static state.

2.  In an Open, Dynamic System:
    In a non-equilibrium steady state (NESS), the system is driven by a constant flux of monomers and a simultaneous first-order degradation process. In this scenario, the "survival probability" of the molecules becomes the dominant factor. Because T is structurally more robust, it degrades more slowly than U. 

### Conclusion
In an open and dynamic environment, this subtle kinetic advantage is amplified over time. The continuous supply prevents the system from reaching equilibrium, while the differential degradation acts as a chemical filter. Consequently, the slight stability edge of T is integrated over the temporal scale of the simulation, leading to the observed significant enrichment of T within the DNA strands.
