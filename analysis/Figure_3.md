<img width="7719" height="4424" alt="sensitivity_v601" src="https://github.com/user-attachments/assets/c70d592c-6d9d-4af2-8ca2-a93f600c6b62" />

## 🔬 Figure 3: Comprehensive 2D Parameter Sweep (Temperature vs. pH Phase Space)

### 1. Technical Overview ("The What")
This 12-panel grid maps the complete thermodynamic and chemical phase space of the UPDSF v6.0.1 framework across a dual-parameter sweep of Temperature ($40\text{--}90^\circ\text{C}$) and pH ($5\text{--}11$). It visualizes global trends in total Chemical Selection Index (CSI), direct structural preferences, pool-composition drift, control controls (zero UV, equal pools), UV destruction fractions, deamination kinetics, and numerical error bounds.

### 2. Scientific Insight ("The So What")
- *The "Hot-Alkaline" Optimal Corner: Across nearly all selection and enrichment metrics, the top-right corner of the phase space—representing maximum temperature ($90^\circ\text{C}$) and maximum alkalinity ($\text{pH } 11$)—stands out as the absolute global optimum for Thymine enrichment. 
- Systematic Thermodynamic Response: The smooth gradients across the heatmaps prove that T-enrichment is not a stochastic artifact, but a systematic, predictable thermodynamic response to early-Earth environmental stressors.
- Isolating Mechanisms: The control panels (`CSI - zero UV control` and `CSI - equal pools`) decouple the drivers of selection, confirming that both environmental degradation kinetics and direct structural biases cooperate synergistically at high pH and temperature.

### 3. Key Metric Highlights
- Peak CSI Hotspot: Reaches its maximum values ($>2.60\times$) exclusively at the extreme boundary of $T = 90^\circ\text{C}$ and $\text{pH } 11$.
- Residue Retention Phase Space: Shows that structural integrity remains robust (>80% retention) across broad regions of the grid before dropping off sharply near extreme thermal/pH limits.
- Mass Conservation Uniformity: The error-bound panel (`log10 max |consumed identity rel err|`) maintains a uniform baseline around $10^{-12}$ to $10^{-15}$ across the entire 2D matrix, verifying numerical stability in every simulated condition.

### 4. Prebiotic Implication
This comprehensive parameter sweep provides robust justification for alkaline hydrothermal vent scenarios as cradles of prebiotic chemistry. By mapping the entire phase space, the model demonstrates that the conditions favoring Thymine over Uracil are strictly bounded to high-temperature, alkaline niches—neatly matching the geochemical realities of Hadean submarine hydrothermal systems.
