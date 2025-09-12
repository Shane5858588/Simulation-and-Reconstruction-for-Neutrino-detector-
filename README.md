# Cherenkov Track Simulation & Reconstruction (Poisson sampling)

**Notebook:** `sim_test1.ipynb`  
**Last updated:** 2025-09-12

This notebook simulates a straight charged‑particle track in a homogeneous medium, produces Cherenkov‑like light seen by a 2‑D array of DOMs, Poisson‑samples the light yield to create realistic hit counts, and then **reconstructs the track direction and impact parameters** from hit times and intensities. It concludes with a **resolution study** and several diagnostic plots.

---

## Contents
- [Overview](#overview)
- [Core Ideas](#core-ideas)
- [Dependencies](#dependencies)
- [Quickstart](#quickstart)
- [Key Parameters](#key-parameters)
- [Pipeline](#pipeline)
- [Outputs & Plots](#outputs--plots)
- [How to Read the Figures](#how-to-read-the-figures)
- [Functions Reference](#functions-reference)
- [Reproducing the Study](#reproducing-the-study)
- [Troubleshooting](#troubleshooting)
- [License](#license)

## Overview

- **Geometry:** A straight track parameterized by zenith/azimuth and an impact vector. DOMs are arranged on a 2‑D grid (configurable spacing).
- **Photon model:** Intensity follows an attenuated \(1/d^2\) law with exponential absorption `L_att` and angular factors from the Cherenkov cone.
- **Timing model:** Arrival times are computed from closest‑approach geometry and Cherenkov kinematics (uses `n_ref`, `sin_C`, `tan_C`).
- **Noise/realism:** Continuous intensities are **Poisson‑sampled** to discrete counts to emulate detection statistics.
- **Reconstruction:** Uses a **weighted PCA** for a physics‑informed seed, then a **χ²/likelihood minimization** (via `iminuit.Minuit`) to fit direction and impact.
- **Evaluation:** Generates distributions and correlations of angular error vs. observables (hit count, total signal, reduced χ²).

## Core Ideas

1. **Closest‑approach geometry:** `perp_distance(x0, u, x_dom)` gives the perpendicular miss distance between track and sensor.  
2. **Light yield:** `intensity(r_perp, L_att, sC, I0, d_floor)` encodes attenuation and geometric falloff.  
3. **Timing:** `arrival_time(s0, r_perp, c, n_ref, sin_C, tan_C)` maps geometry to first‑photon times.  
4. **Poissonization:** `generate_poisson_data(hits, E_total, K_min, random_seed)` converts continuous yields to integer counts.  
5. **Reconstruction:** `weighted_pca(...)` for seeding, then `chi2(...)` and `reconstruct_from_poisson_minuit(...)` refine parameters.  
6. **Resolution study:** `resolution_study(...)` loops over many synthetic events and aggregates metrics/plots.

## Dependencies

Install these Python packages (Python 3.11.13):  
- `iminuit`
- `matplotlib`
- `numpy`

```bash
pip install numpy matplotlib iminuit
# optional (if you run as a notebook):
pip install jupyter
```

> **Note:** `iminuit` provides a robust minimizer; if it's not installed you can switch to `scipy.optimize` with small code edits (see `reconstruct_from_poisson_minuit`).

## Quickstart

1. **Open the notebook**

   ```bash
   jupyter notebook sim_test1.ipynb
   ```

2. **Run all cells** (Kernel → Restart & Run All).  
   This will:
   - Build geometry and a truth track
   - Simulate continuous intensities and first‑photon times
   - Poisson‑sample hits
   - Reconstruct and produce resolution diagnostics

3. *(Optional)* **Export to a script**

   ```bash
   jupyter nbconvert --to script sim_test1.ipynb
   python sim_test1.py
   ```

## Key Parameters

Common knobs you’ll find near the top of the notebook:

- `n_ref = 1.33` — refractive index of the medium  
- `L_att = 100.0` — attenuation length (m)  
- Cherenkov geometry helpers: `theta_C`, `sin_C`, `tan_C`  
- Grid layout: `spacing` for inter‑DOM distance  
- Event brightness: `E_total`, per‑DOM cutoff `K_min`  
- Minimum distance clamp: `d_floor` to avoid singularities  
- Random seed: `random_seed` for reproducibility

## Pipeline

1. **Truth Setup:** choose `(theta, phi)` and a starting point `x0`; build unit vector `u`.  
2. **Project Geometry:** compute perpendicular distances and along‑track positions for each DOM.  
3. **Signal Model:** compute intensity and timing at each DOM; drop sensors below threshold.  
4. **Poisson Sampling:** draw counts for each hit DOM to mimic detected photoelectrons.  
5. **Initial Guess:** `weighted_pca` on space–time points to estimate direction/anchor.  
6. **Fit:** minimize `chi2` with `iminuit` to estimate angles and impact vector.  
7. **Metrics:** compare to truth → **angular error**; compute **reduced χ²**.  
8. **Plots:** generate diagnostic figures listed below.

## Outputs & Plots

This notebook produces the following figures:

- **Angular Error Distribution** — Histogram of the reconstruction error in degrees between the true track and the reconstructed track. Also shows a CDF curve for quick percentile reads (e.g., P50, P68).
- **Angular Error vs. Detected Hits** — Scatter plot showing how resolution improves (error drops) as the number of hit DOMs increases.
- **Angular Error vs. ΣK (total counts)** — Scatter plot of error against the summed Poisson counts across hit DOMs; points are color‑mapped by reduced χ² to visualize fit quality.
- **Reduced χ² Distribution & Error vs. Reduced χ²** — A histogram of reduced χ² to diagnose fit quality and a scatter of reduced χ² vs. angular error to expose outliers / failure modes.


The notebook displays figures inline. To save PNGs, add lines like:

```python
plt.tight_layout()
plt.savefig("fig_ang_err_hist.png", dpi=200)
```

## How to Read the Figures

**How to read the key figures**

- *Angular Error histogram*: A tight peak near 0° with a thin tail means the reconstructor is accurate and stable. Report median (P50), P68, P90 as headline numbers.
- *Error vs. Detected Hits*: Downward trend indicates more geometry/time constraints → better fits. A wide vertical spread at low hits signals under‑constrained events.
- *Error vs. ΣK*: Larger total light generally improves timing and weighting, thus lower error. Coloring by reduced χ² lets you see that high‑χ² points often drive the tail.
- *Reduced χ² plots*: A peak around 1 suggests the timing/variance model is well‑scaled. A heavy tail >~2 implies mis‑modeling (e.g., wrong σ, scattered light) or poor geometry.


## Functions Reference

Below is a lightweight index of functions discovered in the notebook (arguments only). For details, open the corresponding cell.

- `perp_distance(x0, u, x_dom)`
- `intensity(r_perp, L_att, sC, I0, d_floor)`
- `arrival_time(s0, r_perp, c, n, sinC, tanC)`
- `generate_poisson_data(hits, E_total, K_min, random_seed)`
- `generate_multiple_poisson_files(hits, n_files, E_total, K_min, start_seed, save_files, output_dir)`
- `angle_deg(u1, u2)`
- `reconstruct_from_poisson_minuit(DOM_positions, poisson_hits, sigma_t, d_floor)`
- `resolution_study(num_events, start_seed)`
- `angle_deg(u1, u2)`
- `u_from_angles(theta, phi)`
- `basis_from_u(u)`
- `model_I_t(theta, phi, b1, b2)`
- `residuals(theta, phi, b1, b2)`
- `chi2(theta, phi, b1, b2)`
- `chi2_mix(theta, phi, b1, b2, t_off, log_a)`
- `weighted_pca(points, w)`

## Reproducing the Study

Use the built‑in **`resolution_study(...)`** helper to sweep events:

```python
# Example
stats = resolution_study(num_events=1000, start_seed=12345)
print(stats["P50_deg"], stats["P68_deg"], stats["P90_deg"])
```

Typical headline metrics to report:
- **Median angular error (P50)**
- **P68 & P90 containment**
- Trend lines in *Error vs. Hits* and *Error vs. ΣK*
- Reduced χ² distribution (peak and tails)

## Troubleshooting

- **`ModuleNotFoundError: iminuit`** — `pip install iminuit` (or switch to SciPy’s `least_squares`/`minimize`).  
- **Unstable fits / big tails** — check `σ_t`, `d_floor`, and thresholds; ensure `weighted_pca` seeding is used.  
- **No points plotted** — thresholds (`K_min`) may be too high; reduce or increase `E_total`.  
- **Reduced χ² far from 1** — revisit timing variance assumptions or include scattering/jitter modeling.

## License

This repository is distributed for research and educational use. Include citation or acknowledgment if you use the code/ideas.

---

*Auto‑generated from `sim_test1.ipynb`.*
