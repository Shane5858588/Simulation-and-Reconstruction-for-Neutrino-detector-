# Muon Track Simulation Reconstruction

This repo contains a **clean, reproducible pipeline** to
1) **simulate** a straight muon-like track in water,  
2) **Poissonize** expected DOM light into discrete photoelectron counts, and  
3) **reconstruct** the track direction using a **time‑only χ²** fit,  
plus utilities to **batch-generate events**, **save/load NPZ files**, and **plot diagnostics**.

The code lives in a single notebook (`sim_test2.ipynb`) that is easy to read top‑to‑bottom.

---

## What’s inside (high level)

- **Simulation block**: detector geometry (3D DOM grid), muon track (Free to change the starting position and angle), analytic Cherenkov timing & intensity, hit selection.
- **Poissonization block**: converts expected intensities → integer counts with thresholds; optionally saves each event to disk.
- **Reconstruction block**: time‑only χ² fit of `(θ, φ, b1, b2, t_off)`; dual‑seed (`u` and `−u`) to curb 180° flips.
- **Plotting block**: Angular Error CDF, Reduced χ² histogram, Angular error vs hits, Zenith/Azimuth, ΣK (energy proxy).
- **Executive Summary**: boss‑friendly KPIs (P50/P68/P90, χ²/ndof stats, flip rate, hits/ΣK means).

---

## Quickstart

> Open and run `sim_test2.ipynb` **top to bottom**. These are the key cells you’ll use most often:

### 1) Generate & Save Poisson Events
```python
OUT_DIR = "/path/to/poisson_runs/runA"
poisson_paths, stats = generate_multiple_poisson_files(
    hits,              # produced by the simulation block
    n_files=200,
    E_total=E_TOTAL,   # total expected counts per event
    K_min=K_MIN,       # DOM detection threshold in counts
    start_seed=0,
    output_dir=OUT_DIR
)
```
| Label | Start `x0` (m)   | Zenith (°)        | Azimuth (°)            | Why this case matters                                      |
| ----: | ---------------- | ----------------- | ---------------------- | ---------------------------------------------------------- |
|     A | (   0,   0,   0) | 120               | 45                     | Baseline (your default truth).                             |
|     B | ( 300,   0,   0) | 120               | 45                     | Shifted inside volume; similar direction.                  |
|     C | (   0, 300,   0) | 60                | 135                    | Up-going-ish, diagonal azimuth.                            |
|     D | (   0,   0, 300) | 30                | 270                    | Shallow zenith (near down-going), axis-aligned az.         |
|     E | (-600,   0,   0) | 100               | 20                     | Starts just outside x-min; through-going.                  |
|     F | (   0,-600,   0) | 150               | 200                    | Starts just outside y-min; near horizontal.                |
|     G | (   0,   0,-600) | 90                | 0                      | Exactly horizontal along +x (max Cherenkov cone symmetry). |
|     H | ( 500, 500, 500) | 80                | 225                    | Corner entry, slanted track.                               |
|     I | (-500, 500,-500) | 140               | 315                    | Opposite corner; near-horizontal.                          |
|     J | (   0,   0,   0) | 60                | 225                    | Mirror of C around vertical; checks azimuth wrapping.      |
|     K | (   0,   0,   0) | **60**            | **225+180** (=405→225) | **Antipode of J** → probes 180° flip handling.             |
|     L | (   0,   0,   0) | **180−120** (=60) | **45+180** (=225)      | **Antipode of A** → same line, reverse direction.          |

### 2) Load Events
```python
import os, glob, numpy as np

poisson_paths = sorted(glob.glob(os.path.join(OUT_DIR, "poisson_evt_*.npz")))
events = []
for pth in poisson_paths:
    with np.load(pth, allow_pickle=True) as f:
        ph = f["poisson_hits"]
    # list[(idx, t_ns, K)]
    events.append([(int(i), float(t), int(k)) for (i,t,k) in np.array(ph, dtype=object).tolist()])
```

### 3) Reconstruct
```python
results = []
for ph in events:
    r = reconstruct_one(ph)   # time-only χ² fit
    if r.get("success"):
        results.append(r)

print(f"Reconstructed {len(results)} / {len(events)} events")
```

### 4) Plots that matter
- **Angular Error CDF** → quote P50 / P68 / P90 (resolution & containment)
- **Reduced χ² (time-only)** → should peak near 1 if model & σₜ are consistent
- **Angular Error vs Hits** → resolution improves ~ 1/√N
- **Zenith / Azimuth** → should peak at truth angles
- **ΣK** → brightness proxy; brighter events reconstruct better

### 5) Executive Summary (KPIs)
The notebook includes a one‑cell summary that prints:
- P50 / P68 / P90 angular error  
- χ²/ndof mean & median (time‑only)  
- Flip rate (±u ambiguity)  
- Mean detected hits and ΣK

---

## Installation

Requires **Python 3.9+**.

```bash
pip install numpy matplotlib scipy
# optional, if using iminuit flavor:
pip install iminuit
```

---

## Configuration knobs (most useful)

| Group         | Name            | Meaning |
|---------------|-----------------|---------|
| Geometry      | `n_dom_x/y/z`   | DOM grid dimensions |
|               | `spacing`       | Inter‑DOM spacing (m) |
| Physics       | `n_ref`         | Refractive index (water) |
|               | `L_att`         | Attenuation length (m) |
|               | `theta_C`       | Cherenkov angle (derived from `n_ref`) |
| Simulation    | `I_0`           | Intensity scale factor |
|               | `I_thresh`      | Intensity cut for accepting a DOM hit |
|               | `σ_t` (in code: `sigma_t`) | Timing jitter (ns) |
| Poisson       | `E_TOTAL`       | Target expected total counts per event |
|               | `K_MIN`         | Per‑DOM count threshold |
| Reconstruction| `dual-seed`     | Fit both `u` and `−u` to reduce flips |
|               | `t_off`         | Global timing offset parameter |

---

## Functions (what they do)

### Simulation
- **`perp_distance(x0, u, x_dom)`** → Perpendicular distance $r_\perp$ from DOM to the infinite track.  
- **`intensity(r_perp)`** → Expected relative intensity $\propto e^{-d/L_{\text{att}}} / d^2$, with $d=r_\perp/\sin\theta_C$ and a near-field clamp.  
- **`arrival_time(s0, r_perp)`** → First-photon Cherenkov arrival time at a DOM (ns):

$$
t_{\text{hit}}
=\frac{s_0}{c}
-\frac{r_\perp}{c\,\tan\theta_C}
+\frac{n}{c}\,\frac{r_\perp}{\sin\theta_C}
+t_{\text{off}}
$$

### Poissonization
- **`generate_poisson_data(hits, E_total, K_min, seed)`** → Converts continuous expectations to integer counts: $K \sim \text{Poisson}(\lambda)$, scaled so $\sum \lambda = E_{\text{total}}$. Drops DOMs with $K < K_{\min}$.

- **`generate_multiple_poisson_files(hits, n_files, E_total, K_min, start_seed, output_dir)`** → Repeats `generate_poisson_data` and writes each event as `poisson_evt_XXXX.npz` under `output_dir`.


### Reconstruction
- **`reconstruct_one(poisson_hits)`** → Time‑only χ² fit for `(θ, φ, b1, b2, t_off)`. Uses dual‑seed (`u` and `−u`) to mitigate 180° flips. Returns a dict including:
  - `u_hat`, `theta_hat_deg`, `phi_hat_deg`  
  - `chi2_time`, `ndof_time` (**use these** for reduced χ²)  
  - `detected_hits`, `sumK`

### Loading
- NPZ reader (simple loop) → Loads `poisson_evt_*.npz` files and converts `poisson_hits` back into `list[(idx, t_ns, K)]`.

### Plotting (notebook cells)
- Angular Error CDF; Reduced χ² histogram; Angular error vs hits; Zenith/Azimuth histograms; ΣK distribution.

### Summary (notebook cell)
- P50/P68/P90 angular error; χ²/ndof mean/median; flip rate; mean hits & ΣK.

---

## Data formats

### Input/Intermediate (in‑memory)
- **`hits`**: `list[(dom_idx, t_ns, intensity)]` (continuous simulation)  
- **`poisson_hits`**: `list[(dom_idx, t_ns, K)]` (discrete counts)  

### On disk
- **`poisson_evt_XXXX.npz`** contains at least:
  - `poisson_hits` (object array convertible to `list[(idx, t, K)]`)
  - optional metadata: `file_id`, `n_detected_hits`, `total_counts`, `random_seed`

---

## Interpreting the plots

- **Angular Error CDF** → Resolution & containment (quote P50/P68/P90).  
- **Reduced χ² (time‑only)** → Fit quality; should center near 1. Tails/bimodality can indicate model mismatch or 180° flips.  
- **Angular Error vs Hits** → Expect ~ 1/√N trend (more hits → smaller error).  
- **Zenith/Azimuth** → Should peak at truth angles for closure tests.  
- **ΣK** → Brighter events tend to reconstruct better.

---

## Troubleshooting & tips

- **Use the *time‑only* χ²** and its `ndof = N_hits − 5` when you plot reduced χ². Don’t mix with Poisson C‑stat.  
- **Angles wrap**: compare to truth using `zenith % 180` and `azimuth % 360`.  
- **Flip ambiguity**: horizontal or low‑hit events flip more often; keep the dual‑seed and forward‑time check enabled.  
- **Consistent paths**: save & load using `os.path.join(OUT_DIR, "poisson_evt_*.npz")`.

---

## Environment

```bash
python --version    # 3.9+ recommended
pip install numpy matplotlib scipy
# optional (alternative fitter)
pip install iminuit
```

---

## Citation / Acknowledgment

If this helps your work, a citation or acknowledgment is appreciated.  
This pipeline was designed to benchmark **time‑only χ²** reconstruction, study **resolution vs statistics**, and validate **Poissonized timing** in water/ice Cherenkov‑style arrays.
