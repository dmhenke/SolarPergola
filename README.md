# Solar Pergola Structural Sizing & Regulatory Compliance Engine

> **Live Web Report & Interactive Sizing Engine:** [https://dmhenke.github.io/SolarPergola/](https://dmhenke.github.io/SolarPergola/)

Preliminary structural sizing and dimensional limit analysis for freestanding 4-post steel solar pergolas in **Salt Lake City, UT (ZIP 84108)** and surrounding jurisdictions under:
* **2024 International Building Code (IBC)**
* **ASCE 7-22** (Minimum Design Loads and Associated Criteria for Buildings and Other Structures)
* **AISC 360-22** (Specification for Structural Steel Buildings)
* **Aluminum Design Manual (ADM 2020)**
* **Utah Code Title 15A Chapter 3 §107** (State Building Code Amendments for Snow & Seismic Loads)
* **Salt Lake City Code Title 21A.40.050** (Accessory Structure Zoning Regulations)

---

## 1. Dimensional Coordinate System

To ensure seamless translation to 3D CAD/BIM modeling and parametric generators:

```
              High Side Beam (Girder)
           ●===========================●
          /|                          /|
         / |                         / |
        /  |                        /  |
       /   |                       /   |   ▲
      /    |                      /    |   │
     /     |                     /     |   │ Slope / Rafter Span (W)
    /      |                    /      |   │ [12 ft, 16 ft, 20 ft]
   /       ●===================/=======●   ▼
  /       /                   /       /
 /       /                   /       /
●===========================●       /
|       /                   |      /
|      /                    |     /
|     /                     |    /
|    /                      |   /
|   /                       |  /
●==/========================● /
| /                         |/
●                           ●
◄───────────────────────────►
  Post Spacing / Girder Span (L)
     [16 ft, 20 ft, 24 ft]
```

* **$W$ (Depth)**: Dimension along the roof pitch / slope. Spanned by rafters.
* **$L$ (Width)**: Dimension along the horizontal perimeter. Spanned by girders / perimeter beams between uprights (post spacing).
* **$\theta$ (Pitch)**: Roof angle ($7^\circ, 10^\circ, 15^\circ, 20^\circ, 25^\circ, 30^\circ, 35^\circ$).
* **Clear Height**: $8.0\text{ ft}$ standard clear clearance under low-edge steel framing.

---

## 2. Architecture & Computational Pipeline

All structural mechanics, code standards, and section databases are centralized in `solar_pergola_engine.R`. Scripts run independently without circular dependencies:

```
                      ┌─────────────────────────┐
                      │ solar_pergola_engine.R  │  (Central Core)
                      └────────────┬────────────┘
         ┌─────────────────────────┼─────────────────────────┐
         │                         │                         │
         ▼                         ▼                         ▼
┌──────────────────────┐  ┌──────────────────┐  ┌──────────────────┐
│solar_pergola_sizing.R│  │max_post_spacing.R│  │  hss6x2_beam.R   │
└──────────┬───────────┘  └────────┬─────────┘  └────────┬─────────┘
           │                       │                     │
           ▼                       ▼                     ▼
 frame_sizing.csv         max_post_spacing.csv   beam_max_span.csv
 rail_module_checks.csv                          upright_max_spacing.csv
                                                 hss6x2_properties.csv
                                   │
                                   ▼
                           ┌───────────────┐
                           │ max_span_v2.R │
                           └───────┬───────┘
                                   │
                                   ▼
                            max_span_v2.csv
```

---

## 3. Key Engineering & Regulatory Features

### A. Elevation-Harmonized Snow Loads (Utah Code §15A-3-107)
Ground snow load ($P_g$) is computed via `get_snow_loads(elevation_ft)` using the Utah Snow Load Study (USU) model calibrated for Salt Lake County:
* **Valley Floor** ($4,239\text{--}4,500\text{ ft}$): $P_g = 28\text{--}35\text{ psf}$ ($P_f \approx 24\text{--}29\text{ psf}$).
* **East Bench / Foothill** ($5,500\text{--}5,884\text{ ft}$ in ZIP 84108): $P_g = 90.1\text{--}95\text{ psf}$ ($P_f \approx 75.7\text{--}79.8\text{ psf}$).
* **Seismic-Snow Co-Load**: Implements Utah Code §15A-3-107: $S_{seis} = [0.15 + 0.016(A - 5)] P_{roof}$ for $P_f > 30\text{ psf}$ (where $A$ is elevation in thousands of feet), replacing the generic ASCE 7-22 $0.20$ flat factor.

### B. Unbraced Cantilever vs. Knee-Braced Framing
The engine models two distinct structural behaviors:
1. **Unbraced Cantilever Column (Flagpole Mode)**:
   * Ordinary Cantilever Column ($R = 1.25, K = 2.1$).
   * Full base moment transfer into foundation.
   * Result: S-shapes fail slenderness ($KL/r_y > 300$). Pier depths require **$8.6\text{ to }14.3\text{ ft}$**.
2. **Knee-Braced Frame Mode (Brace Arm $a_b = 2.5\text{ ft}$)**:
   * Intermediate lateral support drops effective unbraced column length to $KL = 1.0 \times (H_{col} - a_b)$.
   * Reduces column base moment to zero (pinned/shear-only base).
   * Result: **S-shapes (S4x7.7 to S8x18.4) fully qualify as posts**. Drilled pier depths drop from $10\text{--}14\text{ ft}$ down to **$3.4\text{ to }5.5\text{ ft}$** (saving substantial drilling and concrete costs).

### C. Regulatory & Equipment Sanity Gating
All generated tables evaluate three automatic go/no-go gates:
1. **Zoning Max Height Gate** (`zoning_ok`): Checks if total structure apex height $\le 15.0\text{ ft}$ (Salt Lake City §21A.40.050 accessory structure limit).
2. **PV Module Wind Uplift Gate** (`mod_rear_ok`): Compares ASCE 7-22 C&C negative wind pressure against module certified rear design rating ($33.4\text{ psf}$ / $2400\text{ Pa}$ test load with $1.5$ safety factor).
3. **Column Slenderness Gate**: Rejects any column where $KL/r > 300$ (AISC 360-22 Section E2).

---

## 4. Summary of Output Data Files

| File | Description | Key Columns |
| :--- | :--- | :--- |
| `frame_sizing.csv` | Full pergola system sizing across sizes (S, M, L) and pitches ($7^\circ\text{--}35^\circ$). | `size`, `pitch`, `top_ht_ft`, `zoning_ok`, `mod_rear_ok`, `status`, `rafter_W`, `girder_W`, `col_cant_HSS`, `pier_cant_depth_ft`, `col_braced_HSS`, `col_braced_S`, `pier_braced_depth_ft`. |
| `rail_module_checks.csv` | Solar rail spans (IronRidge XR10, XR100, XR1000) and module front/rear pressure checks. | `size`, `pitch`, `XR10_in`, `XR100_in`, `XR1000_in`, `mod_front_psf`, `mod_rear_psf`, `mod_rear_ok`, `status`. |
| `max_post_spacing.csv` | Maximum allowable post-to-post clear spacing ($L$) for S-shapes and Square HSS. | `W_ft`, `pitch`, `shape`, `girder_maxL`, `post_cant_maxL`, `post_braced_maxL`, `zoning_ok`. |
| `beam_max_span.csv` | Master lookup of maximum girder spans across all shapes including rectangular HSS 6x2. | `W_ft`, `pitch`, `shape`, `beam_maxL`, `top_ht_ft`, `zoning_ok`. |
| `upright_max_spacing.csv` | Maximum post spacing governed by upright capacity (cantilever vs. knee-braced). | `W_ft`, `pitch`, `shape`, `post_cant_maxL`, `post_braced_maxL`, `zoning_ok`. |
| `hss6x2_properties.csv` | Exact integrated cross-sectional properties for HSS 6x2 rectangular tubes. | `shape`, `t`, `wt`, `A`, `Ix`, `Sx`, `Zx`, `ry`, `J`. |
| `max_span_v2.csv` | Parametric 4-post clear span limits, upright selection, total steel weight, and estimated material cost. | `case`, `pitch`, `beam`, `max_span_ft`, `upright_cant`, `upright_braced`, `top_ht_ft`, `zoning_ok`, `steel_lb`, `cost_mid`. |

---

## 5. How to Run the Scripts

Ensure `R` is installed. Run commands from the project root:

```bash
# 1. Run main sizing matrix for standard S, M, L pergolas
Rscript solar_pergola_sizing.R

# 2. Run post spacing limits for S-shapes and Square HSS
Rscript max_post_spacing.R

# 3. Run HSS 6x2 perimeter beam and combine master span tables
Rscript hss6x2_beam.R

# 4. Run parametric clear span study (default elevation 5500 ft)
Rscript max_span_v2.R

# (Optional) Run parametric clear span study for valley floor (4300 ft elevation)
Rscript max_span_v2.R 4300
```
