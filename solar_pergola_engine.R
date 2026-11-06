## =====================================================================
## Solar Pergola Structural & Regulatory Engine (Core Module)
## Codes: 2024 IBC / ASCE 7-22, AISC 360-22, ADM 2020, Utah Code §15A-3-107
##
## This engine provides:
##  1. AISC 360-22 shape databases (W, Square HSS, Rect HSS 6x2, S-shapes)
##     with exact section properties (including true web thickness tw).
##  2. Flexural, axial, shear, and LTB capacities (LRFD & ASD).
##  3. ASCE 7-22 wind and C&C / MWFRS monoslope open roof pressure coefficients.
##  4. Elevation-harmonized snow loading per Utah Code §15A-3-107.
##  5. Dual-mode column & girder analysis: Unbraced Cantilever vs. Knee-Braced.
##  6. Drilled pier foundation depth sizing (2024 IBC §1807.3.2.1 Formula 18-1).
##  7. Regulatory & equipment sanity gating (Zoning height, Module UL rating).
## =====================================================================

options(stringsAsFactors = FALSE)

## ---------------------------------------------------------------------
## 1. MATERIAL & SITE DEFAULTS
## ---------------------------------------------------------------------
Es <- 29000  # Steel Young's modulus, ksi
E_al <- 10100 # Aluminum Young's modulus, ksi
Fy_al <- 34   # Aluminum yield strength, ksi (6005-T5/T61)

default_site <- list(
  elevation_ft = 5500,  # Salt Lake City East Bench (ZIP 84108)
  Ce           = 1.0,   # Exposure C/B, partially exposed
  Ct           = 1.2,   # Unheated open-air structure (ASCE 7-22 Table 7.3-2)
  V            = 115,   # mph, basic wind speed (SLC residential criterion)
  Kzt          = 1.0,   # Topographic factor (verify on hilltops)
  Ke           = 1.0,   # Ground elevation factor (conservative)
  Kd           = 0.85,  # Wind directionality factor
  G            = 0.85,  # Gust effect factor
  SDS          = 1.10,  # Design spectral acceleration (ZIP 84108)
  R_cant       = 1.25,  # Steel ordinary cantilever column system
  Omega0_cant  = 1.25,
  Cd_cant      = 1.25,
  R_braced     = 3.25,  # Ordinary steel braced / moment frame system
  Omega0_braced= 2.0,
  Cd_braced    = 3.25,
  rho          = 1.3,   # Redundancy factor
  frost_in     = 30,    # Frost depth, inches
  S_lat        = 100,   # Allowable lateral soil bearing, psf/ft (IBC Table 1806.2)
  q_bear       = 1500,  # Allowable vertical soil bearing, psf
  clear_ft     = 8.0,   # Clear height under framing low edge, ft
  zoning_max_ht= 15.0   # Max height for detached accessory structure (SLC §21A.40.050)
)

default_pv <- list(
  L_in              = 93.8,
  W_in              = 44.6,
  wt_lb             = 74.1,
  front_design_psf  = 5400 / 1.5 * 0.020885,  # 75.2 psf (5400 Pa test / 1.5)
  rear_design_psf   = 2400 / 1.5 * 0.020885   # 33.4 psf (2400 Pa test / 1.5)
)
default_pv$area_ft2    <- default_pv$L_in * default_pv$W_in / 144
default_pv$D_psf       <- default_pv$wt_lb / default_pv$area_ft2 + 0.25 # + hardware
default_pv$rail_trib_ft<- default_pv$L_in / 2 / 12  # portrait, 2 rails/module

## IronRidge rail properties (IronRidge RM Design Guide)
rails <- data.frame(
  rail = c("XR10", "XR100", "XR1000"),
  Ix   = c(0.124, 0.390, 0.843),
  Sx   = c(0.136, 0.297, 0.530),
  wt   = c(0.436, 0.680, 0.945)
)
rails$Ma_lbft <- rails$Sx * Fy_al / 1.65 * 1000 / 12 # ASD allowable moment

## ---------------------------------------------------------------------
## 2. STRUCTURAL STEEL SHAPES DATABASE
## ---------------------------------------------------------------------
## W-shapes (AISC Manual v16)
W_shapes <- read.table(header = TRUE, text = "
shape   wt   A    d     Ix   Sx   Zx   Iy    ry    rts   J      ho
S4x7.7  7.7  2.26 4.00  6.08 3.04 3.51 0.764 0.581 0.725 0.0703 3.71
W4x13   13   3.83 4.16  11.3 5.46 6.28 3.86  1.00  1.16  0.151  3.81
W6x9    9    2.68 5.90  16.4 5.56 6.23 2.20  0.905 1.06  0.0405 5.69
W6x12   12   3.55 6.03  22.1 7.31 8.30 2.99  0.918 1.08  0.0903 5.75
W8x10   10   2.96 7.89  30.8 7.81 8.87 2.09  0.841 1.01  0.0426 7.69
W8x13   13   3.84 7.99  39.6 9.91 11.4 2.73  0.843 1.03  0.0870 7.74
W8x18   18   5.26 8.14  61.9 15.2 17.0 7.97  1.23  1.43  0.172  7.81
W10x12  12   3.54 9.87  53.8 10.9 12.6 2.18  0.785 0.983 0.0547 9.66
W10x15  15   4.41 9.99  68.9 13.8 16.0 2.89  0.810 1.01  0.104  9.72
W10x19  19   5.62 10.2  96.3 18.8 21.6 4.29  0.874 1.06  0.233  9.85
W10x22  22   6.49 10.2  118  23.2 26.0 11.4  1.33  1.55  0.239  9.86
W12x16  16   4.71 12.0  103  17.1 20.1 2.82  0.773 0.994 0.103  11.7
W12x19  19   5.57 12.2  130  21.3 24.7 3.76  0.822 1.03  0.180  11.9
W12x22  22   6.48 12.3  156  25.4 29.3 4.66  0.848 1.04  0.293  11.9
W12x26  26   7.65 12.2  204  33.4 37.2 17.3  1.51  1.75  0.300  11.8
W14x22  22   6.49 13.7  199  29.0 33.2 7.00  1.04  1.27  0.208  13.4
W14x26  26   7.69 13.9  245  35.3 40.2 8.91  1.08  1.31  0.358  13.5
W16x26  26   7.68 15.7  301  38.4 44.2 9.59  1.12  1.36  0.262  15.4
W16x31  31   9.13 15.9  375  47.2 54.0 12.4  1.17  1.39  0.461  15.5
W18x35  35   10.3 17.7  510  57.6 66.5 15.3  1.22  1.51  0.506  17.3
")
W_shapes$Fy <- ifelse(substr(W_shapes$shape, 1, 1) == "S", 36, 50)

## S-shapes with EXACT AISC web thicknesses (tw) and properties
S_shapes <- read.table(header = TRUE, text = "
shape    wt    A    d    bf   tw    tf    Ix   Sx   Zx   Iy    Sy    Zy    ry    rts   J      ho
S4x7.7   7.7   2.26 4.00 2.66 0.193 0.293 6.08 3.04 3.51 0.764 0.574 0.966 0.581 0.719 0.0703 3.71
S4x9.5   9.5   2.79 4.00 2.80 0.326 0.293 6.79 3.39 4.01 0.903 0.646 1.12  0.569 0.717 0.119  3.71
S5x10    10.0  2.93 5.00 3.00 0.214 0.326 12.3 4.90 5.66 1.19  0.795 1.36  0.638 0.795 0.114  4.67
S6x12.5  12.5  3.66 6.00 3.33 0.232 0.359 22.0 7.33 8.45 1.80  1.08  1.86  0.702 0.881 0.165  5.64
S6x17.25 17.25 5.06 6.00 3.57 0.465 0.359 26.2 8.74 10.5 2.31  1.30  2.36  0.676 0.870 0.369  5.64
S8x18.4  18.4  5.40 8.00 4.00 0.271 0.425 57.6 14.4 16.5 3.73  1.87  3.27  0.831 1.04  0.335  7.57
")
S_shapes$Fy <- 36

## Square HSS shapes (AISC Manual v16, A500 Gr B/C, Fy = 46 ksi)
HSS_sq <- read.table(header = TRUE, text = "
shape        b  tdes  wt    A     I     S     Z     r
HSS4x4x3/16  4  0.174 9.42  2.58  6.21  3.10  3.67  1.55
HSS4x4x1/4   4  0.233 12.21 3.37  7.80  3.90  4.69  1.52
HSS4x4x5/16  4  0.291 14.83 4.10  9.14  4.57  5.62  1.49
HSS4x4x1/2   4  0.465 21.63 6.02  11.9  5.97  7.70  1.41
HSS5x5x1/4   5  0.233 15.62 4.30  16.0  6.41  7.61  1.93
HSS5x5x5/16  5  0.291 19.08 5.26  19.0  7.62  9.20  1.90
HSS5x5x1/2   5  0.465 28.43 7.88  26.0  10.4  13.1  1.82
HSS6x6x1/4   6  0.233 19.02 5.24  28.6  9.54  11.2  2.34
HSS6x6x5/16  6  0.291 23.34 6.43  34.3  11.4  13.6  2.31
HSS6x6x3/8   6  0.349 27.48 7.58  39.5  13.2  15.8  2.28
HSS6x6x1/2   6  0.465 35.24 9.74  48.3  16.1  20.0  2.23
HSS6x6x5/8   6  0.581 42.30 11.7  55.2  18.4  23.2  2.17
HSS8x8x1/4   8  0.233 25.82 7.10  70.7  17.7  20.5  3.15
HSS8x8x5/16  8  0.291 31.84 8.76  85.6  21.4  25.1  3.13
HSS8x8x3/8   8  0.349 37.69 10.4  100   25.0  29.4  3.10
HSS8x8x1/2   8  0.465 48.85 13.5  125   31.2  37.5  3.04
HSS8x8x5/8   8  0.581 59.32 16.4  146   36.5  44.7  2.99
HSS10x10x5/16 10 0.291 40.35 11.1 172   34.3  39.6  3.94
HSS10x10x3/8 10 0.349 48.24 13.2  202   40.4  46.7  3.91
HSS10x10x1/2 10 0.465 62.46 17.2  256   51.2  60.7  3.86
HSS12x12x3/8 12 0.349 58.45 16.0  357   59.6  68.9  4.72
HSS12x12x1/2 12 0.465 76.07 21.0  459   76.5  89.6  4.68
")
HSS_sq$Fy <- 46

## Custom Rectangular HSS Section Properties Integrator
rect_props <- function(B, Hh, t) {
  n <- 1200; ys <- seq(-Hh / 2, Hh / 2, length.out = n); dy <- ys[2] - ys[1]
  ro <- 2 * t; ri <- t
  width_at <- function(y, b, h, r) {
    ay <- abs(y); if (ay > h / 2) return(0)
    if (ay <= h / 2 - r) return(b)
    dy <- ay - (h / 2 - r); b - 2 * (r - sqrt(max(r^2 - dy^2, 0)))
  }
  wo <- sapply(ys, width_at, b = B, h = Hh, r = ro)
  wi <- sapply(ys, width_at, b = B - 2 * t, h = Hh - 2 * t, r = ri)
  w <- wo - wi
  A <- sum(w) * dy; Ix <- sum(w * ys^2) * dy; Zx <- 2 * sum(w[ys > 0] * ys[ys > 0]) * dy
  xs <- seq(-B / 2, B / 2, length.out = n); dx <- xs[2] - xs[1]
  ho <- sapply(xs, width_at, b = Hh, h = B, r = ro)
  hi <- sapply(xs, width_at, b = Hh - 2 * t, h = B - 2 * t, r = ri)
  Iy <- sum((ho - hi) * xs^2) * dx
  J  <- 2 * t^2 * (B - t)^2 * (Hh - t)^2 / (B * t + Hh * t - 2 * t^2)
  list(A = A, Ix = Ix, Sx = Ix / (Hh / 2), Zx = Zx, Iy = Iy, ry = sqrt(Iy / A),
       J = J, wt = A * 3.40)
}

walls_R6x2 <- data.frame(nom = c("1/8", "3/16", "1/4", "5/16", "3/8"),
                         t = c(0.116, 0.174, 0.233, 0.291, 0.349))
HSS_6x2 <- do.call(rbind, lapply(seq_len(nrow(walls_R6x2)), function(i) {
  p <- rect_props(2, 6, walls_R6x2$t[i])
  data.frame(shape = paste0("HSS6x2x", walls_R6x2$nom[i]), t = walls_R6x2$t[i],
             wt = round(p$wt, 2), A = round(p$A, 2), Ix = round(p$Ix, 2),
             Sx = round(p$Sx, 2), Zx = round(p$Zx, 2), ry = round(p$ry, 3),
             J = round(p$J, 2))
}))
HSS_6x2$Fy <- 46

## ---------------------------------------------------------------------
## 3. HARMONIZED ELEVATION-BASED SNOW LOAD (Utah Code §15A-3-107)
## ---------------------------------------------------------------------
get_snow_loads <- function(elevation_ft = 5500, Ce = 1.0, Ct = 1.2, override_pg = NULL) {
  base_elev <- 4239  # ft MSL, Salt Lake City base datum
  base_pg   <- 28.0  # psf, Salt Lake City code minimum
  
  if (!is.null(override_pg)) {
    pg <- override_pg
  } else if (elevation_ft <= base_elev) {
    pg <- base_pg
  } else {
    # Utah Snow Study (utahsnowload.usu.edu) interpolation for Wasatch Front
    known_elev <- c(4239, 4500, 5000, 5515, 5884)
    known_pg   <- c(28.0, 35.0, 60.0, 91.0, 95.0)
    pg <- approx(known_elev, known_pg, xout = elevation_ft, rule = 2)$y
  }
  
  # ASCE 7-22 §7.3 Flat roof snow pf (strength level)
  pf <- 0.7 * Ce * Ct * pg
  # ASD flat roof snow (ASCE 7-22 load combination factor 0.7 on S)
  pf_asd <- 0.7 * pf
  
  # Utah Code §15A-3-107 Seismic snow co-load factor
  elev_kft <- elevation_ft / 1000
  seismic_coeff <- if (pf > 30) min(max(0.15 + 0.016 * (elev_kft - 5), 0.15), 0.25) else 0.0
  
  list(elevation_ft = elevation_ft, pg = round(pg, 1), pf = round(pf, 1),
       pf_asd = round(pf_asd, 1), seismic_coeff = round(seismic_coeff, 3))
}

## Roof slope reduction Cs (slippery PV surface, Ct = 1.2)
Cs_slip <- function(th) ifelse(th <= 15, 1.0, pmax(0, 1.0 - (th - 15) / 55))

## ---------------------------------------------------------------------
## 4. ASCE 7-22 WIND LOADS (Monoslope Free Roof)
## ---------------------------------------------------------------------
Kz <- function(z) 2.41 * (max(z, 15) / 2460)^(2 / 9.8) # Exposure C
qKd <- function(h, site = default_site) {
  0.00256 * Kz(h) * site$Kzt * site$Ke * site$V^2 * site$Kd
}

## C&C Pressure Coefficients (Figure 30.7-1, Monoslope Free Roof, Clear Flow)
cc_tab <- data.frame(
  ang = c(0, 7.5, 15, 30, 45),
  b1p = c(2.4, 3.2, 3.6, 5.2, 5.2), b1n = c(-3.3, -4.2, -3.8, -5.0, -4.6),
  b2p = c(1.8, 2.4, 2.7, 3.9, 3.9), b2n = c(-1.7, -2.1, -2.9, -3.8, -3.5),
  b3p = c(1.2, 1.6, 1.8, 2.6, 2.6), b3n = c(-1.1, -1.4, -1.9, -2.5, -2.3)
)
CN_cc <- function(th, EA, a) {
  band <- if (EA <= a^2) 1 else if (EA <= 4 * a^2) 2 else 3
  p <- approx(cc_tab$ang, cc_tab[[paste0("b", band, "p")]], th)$y
  n <- approx(cc_tab$ang, cc_tab[[paste0("b", band, "n")]], th)$y
  c(pos = p, neg = n)
}

## MWFRS Pressure Coefficients (Figure 27.3-4, Clear Flow)
mw_tab <- data.frame(
  ang = c(0, 7.5, 15, 22.5, 30, 37.5, 45),
  A0W = c(1.2, -0.6, -0.9, -1.5, -1.8, -1.8, -1.6), A0L = c(0.3, -1.0, -1.3, -1.6, -1.8, -1.8, -1.8),
  B0W = c(-1.1, -1.4, -1.9, -2.4, -2.5, -2.4, -2.3), B0L = c(-0.1, 0.0, 0.0, -0.3, -0.5, -0.6, -0.7),
  A1W = c(1.2, 0.9, 1.3, 1.7, 2.1, 2.1, 2.2),       A1L = c(0.3, 1.5, 1.6, 1.8, 2.1, 2.2, 2.5),
  B1W = c(-1.1, 1.6, 1.8, 2.2, 2.6, 2.7, 2.6),      B1L = c(-0.1, 0.3, 0.6, 0.7, 1.0, 1.1, 1.4)
)
CN_mw <- function(th) sapply(mw_tab[-1], function(col) approx(mw_tab$ang, col, min(th, 45))$y)

## ---------------------------------------------------------------------
## 5. AISC 360-22 CAPACITY EQUATIONS
## ---------------------------------------------------------------------
## W & S Shapes Flexural Capacity (AISC F2) - phi*Mn (kip-ft)
phiMn_W <- function(s, Lb_ft) {
  Fy <- s$Fy; Mp <- Fy * s$Zx; Lb <- Lb_ft * 12
  Lp <- 1.76 * s$ry * sqrt(Es / Fy); c <- 1
  Lr <- 1.95 * s$rts * Es / (0.7 * Fy) *
        sqrt(s$J * c / (s$Sx * s$ho) + sqrt((s$J * c / (s$Sx * s$ho))^2 +
             6.76 * (0.7 * Fy / Es)^2))
  Mn <- if (Lb <= Lp) Mp else if (Lb <= Lr)
          Mp - (Mp - 0.7 * Fy * s$Sx) * (Lb - Lp) / (Lr - Lp)
        else min(Mp, (pi^2 * Es / (Lb / s$rts)^2 *
                 sqrt(1 + 0.078 * s$J * c / (s$Sx * s$ho) * (Lb / s$rts)^2)) * s$Sx)
  0.9 * Mn / 12
}

## S-shape Weak-Axis Flexure (AISC F6) - phi*Mny (kip-ft)
phiMny_S <- function(s) {
  0.9 * min(s$Fy * s$Zy, 1.6 * s$Fy * s$Sy) / 12
}

## Square HSS Flexural Capacity (AISC F7) - phi*Mn (kip-ft)
phiMn_H <- function(s) {
  lam <- (s$b - 3 * s$tdes) / s$tdes
  lp  <- 1.12 * sqrt(Es / s$Fy); lr <- 1.40 * sqrt(Es / s$Fy)
  Mp  <- s$Fy * s$Z
  Mn  <- if (lam <= lp) Mp else if (lam <= lr)
           min(Mp, Mp - (Mp - s$Fy * s$S) * (3.57 * lam * sqrt(s$Fy / Es) - 4))
         else NA
  if (is.na(Mn)) return(NA)
  0.9 * Mn / 12
}

## Rectangular HSS 6x2 Flexural Capacity (AISC F7 including LTB) - phi*Mn (kip-ft)
phiMn_R <- function(s, Lb_ft) {
  Fy <- s$Fy; Mp <- Fy * s$Zx
  bt <- (2 - 3 * s$t) / s$t; ht <- (6 - 3 * s$t) / s$t
  lpf <- 1.12 * sqrt(Es / Fy); lrf <- 1.40 * sqrt(Es / Fy)
  lpw <- 2.42 * sqrt(Es / Fy); lrw <- 5.70 * sqrt(Es / Fy)
  Mn <- Mp
  if (bt > lpf) Mn <- min(Mn, Mp - (Mp - Fy * s$Sx) * (3.57 * bt * sqrt(Fy / Es) - 4))
  if (ht > lpw) Mn <- min(Mn, Mp - (Mp - Fy * s$Sx) * (0.305 * ht * sqrt(Fy / Es) - 0.738))
  Lb <- Lb_ft * 12
  Lp <- 0.13 * Es * s$ry * sqrt(s$J * s$A) / Mp
  Lr <- 2 * Es * s$ry * sqrt(s$J * s$A) / (0.7 * Fy * s$Sx)
  if (Lb > Lp) {
    ltb <- if (Lb <= Lr) Mp - (Mp - 0.7 * Fy * s$Sx) * (Lb - Lp) / (Lr - Lp)
           else 2 * Es * sqrt(s$J * s$A) / (Lb / s$ry)
    Mn <- min(Mn, ltb)
  }
  0.9 * Mn / 12
}

## Axial Compression Capacity (AISC E3) - phi*Pn (kips)
phiPn_col <- function(A, r, Fy, KL_ft) {
  sl <- KL_ft * 12 / r
  if (sl > 300) return(NA)  # AISC slenderness cap
  Fe <- pi^2 * Es / sl^2
  Fcr <- if (Fy / Fe <= 2.25) 0.658^(Fy / Fe) * Fy else 0.877 * Fe
  0.9 * Fcr * A
}

## Shear Capacity with exact web thickness (AISC G2) - phi*Vn (kips)
phiVn_shape <- function(d, tw, Fy, type = c("W", "S", "HSS")) {
  type <- match.arg(type)
  if (type %in% c("W", "S")) {
    1.0 * 0.6 * Fy * d * tw
  } else {
    0.9 * 0.6 * Fy * 2 * d * tw
  }
}

## ---------------------------------------------------------------------
## 6. DRILLED PIER FOUNDATION DEPTH (2024 IBC §1807.3.2.1 Formula 18-1)
## ---------------------------------------------------------------------
calc_pier_depth <- function(P_lat_lb, h_lever_ft, b_pier_ft, S_lat = 100, iso_factor = 1.0, frost_in = 30) {
  # Formula 18-1: d = 0.5 * A * (1 + sqrt(1 + 4.36 * h / A))
  # where A = 2.34 * P / (S1 * b)
  d <- 3.0
  for (k in 1:60) {
    S1 <- min(S_lat * iso_factor * d / 3, 15 * S_lat * iso_factor)
    A  <- 2.34 * P_lat_lb / (S1 * b_pier_ft)
    dn <- 0.5 * A * (1 + sqrt(1 + 4.36 * h_lever_ft / A))
    if (abs(dn - d) < 0.01) break
    d <- dn
  }
  max(d, (frost_in + 6) / 12)
}

## ---------------------------------------------------------------------
## 7. SANITY GATING & COMPLIANCE ASSESSMENT
## ---------------------------------------------------------------------
evaluate_compliance <- function(top_ht_ft, mod_rear_psf, mod_front_psf,
                                max_zoning_ht = 15.0, pv = default_pv) {
  zoning_ok  <- top_ht_ft <= max_zoning_ht
  mod_rear_ok<- mod_rear_psf <= pv$rear_design_psf
  mod_front_ok<- mod_front_psf <= pv$front_design_psf
  
  status <- if (!mod_rear_ok) {
    "FAIL: MODULE_UPLIFT"
  } else if (!zoning_ok) {
    "FAIL: EXCEEDS_ZONING_HT"
  } else if (!mod_front_ok) {
    "FAIL: MODULE_SNOW_LIMIT"
  } else {
    "COMPLIANT"
  }
  
  list(zoning_ok = zoning_ok, mod_rear_ok = mod_rear_ok,
       mod_front_ok = mod_front_ok, status = status)
}
