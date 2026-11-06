## =====================================================================
## Maximum Post Spacing Study - Solar Pergola, Salt Lake City (ZIP 84108)
## Codes: 2024 IBC / ASCE 7-22, AISC 360-22, Utah Code §15A-3-107
##
## Refactored:
##  - Sourced directly from solar_pergola_engine.R (no circular dependencies).
##  - Uses elevation-harmonized snow load.
##  - Uses exact AISC tw web thickness for S-shape shear.
##  - Computes BOTH Unbraced Cantilever post spacing and Knee-Braced post spacing.
##    (Demonstrates how knee braces unlock S-shapes as uprights).
## =====================================================================

source("solar_pergola_engine.R")

options(stringsAsFactors = FALSE)
site <- default_site
pv   <- default_pv

# Elevation-harmonized snow load (5,500 ft East Bench)
snow_info <- get_snow_loads(site$elevation_ft, site$Ce, site$Ct)
pf <- snow_info$pf

Lgrid  <- seq(6, 40, by = 0.5)
widths <- c(12, 16, 20)
pitches <- c(7, 10, 15, 20, 25, 30, 35)
s_raf  <- 3.0   # ft, rafter spacing (XR10)
D_fr   <- 3.0   # psf framing allowance (rafters)

geom <- function(Wd, th) {
  t <- th * pi / 180; rise <- Wd * tan(t)
  h <- site$clear_ft + 1.5 + rise / 2
  qk <- qKd(h, site)
  D_h <- (pv$D_psf + rails$wt[1] / pv$rail_trib_ft) / cos(t) + D_fr
  S_h <- Cs_slip(th) * pf
  pmw <- qk * site$G * CN_mw(th)
  ks  <- c("A0", "B0", "A1", "B1")
  rW <- sapply(ks, function(k) Wd * (3 * pmw[paste0(k, "W")] + pmw[paste0(k, "L")]) / 8 / cos(t))
  rL <- sapply(ks, function(k) Wd * (pmw[paste0(k, "W")] + 3 * pmw[paste0(k, "L")]) / 8 / cos(t))
  list(t = t, rise = rise, qk = qk, D_h = D_h, S_h = S_h, pmw = pmw,
       wdn = max(c(rW, rL), 0), wup = min(c(rW, rL), 0),
       Hcol = site$clear_ft + 1.0 + rise)
}

## Girder (perimeter beam) capacity check
girder_ok <- function(s, L, g, Wd, type) {
  wg_g <- max((Wd / 2) * (1.2 * g$D_h + 1.0 * g$S_h) + 0.5 * g$wdn,
              (Wd / 2) * (1.2 * g$D_h + 0.3 * g$S_h) + 1.0 * g$wdn) + 1.2 * s$wt
  wg_u <- (Wd / 2) * 0.9 * g$D_h + 0.9 * s$wt + g$wup
  wS <- (Wd / 2) * 0.7 * g$S_h
  wT <- (Wd / 2) * (g$D_h + 0.7 * g$S_h) + s$wt
  Mg <- wg_g * L^2 / 8 / 1000
  Mu <- abs(min(wg_u, 0)) * L^2 / 8 / 1000
  
  if (type == "S") {
    capg <- phiMn_W(s, s_raf); capu <- phiMn_W(s, L); I <- s$Ix
    Vcap <- phiVn_shape(s$d, s$tw, s$Fy, "S")
  } else {
    capg <- capu <- phiMn_H(s); I <- s$I
    Vcap <- phiVn_shape(s$b, s$tdes, s$Fy, "HSS")
  }
  Vu <- wg_g * L / 2 / 1000
  d1 <- 5 * (wS / 12) * (L * 12)^4 / (384 * Es * 1000 * I)
  d2 <- 5 * (wT / 12) * (L * 12)^4 / (384 * Es * 1000 * I)
  
  Mg <= capg && Mu <= capu && Vu <= Vcap && d1 <= L * 12 / 180 && d2 <= L * 12 / 120
}

## 1. Unbraced Cantilever Post Check (R = 1.25, K = 2.1)
post_cant_ok <- function(s, L, g, Wd, type, girder_wt = 20) {
  area <- Wd * L; Hc <- g$Hcol
  D_tot <- g$D_h + 2 * girder_wt * L / area
  Wseis <- area * (D_tot + snow_info$seismic_coeff * pf)
  V <- site$rho * (site$SDS / site$R_cant) * Wseis / 4
  Hw <- max(abs(g$pmw)) * 0.75 * sin(g$t) * area / cos(g$t) +
        g$qk * 0.85 * 2.0 * 1.5 * max(Wd, L)
  Vw <- Hw / 4
  Pg <- area / 4 * (1.2 * D_tot + 1.0 * g$S_h)
  PE <- area / 4 * (1.2 * D_tot + 0.15 * g$S_h) + 0.2 * site$SDS * area / 4 * D_tot
  PW <- area / 4 * (1.2 * D_tot + 0.3 * g$S_h)
  ME <- V * Hc / 1000; MW <- Vw * Hc / 1000
  KL <- 2.1 * Hc
  
  if (type == "S") {
    rmin <- s$ry; Imin <- s$Iy; Imaj <- s$Ix
    sl <- KL * 12 / rmin
    if (sl > 300) return(FALSE)
    Pc <- phiPn_col(s$A, rmin, s$Fy, KL)
    if (is.na(Pc)) return(FALSE)
    Mx <- phiMn_W(s, Hc); My <- phiMny_S(s)
    Pex <- pi^2 * Es * Imaj / (KL * 12)^2; Pey <- pi^2 * Es * Imin / (KL * 12)^2
  } else {
    Pc <- phiPn_col(s$A, s$r, s$Fy, KL)
    if (is.na(Pc)) return(FALSE)
    Mx <- My <- phiMn_H(s)
    Pex <- Pey <- pi^2 * Es * s$I / (KL * 12)^2; Imin <- Imaj <- s$I
  }
  h1 <- function(P, mx, my) {
    Bx <- 1 / max(0.02, 1 - P / 1000 / Pex); By <- 1 / max(0.02, 1 - P / 1000 / Pey)
    r <- P / 1000 / Pc
    term <- Bx * mx / Mx + By * my / My
    if (r >= 0.2) r + 8 / 9 * term else r / 2 + term
  }
  ok <- h1(Pg, 0, 0) <= 1 &&
        h1(PE, ME, 0.3 * ME) <= 1 && h1(PE, 0.3 * ME, ME) <= 1 &&
        h1(PW, MW, 0) <= 1 && h1(PW, 0, MW) <= 1
  drift <- V / site$rho * (Hc * 12)^3 / (3 * Es * 1000 * Imin) * site$Cd_cant
  ok && drift <= 0.025 * Hc * 12
}

## 2. Knee-Braced Post Check (ab = 2.5 ft, R = 3.25, K = 1.0)
post_braced_ok <- function(s, L, g, Wd, type, girder_wt = 20) {
  area <- Wd * L; Hc <- g$Hcol
  D_tot <- g$D_h + 2 * girder_wt * L / area
  Wseis <- area * (D_tot + snow_info$seismic_coeff * pf)
  V <- site$rho * (site$SDS / site$R_braced) * Wseis / 4
  Hw <- max(abs(g$pmw)) * 0.75 * sin(g$t) * area / cos(g$t) +
        g$qk * 0.85 * 2.0 * 1.5 * max(Wd, L)
  Vw <- Hw / 4
  Pg <- area / 4 * (1.2 * D_tot + 1.0 * g$S_h)
  PE <- area / 4 * (1.2 * D_tot + 0.15 * g$S_h) + 0.2 * site$SDS * area / 4 * D_tot
  
  ab_ft <- 2.5
  H_unbrc <- max(Hc - ab_ft, site$clear_ft)
  ME_b <- V * ab_ft / 1000
  MW_b <- Vw * ab_ft / 1000
  KL <- 1.0 * H_unbrc
  
  if (type == "S") {
    rmin <- s$ry; Imin <- s$Iy; Imaj <- s$Ix
    sl <- KL * 12 / rmin
    if (sl > 300) return(FALSE)
    Pc <- phiPn_col(s$A, rmin, s$Fy, KL)
    if (is.na(Pc)) return(FALSE)
    Mx <- phiMn_W(s, Lb_ft = H_unbrc); My <- phiMny_S(s)
    Pex <- pi^2 * Es * Imaj / (KL * 12)^2; Pey <- pi^2 * Es * Imin / (KL * 12)^2
  } else {
    Pc <- phiPn_col(s$A, s$r, s$Fy, KL)
    if (is.na(Pc)) return(FALSE)
    Mx <- My <- phiMn_H(s)
    Pex <- Pey <- pi^2 * Es * s$I / (KL * 12)^2; Imin <- Imaj <- s$I
  }
  h1 <- function(P, mx, my) {
    Bx <- 1 / max(0.02, 1 - P / 1000 / Pex); By <- 1 / max(0.02, 1 - P / 1000 / Pey)
    r <- P / 1000 / Pc
    term <- Bx * mx / Mx + By * my / My
    if (r >= 0.2) r + 8 / 9 * term else r / 2 + term
  }
  ok <- h1(Pg, 0, 0) <= 1 &&
        h1(PE, ME_b, 0.3 * ME_b) <= 1 &&
        h1(PW <- area / 4 * (1.2 * D_tot + 0.3 * g$S_h), MW_b, 0) <= 1
  drift <- V / site$rho * (H_unbrc * 12)^3 / (3 * Es * 1000 * Imin) * site$Cd_braced
  ok && drift <= 0.025 * H_unbrc * 12
}

find_maxL <- function(fun, s, g, Wd, type) {
  okL <- Lgrid[vapply(Lgrid, function(L) fun(s, L, g, Wd, type), logical(1))]
  if (length(okL)) max(okL) else NA
}

out <- list()
for (Wd in widths) for (th in pitches) {
  g <- geom(Wd, th)
  top_ht <- site$clear_ft + 1.9 + g$rise
  zoning_ok <- top_ht <= site$zoning_max_ht
  
  for (i in seq_len(nrow(S_shapes))) {
    s <- S_shapes[i, ]
    out[[length(out) + 1]] <- data.frame(
      W_ft = Wd, pitch = th, shape = s$shape, type = "S",
      girder_maxL = find_maxL(girder_ok, s, g, Wd, "S"),
      post_cant_maxL = find_maxL(post_cant_ok, s, g, Wd, "S"),
      post_braced_maxL = find_maxL(post_braced_ok, s, g, Wd, "S"),
      post_ht = round(g$Hcol, 1), top_ht_ft = round(top_ht, 1),
      zoning_ok = zoning_ok
    )
  }
  for (i in seq_len(nrow(HSS_sq))) {
    s <- HSS_sq[i, ]
    out[[length(out) + 1]] <- data.frame(
      W_ft = Wd, pitch = th, shape = s$shape, type = "HSS",
      girder_maxL = find_maxL(girder_ok, s, g, Wd, "H"),
      post_cant_maxL = find_maxL(post_cant_ok, s, g, Wd, "H"),
      post_braced_maxL = find_maxL(post_braced_ok, s, g, Wd, "H"),
      post_ht = round(g$Hcol, 1), top_ht_ft = round(top_ht, 1),
      zoning_ok = zoning_ok
    )
  }
}

res <- do.call(rbind, out)
write.csv(res, "max_post_spacing.csv", row.names = FALSE)
cat("Completed max_post_spacing.R. Results written to max_post_spacing.csv.\n")
