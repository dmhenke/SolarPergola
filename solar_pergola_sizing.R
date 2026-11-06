## =====================================================================
## Solar Pergola Preliminary Sizing  -  Freestanding 4-Post Steel Frame
## Site envelope: Salt Lake City, UT (ZIP 84108)
## Codes: 2024 IBC / ASCE 7-22, AISC 360-22, ADM 2020, Utah Code §15A-3-107
##
## Refactored:
##  - Sourced directly from solar_pergola_engine.R (no circular dependencies).
##  - Harmonized elevation-based snow load.
##  - Compares Unbraced Cantilever vs Knee-Braced structural framing.
##  - Enforces local zoning height gates and PV module uplift limits.
## =====================================================================

source("solar_pergola_engine.R")

options(width = 160, stringsAsFactors = FALSE)
out_dir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."

site <- default_site
pv   <- default_pv

# Compute site snow loads based on elevation (5,500 ft East Bench)
snow_info <- get_snow_loads(site$elevation_ft, site$Ce, site$Ct)
pf <- snow_info$pf
cat(sprintf("Site Elevation: %d ft | Pg: %.1f psf | Pf (strength): %.1f psf | Pf (ASD): %.1f psf\n",
            site$elevation_ft, snow_info$pg, snow_info$pf, snow_info$pf_asd))

sizes <- data.frame(
  name = c("S 12x16", "M 16x20", "L 20x24"),
  W = c(12, 16, 20),  # Slope direction (rafter span), ft
  L = c(16, 20, 24)   # Girder span / post spacing, ft
)
pitches <- c(7, 10, 15, 20, 25, 30, 35)

## ---------------------------------------------------------------------
## BEAM PICKER
## ---------------------------------------------------------------------
pick_beam <- function(span, w_grav_lrfd, w_up_lrfd, w_svc_S, w_svc_tot,
                      Lb_grav, kind = c("W", "HSS")) {
  kind <- match.arg(kind)
  Mu_g <- w_grav_lrfd * span^2 / 8 / 1000
  Mu_u <- abs(w_up_lrfd) * span^2 / 8 / 1000
  tab <- if (kind == "W") W_shapes[order(W_shapes$wt), ] else HSS_sq[order(HSS_sq$wt), ]
  
  for (i in seq_len(nrow(tab))) {
    s <- tab[i, ]
    sw <- s$wt
    Mu_g2 <- Mu_g + 1.2 * sw * span^2 / 8 / 1000
    if (kind == "W") {
      ok_g <- Mu_g2 <= phiMn_W(s, Lb_grav)
      ok_u <- Mu_u <= phiMn_W(s, span) # bottom flange unbraced in uplift
      I <- s$Ix
    } else {
      cap <- phiMn_H(s)
      if (is.na(cap)) next
      ok_g <- Mu_g2 <= cap
      ok_u <- Mu_u <= cap
      I <- s$I
    }
    dS <- 5 * (w_svc_S / 12) * (span * 12)^4 / (384 * Es * 1000 * I)
    dT <- 5 * ((w_svc_tot + sw) / 12) * (span * 12)^4 / (384 * Es * 1000 * I)
    if (ok_g && ok_u && dS <= span * 12 / 180 && dT <= span * 12 / 120) {
      return(list(shape = s$shape, wt = sw, Mu = round(Mu_g2, 1), defl_S = round(dS, 2)))
    }
  }
  list(shape = "NONE", wt = NA, Mu = round(Mu_g, 1), defl_S = NA)
}

## ---------------------------------------------------------------------
## RAIL & MODULE CHECKS
## ---------------------------------------------------------------------
rail_span <- function(th, rail, a_ft, h_ft) {
  qk <- qKd(h_ft, site)
  c2 <- cos(th * pi / 180)
  best <- NA
  for (Lin in seq(96, 12, by = -2)) {
    Lft <- Lin / 12
    EA  <- max(Lft * pv$rail_trib_ft, Lft^2 / 3)
    cn  <- CN_cc(th, EA, a_ft)
    pW_dn <- qk * site$G * cn["pos"]
    pW_up <- qk * site$G * cn["neg"]
    D  <- (pv$D_psf + rail$wt / pv$rail_trib_ft) * c2
    S  <- 0.7 * (Cs_slip(th) * pf) * c2^2
    w_dn <- max(D + S, D + 0.6 * pW_dn, D + 0.75 * 0.6 * pW_dn + 0.75 * S)
    w_up <- 0.6 * D + 0.6 * pW_up
    w <- max(w_dn, abs(w_up)) * pv$rail_trib_ft
    M_ok <- (w * Lft^2 / 10) <= rail$Ma_lbft
    d_lim <- Lin / 180 # tighter deflection for PV integrity
    defl <- (w / 12) * Lin^4 / (145 * E_al * 1000 * rail$Ix)
    if (M_ok && defl <= d_lim) { best <- Lin; break }
  }
  best
}

module_check <- function(th, a_ft, h_ft) {
  qk <- qKd(h_ft, site)
  c2 <- cos(th * pi / 180)
  cn <- CN_cc(th, pv$area_ft2, a_ft)
  front <- max(pv$D_psf * c2 + 0.7 * (Cs_slip(th) * pf) * c2^2,
               pv$D_psf * c2 + 0.75 * (0.6 * qk * site$G * cn["pos"] + 0.7 * (Cs_slip(th) * pf) * c2^2))
  rear  <- abs(0.6 * pv$D_psf * c2 + 0.6 * qk * site$G * cn["neg"])
  c(front = unname(front), rear = unname(rear),
    front_ok = front <= pv$front_design_psf, rear_ok = rear <= pv$rear_design_psf)
}

## ---------------------------------------------------------------------
## RUN ANALYSIS
## ---------------------------------------------------------------------
results <- list()
railtab <- list()

for (i in seq_len(nrow(sizes))) for (th in pitches) {
  sz <- sizes[i, ]
  t  <- th * pi / 180
  rise   <- sz$W * tan(t)
  h_mean <- site$clear_ft + 1.5 + rise / 2
  least  <- min(sz$W, sz$L)
  a_ft   <- max(min(0.1 * least, 0.4 * h_mean), 0.04 * least, 3)
  qk     <- qKd(h_mean, site)
  
  # Rails & Modules
  rs <- sapply(seq_len(nrow(rails)), function(r) rail_span(th, rails[r, ], a_ft, h_mean))
  mc <- module_check(th, a_ft, h_mean)
  
  top_ht <- site$clear_ft + 1.9 + rise
  gate <- evaluate_compliance(top_ht, mc["rear"], mc["front"], site$zoning_max_ht, pv)
  
  railtab[[length(railtab) + 1]] <- data.frame(
    size = sz$name, pitch = th,
    XR10_in = rs[1], XR100_in = rs[2], XR1000_in = rs[3],
    mod_front_psf = round(mc["front"], 1), mod_rear_psf = round(mc["rear"], 1),
    mod_front_ok = as.logical(mc["front_ok"]), mod_rear_ok = as.logical(mc["rear_ok"]),
    zoning_ok = gate$zoning_ok, status = gate$status
  )
  
  # Rafters (spans W)
  s_in  <- min(rs[1], 72)
  s_ft  <- s_in / 12
  D_pv  <- pv$D_psf + rails$wt[1] / pv$rail_trib_ft
  D_h   <- D_pv / cos(t)
  S_h   <- Cs_slip(th) * pf
  cnr   <- CN_cc(th, max(sz$W / cos(t) * s_ft, (sz$W / cos(t))^2 / 3), a_ft)
  Wdn   <- qk * site$G * cnr["pos"]
  Wup   <- qk * site$G * cnr["neg"]
  
  wr_g  <- s_ft * max(1.2 * D_h + 1.0 * S_h + 0.5 * Wdn,
                      1.2 * D_h + 1.0 * Wdn + 0.3 * S_h)
  wr_u  <- s_ft * (0.9 * D_h + 1.0 * Wup / cos(t))
  wr_S  <- s_ft * 0.7 * S_h
  wr_T  <- s_ft * (D_h + 0.7 * S_h)
  rafW  <- pick_beam(sz$W, wr_g, wr_u, wr_S, wr_T, Lb_grav = pv$rail_trib_ft, "W")
  rafH  <- pick_beam(sz$W, wr_g, wr_u, wr_S, wr_T, Lb_grav = pv$rail_trib_ft, "HSS")
  
  # XR1000 alternative rafter spacing
  s2_ft <- min(rs[3], 72) / 12
  k2    <- s2_ft / s_ft
  rafW2 <- pick_beam(sz$W, wr_g * k2, wr_u * k2, wr_S * k2, wr_T * k2, pv$rail_trib_ft, "W")
  rafH2 <- pick_beam(sz$W, wr_g * k2, wr_u * k2, wr_S * k2, wr_T * k2, pv$rail_trib_ft, "HSS")
  
  # Girder (spans L)
  cm  <- CN_mw(th)
  pmw <- qk * site$G * cm
  rW  <- sapply(c("A0", "B0", "A1", "B1"), function(k)
           sz$W * (3 * pmw[paste0(k, "W")] + pmw[paste0(k, "L")]) / 8 / cos(t))
  rL  <- sapply(c("A0", "B0", "A1", "B1"), function(k)
           sz$W * (pmw[paste0(k, "W")] + 3 * pmw[paste0(k, "L")]) / 8 / cos(t))
  wg_wdn <- max(c(rW, rL), 0)
  wg_wup <- min(c(rW, rL), 0)
  Dfr_h  <- D_h + rafW$wt / s_ft
  wg_g   <- max((sz$W / 2) * (1.2 * Dfr_h + 1.0 * S_h) + 0.5 * wg_wdn,
                (sz$W / 2) * (1.2 * Dfr_h + 0.3 * S_h) + 1.0 * wg_wdn)
  wg_u   <- (sz$W / 2) * 0.9 * Dfr_h + wg_wup
  wg_S   <- (sz$W / 2) * 0.7 * S_h
  wg_T   <- (sz$W / 2) * (Dfr_h + 0.7 * S_h)
  girW   <- pick_beam(sz$L, wg_g, wg_u, wg_S, wg_T, Lb_grav = s_ft, "W")
  girH   <- pick_beam(sz$L, wg_g, wg_u, wg_S, wg_T, Lb_grav = s_ft, "HSS")
  
  # Lateral & Uprights
  area  <- sz$W * sz$L
  D_tot <- Dfr_h + 2 * girW$wt * sz$L / area
  # Utah Code 15A-3-107 seismic snow factor
  Wseis <- area * (D_tot + snow_info$seismic_coeff * pf)
  
  # 1. Unbraced Cantilever Column Mode (R = 1.25)
  Cs_cant   <- site$SDS / site$R_cant
  Vcol_E_c  <- site$rho * Cs_cant * Wseis / 4
  pn_max    <- max(abs(pmw)) * 0.75
  Hwind     <- pn_max * sin(t) * area / cos(t) + qk * 0.85 * 2.0 * 1.5 * max(sz$W, sz$L)
  Vcol_W    <- Hwind / 4
  Hcol      <- site$clear_ft + 1.0 + rise
  Mu_E_c    <- Vcol_E_c * Hcol / 1000
  Mu_W_c    <- Vcol_W * Hcol / 1000
  Pu_g      <- area / 4 * (1.2 * D_tot + 1.0 * S_h)
  Pu_E      <- area / 4 * (1.2 * D_tot + 0.15 * S_h) + 0.2 * site$SDS * area / 4 * D_tot
  
  col_cant <- "NONE"; colwt_cant <- NA; drift_c <- NA
  for (j in order(HSS_sq$wt)) {
    s <- HSS_sq[j, ]; cap <- phiMn_H(s); if (is.na(cap)) next
    Pc <- phiPn_col(s$A, s$r, s$Fy, 2.1 * Hcol); if (is.na(Pc)) next
    Pe <- pi^2 * Es * s$I / (2.1 * Hcol * 12)^2
    chk <- function(P, M, M2 = 0) {
      B <- 1 / max(0.05, 1 - P / 1000 / Pe)
      r <- P / 1000 / Pc
      if (r >= 0.2) r + 8 / 9 * (B * M + B * M2) / cap else r / 2 + (B * M + B * M2) / cap
    }
    ok <- chk(Pu_g, 0) <= 1 && chk(Pu_E, Mu_E_c, 0.3 * Mu_E_c) <= 1 &&
          chk(area / 4 * (1.2 * D_tot + 0.3 * S_h), Mu_W_c) <= 1
    dE <- Vcol_E_c / site$rho * (Hcol * 12)^3 / (3 * Es * 1000 * s$I) * site$Cd_cant
    if (ok && dE <= 0.025 * Hcol * 12) { col_cant <- s$shape; colwt_cant <- s$wt; drift_c <- dE; break }
  }
  
  # Cantilever Pier Foundation
  P_asd_c <- max(0.7 * site$Omega0_cant * Vcol_E_c, 0.6 * Vcol_W)
  M_asd_c <- max(0.7 * site$Omega0_cant * Mu_E_c, 0.6 * Mu_W_c) * 1000
  hP_c    <- M_asd_c / P_asd_c
  bear_dia <- sqrt(4 * (area / 4 * (D_tot + 0.7 * S_h)) / (pi * site$q_bear))
  b_ft    <- max(2.0, ceiling(bear_dia * 2) / 2)
  dp_cant <- calc_pier_depth(P_asd_c, hP_c, b_ft, S_lat = site$S_lat, iso_factor = 1.0, frost_in = site$frost_in)
  
  # 2. Knee-Braced Column Mode (Brace arm ab = 2.5 ft, R = 3.25)
  ab_ft     <- 2.5
  Cs_brc    <- site$SDS / site$R_braced
  Vcol_E_b  <- site$rho * Cs_brc * Wseis / 4
  H_unbrc   <- max(Hcol - ab_ft, site$clear_ft)
  Mu_E_b    <- Vcol_E_b * ab_ft / 1000  # localized couple moment at brace
  Mu_W_b    <- Vcol_W * ab_ft / 1000
  
  # Check HSS and S-shapes for braced column
  col_brc <- "NONE"; col_brc_S <- "NONE"
  for (j in order(HSS_sq$wt)) {
    s <- HSS_sq[j, ]; cap <- phiMn_H(s); if (is.na(cap)) next
    Pc <- phiPn_col(s$A, s$r, s$Fy, 1.0 * H_unbrc); if (is.na(Pc)) next
    Pe <- pi^2 * Es * s$I / (1.0 * H_unbrc * 12)^2
    chk <- function(P, M) {
      B <- 1 / max(0.05, 1 - P / 1000 / Pe)
      r <- P / 1000 / Pc
      if (r >= 0.2) r + 8 / 9 * (B * M) / cap else r / 2 + (B * M) / cap
    }
    if (chk(Pu_g, 0) <= 1 && chk(Pu_E, Mu_E_b) <= 1 && chk(Pu_g, Mu_W_b) <= 1) {
      col_brc <- s$shape; break
    }
  }
  for (j in order(S_shapes$wt)) {
    s <- S_shapes[j, ]
    cap <- phiMn_W(s, Lb_ft = H_unbrc)
    Pc  <- phiPn_col(s$A, s$ry, s$Fy, 1.0 * H_unbrc); if (is.na(Pc)) next
    Pe  <- pi^2 * Es * s$Iy / (1.0 * H_unbrc * 12)^2
    chk <- function(P, M) {
      B <- 1 / max(0.05, 1 - P / 1000 / Pe)
      r <- P / 1000 / Pc
      if (r >= 0.2) r + 8 / 9 * (B * M) / cap else r / 2 + (B * M) / cap
    }
    if (chk(Pu_g, 0) <= 1 && chk(Pu_E, Mu_E_b) <= 1 && chk(Pu_g, Mu_W_b) <= 1) {
      col_brc_S <- s$shape; break
    }
  }
  
  # Braced pier: Base is pinned/semi-rigid, lateral force at ground line (h = 0.5 ft)
  P_asd_b <- max(0.7 * Vcol_E_b, 0.6 * Vcol_W)
  dp_brc  <- calc_pier_depth(P_asd_b, h_lever_ft = 0.5, b_pier_ft = 2.0, S_lat = site$S_lat, frost_in = site$frost_in)

  results[[length(results) + 1]] <- data.frame(
    size = sz$name, pitch = th, rise_ft = round(rise, 1),
    top_ht_ft = round(top_ht, 1), zoning_ok = gate$zoning_ok,
    mod_rear_ok = gate$mod_rear_ok, status = gate$status,
    ps_psf = round(S_h, 1), qh_psf = round(qk, 1),
    rafter_sp_in = s_in, rafter_W = rafW$shape, rafter_HSS = rafH$shape,
    rafter2_W = rafW2$shape, rafter2_HSS = rafH2$shape,
    girder_W = girW$shape, girder_HSS = girH$shape,
    # Unbraced cantilever system
    col_cant_HSS = col_cant, pier_cant_depth_ft = round(dp_cant, 1),
    # Knee-braced system
    col_braced_HSS = col_brc, col_braced_S = col_brc_S,
    pier_braced_depth_ft = round(dp_brc, 1), pier_dia_in = b_ft * 12
  )
}

rail_df <- do.call(rbind, railtab)
res_df  <- do.call(rbind, results)

write.csv(rail_df, file.path(out_dir, "rail_module_checks.csv"), row.names = FALSE)
write.csv(res_df,  file.path(out_dir, "frame_sizing.csv"),        row.names = FALSE)

cat("\nAnalysis complete. Results written to frame_sizing.csv and rail_module_checks.csv.\n")
