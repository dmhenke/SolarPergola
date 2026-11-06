## =====================================================================
## Perimeter Beam Check for HSS6x2 (2x6 rectangular tube, 6 in vertical)
## Codes: 2024 IBC / ASCE 7-22, AISC 360-22, Utah Code §15A-3-107
##
## Refactored:
##  - Sourced directly from solar_pergola_engine.R (no circular dependencies).
##  - Uses elevation-harmonized snow load.
##  - Evaluates both unbraced girder spans and knee-braced spans.
##  - Outputs beam_max_span.csv, upright_max_spacing.csv, and hss6x2_properties.csv.
## =====================================================================

source("solar_pergola_engine.R")

options(stringsAsFactors = FALSE)
site <- default_site
pv   <- default_pv

# Elevation-harmonized snow load (5,500 ft East Bench)
snow_info <- get_snow_loads(site$elevation_ft, site$Ce, site$Ct)
pf <- snow_info$pf

widths  <- c(12, 16, 20)
pitches <- c(7, 10, 15, 20, 25, 30, 35)
s_raf   <- 3.0 # ft, rafter spacing
D_fr    <- 3.0 # psf

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

beam_ok_R <- function(s, L, g, Wd, braced = FALSE, ab_ft = 2.5) {
  L_eff <- if (braced) max(L - 2 * ab_ft, 4.0) else L
  wg_g <- max((Wd / 2) * (1.2 * g$D_h + 1.0 * g$S_h) + 0.5 * g$wdn,
              (Wd / 2) * (1.2 * g$D_h + 0.3 * g$S_h) + 1.0 * g$wdn) + 1.2 * s$wt
  wg_u <- (Wd / 2) * 0.9 * g$D_h + 0.9 * s$wt + g$wup
  wS <- (Wd / 2) * 0.7 * g$S_h
  wT <- (Wd / 2) * (g$D_h + 0.7 * g$S_h) + s$wt
  
  Mg <- wg_g * L_eff^2 / 8 / 1000
  Mu <- abs(min(wg_u, 0)) * L_eff^2 / 8 / 1000
  Vu <- wg_g * L / 2 / 1000
  Vcap <- phiVn_shape(6.0, s$t, s$Fy, "HSS")
  
  d1 <- 5 * (wS / 12) * (L_eff * 12)^4 / (384 * Es * 1000 * s$Ix)
  d2 <- 5 * (wT / 12) * (L_eff * 12)^4 / (384 * Es * 1000 * s$Ix)
  
  Mg <= phiMn_R(s, s_raf) && Mu <= phiMn_R(s, L_eff) && Vu <= Vcap &&
    d1 <= L_eff * 12 / 180 && d2 <= L_eff * 12 / 120
}

Lfine <- seq(3, 40, by = 0.5)
res_R <- list()

for (Wd in widths) for (th in pitches) {
  g <- geom(Wd, th)
  top_ht <- site$clear_ft + 1.9 + g$rise
  zoning_ok <- top_ht <= site$zoning_max_ht
  
  for (i in seq_len(nrow(HSS_6x2))) {
    s <- HSS_6x2[i, ]
    okL_unbrc <- Lfine[vapply(Lfine, function(L) beam_ok_R(s, L, g, Wd, braced = FALSE), logical(1))]
    okL_brc   <- Lfine[vapply(Lfine, function(L) beam_ok_R(s, L, g, Wd, braced = TRUE), logical(1))]
    
    res_R[[length(res_R) + 1]] <- data.frame(
      W_ft = Wd, pitch = th, shape = s$shape,
      beam_unbraced_maxL = if (length(okL_unbrc)) max(okL_unbrc) else NA,
      beam_braced_maxL   = if (length(okL_brc)) max(okL_brc) else NA,
      top_ht_ft = round(top_ht, 1), zoning_ok = zoning_ok
    )
  }
}
res_R_df <- do.call(rbind, res_R)

# Read max_post_spacing.csv to integrate master beam and upright outputs
post_data <- read.csv("max_post_spacing.csv")

# Combine HSS6x2 into beam_max_span.csv
beams_all <- rbind(
  data.frame(W_ft = post_data$W_ft, pitch = post_data$pitch, shape = post_data$shape,
             beam_maxL = post_data$girder_maxL, top_ht_ft = post_data$top_ht_ft, zoning_ok = post_data$zoning_ok),
  data.frame(W_ft = res_R_df$W_ft, pitch = res_R_df$pitch, shape = res_R_df$shape,
             beam_maxL = res_R_df$beam_unbraced_maxL, top_ht_ft = res_R_df$top_ht_ft, zoning_ok = res_R_df$zoning_ok)
)

uprights_all <- data.frame(
  W_ft = post_data$W_ft, pitch = post_data$pitch, shape = post_data$shape,
  post_ht = post_data$post_ht, post_cant_maxL = post_data$post_cant_maxL,
  post_braced_maxL = post_data$post_braced_maxL, top_ht_ft = post_data$top_ht_ft,
  zoning_ok = post_data$zoning_ok
)

write.csv(beams_all, "beam_max_span.csv", row.names = FALSE)
write.csv(uprights_all, "upright_max_spacing.csv", row.names = FALSE)
write.csv(HSS_6x2, "hss6x2_properties.csv", row.names = FALSE)

cat("Completed hss6x2_beam.R. Updated beam_max_span.csv, upright_max_spacing.csv, and hss6x2_properties.csv.\n")
