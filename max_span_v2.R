## =====================================================================
## Max Clear Span Study (4-Post Frame) - v2
## Codes: 2024 IBC / ASCE 7-22, AISC 360-22, Utah Code §15A-3-107
##
## Refactored:
##  - Sourced directly from solar_pergola_engine.R (no circular dependencies).
##  - Parameterized elevation / snow basis (defaults to 5,500 ft East Bench;
##    can be toggled for valley floor 4,300 ft).
##  - Evaluates both Cantilever and Knee-Braced upright selections.
##  - Gated by local zoning height (<= 15 ft) and module uplift rating.
## =====================================================================

source("solar_pergola_engine.R")

options(stringsAsFactors = FALSE)
site <- default_site
pv   <- default_pv

# Elevation-based snow loading: Default 5,500 ft (ZIP 84108)
elevation_ft <- if (length(commandArgs(TRUE))) as.numeric(commandArgs(TRUE)[1]) else 5500
snow_info    <- get_snow_loads(elevation_ft, site$Ce, site$Ct)
pf2          <- snow_info$pf
ps2          <- pf2 # conservative horizontal projection

Hmax    <- 16.0  # ft, maximum allowable upright height
s_raf2  <- 3.0   # ft, rafter spacing
D_fr2   <- 2.5   # psf, rafter allowance
pitches <- c(7, 10, 15, 20, 25, 30, 35)
price   <- c(low = 0.75, mid = 1.00, high = 1.25) # $/lb raw material

cat(sprintf("Running Max Span v2 | Elevation: %d ft | Pf: %.1f psf | Zoning Max Ht: %.1f ft\n",
            elevation_ft, pf2, site$zoning_max_ht))

geom2 <- function(Wd, th) {
  t <- th * pi / 180; rise <- Wd * tan(t)
  h <- site$clear_ft + 1.5 + rise / 2
  qk <- qKd(h, site)
  D_h <- (pv$D_psf + rails$wt[1] / pv$rail_trib_ft) / cos(t) + D_fr2
  pmw <- qk * site$G * CN_mw(th)
  ks  <- c("A0", "B0", "A1", "B1")
  rW <- sapply(ks, function(k) Wd * (3 * pmw[paste0(k, "W")] + pmw[paste0(k, "L")]) / 8 / cos(t))
  rL <- sapply(ks, function(k) Wd * (pmw[paste0(k, "W")] + 3 * pmw[paste0(k, "L")]) / 8 / cos(t))
  list(t = t, rise = rise, qk = qk, D_h = D_h, pmw = pmw,
       wdn = max(c(rW, rL), 0), wup = min(c(rW, rL), 0),
       Hcol = site$clear_ft + 1.0 + rise,
       top_ht = site$clear_ft + 1.9 + rise)
}

beam_ok2 <- function(s, L, g, Wd, type) {
  trib <- Wd / 2
  wg <- max(trib * (1.2 * g$D_h + 1.0 * ps2) + 0.5 * g$wdn,
            trib * (1.2 * g$D_h + 0.3 * ps2) + 1.0 * g$wdn) + 1.2 * s$wt
  wu <- trib * 0.9 * g$D_h + 0.9 * s$wt + g$wup
  wS <- trib * 0.7 * ps2
  wT <- trib * (g$D_h + 0.7 * ps2) + s$wt
  Mg <- wg * L^2 / 8 / 1000
  Mu <- abs(min(wu, 0)) * L^2 / 8 / 1000
  Vu <- wg * L / 2 / 1000
  
  if (type == "S") {
    cg <- phiMn_W(s, s_raf2); cu <- phiMn_W(s, L); I <- s$Ix
    Vc <- phiVn_shape(s$d, s$tw, s$Fy, "S")
  } else if (type == "H") {
    cg <- cu <- phiMn_H(s); I <- s$I
    Vc <- phiVn_shape(s$b, s$tdes, s$Fy, "HSS")
  } else if (type == "R") {
    cg <- phiMn_R(s, s_raf2); cu <- phiMn_R(s, L); I <- s$Ix
    Vc <- phiVn_shape(6.0, s$t, s$Fy, "HSS")
  }
  d1 <- 5 * (wS / 12) * (L * 12)^4 / (384 * Es * 1000 * I)
  d2 <- 5 * (wT / 12) * (L * 12)^4 / (384 * Es * 1000 * I)
  Mg <= cg && Mu <= cu && Vu <= Vc && d1 <= L * 12 / 180 && d2 <= L * 12 / 120
}

post_ok_cant <- function(s, L, g, Wd, type, beam_wt) {
  area <- Wd * L; Hc <- g$Hcol
  D_tot <- g$D_h + 2 * beam_wt * L / area
  Wse <- area * (D_tot + snow_info$seismic_coeff * pf2)
  V   <- site$rho * (site$SDS / site$R_cant) * Wse / 4
  Hw  <- max(abs(g$pmw)) * 0.75 * sin(g$t) * area / cos(g$t) + g$qk * 0.85 * 2.0 * 1.5 * max(Wd, L)
  Vw  <- Hw / 4
  Pg  <- area / 4 * (1.2 * D_tot + 1.0 * ps2)
  PE  <- area / 4 * (1.2 * D_tot + 0.15 * ps2) + 0.2 * site$SDS * area / 4 * D_tot
  PW  <- area / 4 * (1.2 * D_tot + 0.3 * ps2)
  ME  <- V * Hc / 1000; MW <- Vw * Hc / 1000; KL <- 2.1 * Hc
  
  if (type == "S") {
    sl <- KL * 12 / s$ry; if (sl > 300) return(FALSE)
    Pc <- phiPn_col(s$A, s$ry, s$Fy, KL); if (is.na(Pc)) return(FALSE)
    Mx <- phiMn_W(s, Hc); My <- phiMny_S(s)
    Pex <- pi^2 * Es * s$Ix / (KL * 12)^2; Pey <- pi^2 * Es * s$Iy / (KL * 12)^2; Imin <- s$Iy
  } else {
    Pc <- phiPn_col(s$A, s$r, s$Fy, KL); if (is.na(Pc)) return(FALSE)
    Mx <- My <- phiMn_H(s)
    Pex <- Pey <- pi^2 * Es * s$I / (KL * 12)^2; Imin <- s$I
  }
  h1 <- function(P, mx, my) {
    Bx <- 1 / max(0.02, 1 - P / 1000 / Pex); By <- 1 / max(0.02, 1 - P / 1000 / Pey)
    r <- P / 1000 / Pc; tm <- Bx * mx / Mx + By * my / My
    if (r >= 0.2) r + 8 / 9 * tm else r / 2 + tm
  }
  ok <- h1(Pg, 0, 0) <= 1 && h1(PE, ME, 0.3 * ME) <= 1 && h1(PE, 0.3 * ME, ME) <= 1 &&
        h1(PW, MW, 0) <= 1 && h1(PW, 0, MW) <= 1
  drift <- V / site$rho * (Hc * 12)^3 / (3 * Es * 1000 * Imin) * site$Cd_cant
  ok && drift <= 0.025 * Hc * 12
}

post_ok_braced <- function(s, L, g, Wd, type, beam_wt) {
  area <- Wd * L; Hc <- g$Hcol
  D_tot <- g$D_h + 2 * beam_wt * L / area
  Wse <- area * (D_tot + snow_info$seismic_coeff * pf2)
  V   <- site$rho * (site$SDS / site$R_braced) * Wse / 4
  Hw  <- max(abs(g$pmw)) * 0.75 * sin(g$t) * area / cos(g$t) + g$qk * 0.85 * 2.0 * 1.5 * max(Wd, L)
  Vw  <- Hw / 4
  Pg  <- area / 4 * (1.2 * D_tot + 1.0 * ps2)
  PE  <- area / 4 * (1.2 * D_tot + 0.15 * ps2) + 0.2 * site$SDS * area / 4 * D_tot
  
  ab_ft   <- 2.5
  H_unbrc <- max(Hc - ab_ft, site$clear_ft)
  ME_b    <- V * ab_ft / 1000
  MW_b    <- Vw * ab_ft / 1000
  KL      <- 1.0 * H_unbrc
  
  if (type == "S") {
    sl <- KL * 12 / s$ry; if (sl > 300) return(FALSE)
    Pc <- phiPn_col(s$A, s$ry, s$Fy, KL); if (is.na(Pc)) return(FALSE)
    Mx <- phiMn_W(s, Lb_ft = H_unbrc); My <- phiMny_S(s)
    Pex <- pi^2 * Es * s$Ix / (KL * 12)^2; Pey <- pi^2 * Es * s$Iy / (KL * 12)^2; Imin <- s$Iy
  } else {
    Pc <- phiPn_col(s$A, s$r, s$Fy, KL); if (is.na(Pc)) return(FALSE)
    Mx <- My <- phiMn_H(s)
    Pex <- Pey <- pi^2 * Es * s$I / (KL * 12)^2; Imin <- s$I
  }
  h1 <- function(P, mx, my) {
    Bx <- 1 / max(0.02, 1 - P / 1000 / Pex); By <- 1 / max(0.02, 1 - P / 1000 / Pey)
    r <- P / 1000 / Pc; tm <- Bx * mx / Mx + By * my / My
    if (r >= 0.2) r + 8 / 9 * tm else r / 2 + tm
  }
  ok <- h1(Pg, 0, 0) <= 1 && h1(PE, ME_b, 0.3 * ME_b) <= 1 &&
        h1(PW <- area / 4 * (1.2 * D_tot + 0.3 * ps2), MW_b, 0) <= 1
  drift <- V / site$rho * (H_unbrc * 12)^3 / (3 * Es * 1000 * Imin) * site$Cd_braced
  ok && drift <= 0.025 * H_unbrc * 12
}

## Candidates
beams_c <- rbind(
  data.frame(shape = S_shapes$shape, wt = S_shapes$wt, type = "S"),
  data.frame(shape = HSS_6x2$shape,  wt = HSS_6x2$wt,  type = "R"),
  data.frame(shape = HSS_sq$shape,   wt = HSS_sq$wt,   type = "H")
)

getS <- function(sh, type) {
  switch(type,
         S = S_shapes[S_shapes$shape == sh, ],
         H = HSS_sq[HSS_sq$shape == sh, ],
         R = HSS_6x2[HSS_6x2$shape == sh, ])
}

posts_c <- rbind(
  data.frame(shape = S_shapes$shape, wt = S_shapes$wt, type = "S"),
  data.frame(shape = HSS_sq$shape,   wt = HSS_sq$wt,   type = "H")
)
posts_c <- posts_c[order(posts_c$wt), ]

Lgrid2 <- seq(4, 45, by = 0.5)

lightest_post <- function(L, g, Wd, bwt, mode = c("cant", "braced")) {
  mode <- match.arg(mode)
  fn <- if (mode == "cant") post_ok_cant else post_ok_braced
  for (j in seq_len(nrow(posts_c))) {
    p <- posts_c[j, ]
    if (fn(getS(p$shape, p$type), L, g, Wd, p$type, bwt)) return(p)
  }
  NULL
}

rows <- list()
cases <- rbind(
  expand.grid(Wd = c(12, 16, 20), th = pitches, square = FALSE),
  expand.grid(Wd = NA, th = pitches, square = TRUE)
)

for (k in seq_len(nrow(cases))) {
  th <- cases$th[k]; sq <- cases$square[k]
  for (b in seq_len(nrow(beams_c))) {
    bc <- beams_c[b, ]; s <- getS(bc$shape, bc$type)
    best <- NA; post_c <- NULL; post_b <- NULL; Wuse <- NA; Hc <- NA; top_h <- NA
    
    for (L in rev(Lgrid2)) {
      Wd <- if (sq) L else cases$Wd[k]
      g  <- geom2(Wd, th)
      if (g$Hcol > Hmax) next
      if (!beam_ok2(s, L, g, Wd, bc$type)) next
      
      pc <- lightest_post(L, g, Wd, bc$wt, "cant")
      pb <- lightest_post(L, g, Wd, bc$wt, "braced")
      if (is.null(pc) && is.null(pb)) next
      
      best <- L; post_c <- pc; post_b <- pb; Wuse <- Wd; Hc <- g$Hcol; top_h <- g$top_ht
      break
    }
    
    if (is.na(best)) {
      rows[[length(rows) + 1]] <- data.frame(
        case = if (sq) "square" else paste0("W=", cases$Wd[k]),
        pitch = th, beam = bc$shape, beam_wt_plf = bc$wt, max_span_ft = NA, W_ft = NA,
        upright_cant = NA, upright_braced = NA, top_ht_ft = NA, zoning_ok = NA,
        steel_lb = NA, cost_mid = NA
      )
      next
    }
    
    rise <- Wuse * tan(th * pi / 180)
    p_use <- if (!is.null(post_b)) post_b else post_c
    lb <- bc$wt * (2 * best + 2 * Wuse / cos(th * pi / 180)) +
          p_use$wt * (2 * (Hc - rise + 2) + 2 * (Hc + 2))
    
    rows[[length(rows) + 1]] <- data.frame(
      case = if (sq) "square" else paste0("W=", cases$Wd[k]),
      pitch = th, beam = bc$shape, beam_wt_plf = bc$wt, max_span_ft = best, W_ft = Wuse,
      upright_cant   = if (!is.null(post_c)) post_c$shape else "NONE (KL/r>300)",
      upright_braced = if (!is.null(post_b)) post_b$shape else "NONE",
      top_ht_ft      = round(top_h, 1),
      zoning_ok      = top_h <= site$zoning_max_ht,
      steel_lb       = round(lb),
      cost_mid       = round(lb * price["mid"])
    )
  }
}

res_v2 <- do.call(rbind, rows)
write.csv(res_v2, "max_span_v2.csv", row.names = FALSE)
cat("Completed max_span_v2.R. Results written to max_span_v2.csv.\n")
