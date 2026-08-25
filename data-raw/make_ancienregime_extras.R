# ---------------------------------------------------------------------
# Generator for the two companion datasets to `ancienregime`:
#
#   * ancienregime_parishes  -- men of the Third Estate linked across the
#     (fictional) capitation rolls of 1750, 1770, and 1780, with parish
#     of residence recorded each time. 240 parishes in the 12 provinces
#     of the ancienregime world: the sparse-outcome setting. Link
#     failures draw a rival from the man's BIRTH province's roll,
#     because the intendance's clerks searched only their own
#     generalite's registers -- which is what makes the rival law
#     local, and the unconditional model wrong in an instructive way.
#
#   * ancienregime_lineages -- father-son pairs: the father recorded in
#     the 1750 capitation and again (by a link) in the 1745 dixieme; the
#     son linked into the rolls of 1770, 1780, and 1790. Three linked
#     measures of the son are what identify DEPENDENT link failures
#     (with probability pi_s, every one of a son's failed links lands on
#     the same rival), and the clerk's miscoding of occupations (the
#     measurement kernel, share mu) is a second error the linkage model
#     alone cannot see.
#
#   * ancienregime_truth -- the parameters both generators used, plus
#     the parish register and the auxiliary tables an analyst would
#     bring (the vingtieme's recalled-residence flows, the occupation
#     margins). Shipping the truth keeps the vignette honest: every
#     "recovered" number is compared against it, and the tests assert
#     the comparison.
#
# The data-generating processes match the estimators' models EXACTLY --
# rho in latent space, the measurement kernel applied to every emission
# whatever its source, the father's correct link T^dstep_f from his
# latent state -- so that recovery in the vignette demonstrates the
# estimators rather than luck.
#
# Run from the package root:  source("data-raw/make_ancienregime_extras.R")
# Idempotent under the fixed seed.
# ---------------------------------------------------------------------

set.seed(1789)

provinces <- c("Alsace", "Aquitaine", "Auvergne", "Bretagne", "Champagne",
               "Dauphine", "Gascogne", "Ile-de-France", "Languedoc",
               "Lorraine", "Normandie", "Provence")
occupations <- c("Vagabond", "Metayer", "Journalier", "Petit Metiers",
                 "Petite Bourgeoisie", "Haute Bourgeoisie")
J_occ <- length(occupations)

draw_from <- function(P, i) {
  # one categorical draw per element of i, row i of P
  cp <- t(apply(P, 1L, cumsum))
  as.integer(rowSums(runif(length(i)) > cp[i, , drop = FALSE]) + 1L)
}

# ---------------------------------------------------------------------
# The parish register: 20 parishes per province, with names assembled
# from pools a period clerk would recognise, and populations that are
# very unequal -- one small city per province, a scatter of bourgs, and
# many hamlets. Unequal populations are the point: they make the rival
# law informative and the tabulation sparse.
# ---------------------------------------------------------------------

saints <- c("St-Jean", "St-Martin", "Ste-Colombe", "St-Aubin", "Ste-Foy",
            "St-Loup", "St-Genest", "Ste-Radegonde", "St-Ouen", "St-Fiacre",
            "Ste-Eulalie", "St-Saturnin", "St-Vaast", "Ste-Menehould",
            "St-Gildas", "St-Pardoux", "Ste-Sigolene", "St-Emilion",
            "St-Flour", "Ste-Ame")
suffixes <- c("le-Vieux", "la-Foret", "sur-Riviere", "les-Vignes",
              "en-Plaine", "le-Haut", "les-Bains", "du-Marais",
              "la-Chapelle", "aux-Bois", "le-Comtal", "des-Pres",
              "sous-Roche", "en-Vallee", "la-Ville", "le-Bourg",
              "les-Moulins", "du-Gue", "la-Lande", "sur-Colline")

n_par_per <- 20L
J_par <- length(provinces) * n_par_per

parish_register <- do.call(rbind, lapply(seq_along(provinces), function(p) {
  nm <- paste(saints, suffixes[sample.int(20L)], sep = "-")
  # one city, four bourgs, fifteen hamlets
  w <- c(30, rep(6, 4), rep(1, 15))[sample.int(20L)]
  data.frame(parish = (p - 1L) * n_par_per + seq_len(n_par_per),
             parish_name = nm, province = provinces[p],
             weight_1750 = w, stringsAsFactors = FALSE)
}))
parish_register$share_1750 <-
  parish_register$weight_1750 / sum(parish_register$weight_1750)

# ---------------------------------------------------------------------
# The decade-long parish-to-parish transition operator T_par: stay with
# high probability, move within the province with a preference for the
# city, cross provinces rarely. One matrix, applied per decade.
# ---------------------------------------------------------------------

stay <- 0.82; within <- 0.15; across <- 0.03
T_par <- matrix(0, J_par, J_par)
for (j in seq_len(J_par)) {
  pr <- parish_register$province[j]
  own <- parish_register$parish[parish_register$province == pr]
  oth <- setdiff(seq_len(J_par), own)
  w_own <- parish_register$weight_1750[own]; w_own[own == j] <- 0
  T_par[j, own] <- within * w_own / sum(w_own)
  T_par[j, oth] <- across * parish_register$weight_1750[oth] /
    sum(parish_register$weight_1750[oth])
  T_par[j, j] <- T_par[j, j] + stay
}
T_par <- T_par / rowSums(T_par)

# Latent parish margins at each roll, from the 1750 shares pushed
# through T. These are the "full rolls" a clerk could consult, and the
# per-province rival laws are their within-province renormalisations.
m_1750 <- parish_register$share_1750
m_1770 <- as.numeric(m_1750 %*% T_par %*% T_par)
m_1780 <- as.numeric(m_1770 %*% T_par)

local_slab <- function(margin, prov) {
  idx <- parish_register$province == prov
  out <- numeric(J_par)
  out[idx] <- margin[idx] / sum(margin[idx])
  out
}

# ---------------------------------------------------------------------
# ancienregime_parishes: N men. True 1750 parish observed; links into
# 1770 and 1780 fail with rates alpha1, alpha2, and a failed link draws
# a rival from the BIRTH province's roll at that date. Failures are
# independent across the two links here -- dependence is the lineages
# dataset's subject.
# ---------------------------------------------------------------------

N_par <- 120000L
alpha1_par <- 0.07; alpha2_par <- 0.12

p1750 <- sample.int(J_par, N_par, replace = TRUE, prob = m_1750)
prov_birth <- parish_register$province[p1750]

lat1770 <- draw_from(T_par %*% T_par, p1750)   # two decades
lat1780 <- draw_from(T_par, lat1770)           # one more

slab70 <- lapply(provinces, function(p) local_slab(m_1770, p))
slab80 <- lapply(provinces, function(p) local_slab(m_1780, p))
names(slab70) <- names(slab80) <- provinces

fail1 <- runif(N_par) < alpha1_par
fail2 <- runif(N_par) < alpha2_par
obs1770 <- lat1770
obs1780 <- lat1780
for (p in provinces) {
  i1 <- fail1 & prov_birth == p
  if (any(i1)) obs1770[i1] <- sample.int(J_par, sum(i1), TRUE, prob = slab70[[p]])
  i2 <- fail2 & prov_birth == p
  if (any(i2)) obs1780[i2] <- sample.int(J_par, sum(i2), TRUE, prob = slab80[[p]])
}

ancienregime_parishes <- data.frame(
  province_birth = prov_birth,
  parish_1750 = p1750,
  parish_1770_linked = obs1770,
  parish_1780_linked = obs1780,
  stringsAsFactors = FALSE
)

# The vingtieme of 1780 asked every household where it resided in 1770.
# Aggregated over the FULL roll (not the linked sample), that recall
# supplies the decade operator an analyst can bring to the linked data
# -- the toy counterpart of the 1940 census's five-year question. We
# tabulate it from an independent large sample of the world so it
# carries sampling noise like any auxiliary table.
N_ving <- 400000L
v70 <- sample.int(J_par, N_ving, TRUE, prob = m_1770)
v80 <- draw_from(T_par, v70)
ving <- aggregate(list(n = rep(1L, N_ving)),
                  by = list(j = v70, l = v80), FUN = sum)
ving <- ving[order(ving$j, ving$l), ]
ving$t <- ving$n / ave(ving$n, ving$j, FUN = sum)
vingtieme_flows <- ving[, c("j", "l", "t", "n")]
rownames(vingtieme_flows) <- NULL

# ---------------------------------------------------------------------
# The occupation world for the lineages: transmission Pi (father ->
# son-at-1770), a decade operator T_occ with slight upward drift, a
# common latent slab for rivals, and the clerk's adjacent-rung kernel.
# ---------------------------------------------------------------------

pi_father <- c(0.05, 0.30, 0.28, 0.20, 0.13, 0.04)

Pi_trans <- matrix(0, J_occ, J_occ)  # rows father, cols son-at-1770
for (i in seq_len(J_occ)) {
  d <- abs(seq_len(J_occ) - i)
  w <- 0.55 * (d == 0) + 0.28 * (d == 1) + 0.12 * (d == 2) + 0.05 * (d >= 3)
  w <- w / sum(w)
  up <- pmin(seq_len(J_occ), J_occ) / J_occ      # mild pull upward
  w <- w * (0.85 + 0.3 * up)
  Pi_trans[i, ] <- w / sum(w)
}

T_occ <- matrix(0, J_occ, J_occ)
for (i in seq_len(J_occ)) {
  T_occ[i, i] <- 0.80
  if (i > 1) T_occ[i, i - 1] <- 0.06
  if (i < J_occ) T_occ[i, i + 1] <- 0.14 else T_occ[i, i] <- T_occ[i, i] + 0.14
  if (i == 1) T_occ[i, i] <- T_occ[i, i] + 0.06
}
T_occ <- T_occ / rowSums(T_occ)

# the common rival slab: the 1770 latent occupation margin
occ_slab <- as.numeric((pi_father %*% Pi_trans))
occ_slab <- occ_slab / sum(occ_slab)

# the clerk's kernel: a miscoded occupation lands on an adjacent rung
K_clerk <- matrix(0, J_occ, J_occ)
for (i in seq_len(J_occ)) {
  lo <- max(1L, i - 1L); hi <- min(J_occ, i + 1L)
  nb <- setdiff(lo:hi, i)
  K_clerk[i, nb] <- 1 / length(nb)
}

# ---------------------------------------------------------------------
# ancienregime_lineages: N father-son pairs. Shared-rival link failures
# for the son's three links; a linked (fallible) reading of the father
# in the 1745 dixieme; clerk miscoding on every recorded occupation.
# ---------------------------------------------------------------------

N_lin <- 100000L
pi_s <- 0.06
b_lin <- c(0.045, 0.075, 0.105)
alpha_lin <- pi_s + (1 - pi_s) * b_lin      # 0.102, 0.130, 0.159
s_lin <- pi_s / alpha_lin[1]                # ~ 0.586
alpha_f <- 0.20
mu_clerk <- 0.06

fa <- sample.int(J_occ, N_lin, replace = TRUE, prob = pi_father)
son1 <- draw_from(Pi_trans, fa)             # son's latent, 1770
son2 <- draw_from(T_occ, son1)              # 1780
son3 <- draw_from(T_occ, son2)              # 1790

riv1 <- sample.int(J_occ, N_lin, replace = TRUE, prob = occ_slab)
riv2 <- draw_from(T_occ, riv1)
riv3 <- draw_from(T_occ, riv2)

shared <- runif(N_lin) < pi_s               # every link fails to the rival
f1 <- shared | (runif(N_lin) < b_lin[1])
f2 <- shared | (runif(N_lin) < b_lin[2])
f3 <- shared | (runif(N_lin) < b_lin[3])

fresh <- function(n) sample.int(J_occ, n, replace = TRUE, prob = occ_slab)
attach_one <- function(f, own, riv) {
  out <- own
  fr <- f & !shared                          # independent failure: fresh draw
  out[fr] <- fresh(sum(fr))
  out[f & shared] <- riv[f & shared]         # shared failure: THE rival
  out
}
y1 <- attach_one(f1, son1, riv1)
y2 <- attach_one(f2, son2, riv2)
y3 <- attach_one(f3, son3, riv3)

# the father's own link into the 1745 dixieme: correct = one T step from
# his latent state, failed = a draw from the father-occupation slab
ff <- runif(N_lin) < alpha_f
xf <- draw_from(T_occ, fa)
xf[ff] <- sample.int(J_occ, sum(ff), replace = TRUE, prob = pi_father)

# the clerk miscodes EVERY recorded occupation with probability mu,
# whatever its source -- the truth, a rival, or a rival
miscode <- function(v) {
  hit <- runif(length(v)) < mu_clerk
  v[hit] <- draw_from(K_clerk, v[hit])
  v
}
X_rec  <- miscode(fa)
Xf_rec <- miscode(xf)
Y1_rec <- miscode(y1)
Y2_rec <- miscode(y2)
Y3_rec <- miscode(y3)

ancienregime_lineages <- data.frame(
  province = sample(provinces, N_lin, replace = TRUE,
                    prob = tabulate(match(prov_birth, provinces), 12L)),
  father_occupation_1750 = occupations[X_rec],
  father_occupation_1745_linked = occupations[Xf_rec],
  son_occupation_1770_linked = occupations[Y1_rec],
  son_occupation_1780_linked = occupations[Y2_rec],
  son_occupation_1790_linked = occupations[Y3_rec],
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------
# The truth object, and the write-out
# ---------------------------------------------------------------------

ancienregime_truth <- list(
  occupations = occupations,
  provinces = provinces,
  parish_register = parish_register[, c("parish", "parish_name",
                                        "province", "share_1750")],
  parishes = list(alpha = c(alpha1_par, alpha2_par),
                  T_decade = T_par,
                  margin_1750 = m_1750, margin_1770 = m_1770,
                  margin_1780 = m_1780),
  vingtieme_flows = vingtieme_flows,
  lineages = list(alpha = alpha_lin, s = s_lin, pi_shared = pi_s,
                  b = b_lin, alpha_f = alpha_f, mu = mu_clerk,
                  Pi = Pi_trans, T_occ = T_occ,
                  occupation_slab = occ_slab,
                  father_margin = pi_father, kernel = K_clerk)
)

usethis_save <- function(obj, name) {
  assign(name, obj)
  save(list = name, file = file.path("data", paste0(name, ".rda")),
       compress = "xz")
  message(name, ": ", format(object.size(obj), units = "MB"))
}
usethis_save(ancienregime_parishes, "ancienregime_parishes")
usethis_save(ancienregime_lineages, "ancienregime_lineages")
usethis_save(ancienregime_truth, "ancienregime_truth")
