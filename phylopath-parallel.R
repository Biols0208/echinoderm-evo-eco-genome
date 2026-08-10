
# Purpose : Test causal hypotheses linking effective population size (Ne),
#           transposable element (TE) activity, and genome size evolution
#           using phylogenetic path analysis (phylopath, van der Bijl 2018).
#
# Variable legend:
#   GS      	= log(genome size)                									[response]
#   Total      	= total LTR + DNA content (% genome)
#   LTR     	= LTR content (% genome)
#   DNA     	= DNA transposon content (% genome)
#   GC      	= GC content (%)
#   SLF     	= LTR solo/(intact+solo) ratio           							[solo-LTR fraction]
#   YLTR5  	    = proportion of TE that are low-divergence LTRs (K2P < 5%) 		    [low-divergence LTR activity]
#   YDNA5  	    = proportion of TE that are low-divergence DNA TEs (K2P < 5%) 		[low-divergence DNA TE activity]
#   DNDS    	= dN/dS ratio                     									[proxy for selection efficiency]
#   HET     	= heterozygosity                  									
#   RM      	= reproductive mode (Lecithotrophy, Planktotrophy) 					[binary]
#   DLC     	= log(depth suitability center)   									[habitat depth]
#   TLC     	= temperature suitability center  									[habitat temperature]
#   SLC     	= salinity suitability center     									[habitat salinity]
# =============================================================================

# 0. Global settings
rm(list = ls())

args <- commandArgs(T)
tree_file <- args[1]
data_file <- args[2]

out_dir <- args[3]
if (!dir.exists(out_dir)) dir.create(out_dir)

# Species excluded due to missing data or phylogenetic placement issues
drop_species <- c("Pisaster_ochraceus", "Ophiura_ophiura", "Lytechinus_pictus", "Ptychodera_flava", "Schizocardium_californicum" )

# Collinearity threshold for Pearson |r|
collinearity_cutoff <- 0.7

# CICc cutoff for model averaging
# Models within avg_cutoff CICc units of the model are averaged
avg_cutoff <- 2

## set bootstrap and thread
set.seed(12345)
bootstrap_set <- 5000
parall_works_num <- 40

sink("log.txt", split = TRUE)

# 1. Load packages
library(phylopath)  # Phylogenetic path analysis (van der Bijl 2018)
library(ape)
library(ggplot2)
library(dplyr)
library(tibble)
library(readr)
library(nlme)       # PGLS for residual diagnostics
library(car)
library(phylolm)
library(future.apply)

# 2. Read and prune phylogeny
tree_raw <- read.tree(tree_file)
my_tree <- drop.tip(tree_raw, drop_species)
cat("\nTree tips after pruning:", length(my_tree$tip.label), "\n")

# 3. Read data and check variables
my_data <- read.csv(data_file, stringsAsFactors = FALSE)
cat("\nColumn names in input data:\n")
print(colnames(my_data))
rownames(my_data) <- my_data$Species

required_vars <- c(
  "Species", "log_genomesize", "Total", "LTR", "DNA",
  "Heterozygosity", "dNdS", "Reproductive_Mode",
  "depth_suit_center_log", "temp_suit_center", "sal_suit_center",
  "LTR_solo_intact_ratio", "k2p5_DNA", "k2p5_LTR"
)

missing_vars <- setdiff(required_vars, colnames(my_data))

if (length(missing_vars) > 0) {
  cat("\nERROR: Missing variables in input data:\n")
  print(missing_vars)
  stop("Some required variables are missing. Please check column names.")
}

cat("\nData rows:", nrow(my_data), "\n")
cat("Data columns:", ncol(my_data), "\n")

# 4. Match tree and data
species_in_both <- intersect(my_tree$tip.label, my_data$Species)
cat("\nSpecies in tree:", length(my_tree$tip.label), "\n")
cat("Species in data:", nrow(my_data), "\n")
cat("Species in both:", length(species_in_both), "\n")

my_tree <- keep.tip(my_tree, species_in_both)
my_data <- my_data[species_in_both, ]
my_data <- my_data[my_tree$tip.label, ]
stopifnot(all(rownames(my_data) == my_tree$tip.label))

# 5. Build analysis data frame
vars_for_complete <- required_vars[required_vars != "Species"]
complete_idx <- complete.cases(my_data[, vars_for_complete])
cat(sprintf(
  "\nComplete cases: %d / %d\n",
  sum(complete_idx), nrow(my_data)
))

my_data_clean <- my_data[complete_idx, ]
my_tree_clean <- keep.tip(my_tree, my_data_clean$Species)
my_data_clean <- my_data_clean[my_tree_clean$tip.label, ]

head(my_data_clean, n = 2)

# Construct the working data frame with short variable names
data_full <- data.frame(
  GS     = my_data_clean$log_genomesize,
  TE     = my_data_clean$Total,
  LTR    = my_data_clean$LTR,
  DNA    = my_data_clean$DNA,
  GC     = my_data_clean$GC,

  HET    = my_data_clean$Heterozygosity, 
  DNDS   = my_data_clean$dNdS,				

  DLC    = my_data_clean$depth_suit_center_log,
  TLC    = my_data_clean$temp_suit_center,
  SLC    = my_data_clean$sal_suit_center,

  SLF    = my_data_clean$LTR_solo_intact_ratio,

  YDNA5  = my_data_clean$k2p5_DNA,
  YLTR5  = my_data_clean$k2p5_LTR,

  RM     = my_data_clean$Reproductive_Mode
)

rownames(data_full) <- my_data_clean$Species
stopifnot(all(rownames(data_full) == my_tree_clean$tip.label))

# 6. Scale continuous variables
log_positive <- function(x, var_name) {
  ok <- !is.na(x)

  if (any(x[ok] <= 0)) {
    stop(
      var_name,
      " contains zero or negative values; do not add an arbitrary constant ",
      "without first examining these observations."
    )
  }

  log(x)
}

logit_open01 <- function(x, var_name) {
  ok <- !is.na(x)

  if (any(x[ok] < 0 | x[ok] > 1)) {
    stop(var_name, " contains values outside [0, 1].")
  }

  if (any(x[ok] == 0 | x[ok] == 1)) {
    stop(
      var_name,
      " contains exact 0 or 1 values and cannot be directly logit-transformed."
    )
  }

  qlogis(x)
}

adjust_prop <- function(x) {
  n <- sum(is.finite(x))
  (x * (n - 1) + 0.5) / n
}

data_full_scaled <- data_full

data_full_scaled$DNDS <- log_positive(data_full_scaled$DNDS, "DNDS")
data_full_scaled$HET <- log_positive(data_full_scaled$HET, "HET")

te_percent_cols <- c("TE", "LTR", "DNA")
for (v in te_percent_cols) {
  data_full_scaled[[v]] <- adjust_prop(
    data_full_scaled[[v]] / 100
  )
}

ratio_cols <- c("SLF", "YDNA5", "YLTR5")

for (v in ratio_cols) {
  data_full_scaled[[v]] <- adjust_prop(
    data_full_scaled[[v]]
  )
}

data_full_scaled$RM <- as.character(data_full_scaled$RM)
data_full_scaled$RM <- ifelse(data_full_scaled$RM %in% c("1", "RM1"),"RM1","RM0")
data_full_scaled$RM <- factor(data_full_scaled$RM, levels = c("RM0", "RM1"))
data_full_scaled$RM_raw <- data_full_scaled$RM

scale_cols <- names(data_full_scaled)[vapply(data_full_scaled, is.numeric, logical(1))]

data_full_scaled[scale_cols] <- lapply(
  data_full_scaled[scale_cols],
  function(x) as.numeric(scale(x))
)

direct_vars <- c("DNDS", "HET", "DLC", "TLC", "SLC", "RM", "GC")
fits_GS_direct <- lapply(direct_vars, function(x) {
  vars <- c("GS", "TE", x)
  dat_x <- data_full_scaled[
    complete.cases(data_full_scaled[, vars]),
    ,
    drop = FALSE
  ]

  tree_x <- drop.tip(
    my_tree_clean,
    setdiff(my_tree_clean$tip.label, rownames(dat_x))
  )

  dat_x <- dat_x[tree_x$tip.label, , drop = FALSE]

  form <- reformulate(
    termlabels = c("TE", x),
    response = "GS"
  )

  phylolm(
    formula = form,
    data = dat_x,
    phy = tree_x,
    model = "lambda"
  )
})

names(fits_GS_direct) <- direct_vars
coef_tables <- lapply(fits_GS_direct, function(m) { summary(m)$coefficients})
coef_tables

model_compare <- data.frame(
  variable = names(fits_GS_direct),
  logLik = sapply(fits_GS_direct, logLik),
  lambda = sapply(fits_GS_direct, function(m) m$optpar),
  TE_coef = sapply(fits_GS_direct, function(m) coef(m)["TE"]),
  direct_coef = sapply(seq_along(fits_GS_direct), function(i) {
    coefs <- coef(fits_GS_direct[[i]])
    coefs[setdiff(names(coefs), c("(Intercept)", "TE"))][1]
  })
)
model_compare

TE_vars <- c("DNDS", "HET", "DLC", "TLC", "SLC", "RM", "GC")
fits_TE_upstream <- lapply(TE_vars, function(x) {
  vars <- c("TE", x)
  dat_x <- data_full_scaled[
    complete.cases(data_full_scaled[, vars]),
    ,
    drop = FALSE
  ]

  tree_x <- drop.tip(
    my_tree_clean,
    setdiff(my_tree_clean$tip.label, rownames(dat_x))
  )

  dat_x <- dat_x[tree_x$tip.label, , drop = FALSE]

  form <- reformulate(
    termlabels = x,
    response = "TE"
  )

  phylolm(
    formula = form,
    data = dat_x,
    phy = tree_x,
    model = "lambda"
  )
})

names(fits_TE_upstream) <- TE_vars
coef_tables <- lapply(fits_TE_upstream, function(m) { summary(m)$coefficients})
coef_tables

model_compare <- data.frame(
  variable = names(fits_TE_upstream),
  logLik = sapply(fits_TE_upstream, logLik),
  lambda = sapply(fits_TE_upstream, function(m) m$optpar),
  direct_coef = sapply(seq_along(fits_TE_upstream), function(i) {
    coefs <- coef(fits_TE_upstream[[i]])
    coefs[setdiff(names(coefs),"(Intercept)")][1]
  })
)
model_compare

DNDS_vars <- c("DLC", "TLC", "SLC", "RM", "GC")
fits_DNDS <- lapply(DNDS_vars, function(x) {
  vars <- c("DNDS", x)
  dat_x <- data_full_scaled[
    complete.cases(data_full_scaled[, vars]),
    ,
    drop = FALSE
  ]

  tree_x <- drop.tip(
    my_tree_clean,
    setdiff(my_tree_clean$tip.label, rownames(dat_x))
  )

  dat_x <- dat_x[tree_x$tip.label, , drop = FALSE]

  form <- reformulate(
    termlabels = x,
    response = "DNDS"
  )

  phylolm(
    formula = form,
    data = dat_x,
    phy = tree_x,
    model = "lambda"
  )
})

names(fits_DNDS) <- DNDS_vars
coef_tables <- lapply(fits_DNDS, function(m) { summary(m)$coefficients})
coef_tables

model_compare <- data.frame(
  variable = names(fits_DNDS),
  logLik = sapply(fits_DNDS, logLik),
  lambda = sapply(fits_DNDS, function(m) m$optpar),
  direct_coef = sapply(seq_along(fits_DNDS), function(i) {
    coefs <- coef(fits_DNDS[[i]])
    coefs[setdiff(names(coefs), "(Intercept)")][1]
  })
)
model_compare

YDNA5_vars <- c("DLC", "TLC", "SLC", "RM", "GC")
fits_YDNA5 <- lapply(YDNA5_vars, function(x) {
  vars <- c("YDNA5", x)
  dat_x <- data_full_scaled[
    complete.cases(data_full_scaled[, vars]),
    ,
    drop = FALSE
  ]

  tree_x <- drop.tip(
    my_tree_clean,
    setdiff(my_tree_clean$tip.label, rownames(dat_x))
  )

  dat_x <- dat_x[tree_x$tip.label, , drop = FALSE]

  form <- reformulate(
    termlabels = x,
    response = "YDNA5"
  )

  phylolm(
    formula = form,
    data = dat_x,
    phy = tree_x,
    model = "lambda"
  )
})

names(fits_YDNA5) <- YDNA5_vars
coef_tables <- lapply(fits_YDNA5, function(m) { summary(m)$coefficients})
coef_tables

model_compare <- data.frame(
  variable = names(fits_YDNA5),
  logLik = sapply(fits_YDNA5, logLik),
  lambda = sapply(fits_YDNA5, function(m) m$optpar),
  direct_coef = sapply(seq_along(fits_YDNA5), function(i) {
    coefs <- coef(fits_YDNA5[[i]])
    coefs[setdiff(names(coefs), "(Intercept)")][1]
  })
)
model_compare

##  --- Hypothesis and Models --- ##
## Does the positive correlation between DNDS and genome size mediate through total TE content? Or are environmental variables better able to explain TE accumulation than DNDS?
models_RM_upstream <- define_model_set(
  A0_RM = c(GS ~ TE + RM, TE ~ TLC),
  A1_DLC_RM = c(GS ~ TE + RM, TE ~ TLC, RM ~ DLC),
  A2_TLC_RM = c(GS ~ TE + RM, TE ~ TLC, RM ~ TLC),
  A3_D_TLC_RM = c(GS ~ TE + RM, TE ~ TLC, RM ~ DLC + TLC),
  .common = c(TLC ~ DLC))

models_env_GS <- define_model_set(
  B0_bb = c(GS ~ TE + RM, TE ~ TLC),
  B1_DLC_dir = c(GS ~ TE + RM + DLC, TE ~ TLC),
  B2_TLC_dir = c(GS ~ TE + RM + TLC, TE ~ TLC),
  B3_D_TLC_dir = c(GS ~ TE + RM + DLC + TLC, TE ~ TLC),
  .common = c(TLC ~ DLC, RM ~ TLC))

models_env_TE <- define_model_set(
  # TE-burden baseline. Total TE content predicts GS, but no upstream correlate of TE is specified.
  C0_TE = c(GS ~ TE + RM + DLC + TLC),
  # Thermal association. TLC is associated with TE abundance.
  C1_TLC_TE = c(GS ~ TE + RM + DLC + TLC, TE ~ TLC),
  # Direct bathymetric association. DLC is associated with TE abundance.
  C2_DLC_TE = c(GS ~ TE + RM + DLC + TLC, TE ~ DLC),
  # Reproductive-mode association. RM is associated with both TE abundance and GS.
  C3_RM_TE = c(GS ~ TE + RM + DLC + TLC, TE ~ RM),
  # Thermal association with residual depth effect. DLC retains information about TE after TLC is included.
  C4_D_TLC_TE = c(GS ~ TE + RM + DLC + TLC, TE ~ TLC + DLC),
  # Thermal and reproductive-mode associations. RM retains information about TE after TLC is included.
  C5_TLC_RM_TE = c(GS ~ TE + RM + DLC + TLC, TE ~ TLC + RM),
  .common = c(TLC ~ DLC, RM ~ TLC))
# .common = c(TLC ~ DLC))

#             Estimate     StdErr      t.value    p.value
#DLC          0.4565210    0.1108911   4.1168418  0.0001377336
#TLC         -0.59171323   0.1075084  -5.503879   1.155248e-06
#RMRM1        1.2549626    0.2178410   5.760910   4.581004e-07
models_DNDS <- define_model_set(
  # Environment/life-history backbone
  D0_bb = c(GS ~ TE + RM + DLC + TLC, TE ~ TLC),
  # Direct DNDS-GS association
  D1_GS_dir = c(GS ~ TE + RM + DLC + TLC + DNDS, TE ~ TLC),
  # TE-mediated DNDS association
  D2_TE_med = c(GS ~ TE + RM + DLC + TLC, TE ~ TLC + DNDS),
  # Partial mediation
  D3_part_med = c(GS ~ TE + RM + DLC + TLC + DNDS, TE ~ TLC + DNDS),
  .common = c(TLC ~ DLC, RM ~ TLC, DNDS ~ DLC + RM + TLC))
# .common = c(TLC ~ DLC, RM ~ TLC))

models_HET <- define_model_set(
  # Environment/life-history backbone
  E0_bb = c(GS ~ TE + RM + DLC + TLC, TE ~ TLC),
  # Direct HET association
  E1_GS_dir = c(GS ~ TE + RM + DLC + TLC + HET, TE ~ TLC),
  # TE-mediated HET association
  E2_TE_med = c(GS ~ TE + RM + DLC + TLC, TE ~ TLC + HET),
  # Partial mediation
  E3_part_med = c(GS ~ TE + RM + DLC + TLC + HET, TE ~ TLC + HET),
  .common = c(TLC ~ DLC, RM ~ TLC, HET ~ RM))
# .common = c(TLC ~ DLC, RM ~ TLC))

# Is the overall pattern primarily derived from LTRs, DNA transposons? Do environmental factors influence different TE types?

models_DNA_LTR <- define_model_set(
  # Turnover baseline
  F0_turnover = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF, DNA ~ YDNA5),
  # TLC associated only with LTR
  F1_TLC_LTR = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5),
  # DLC associated only with LTR
  F2_DLC_LTR = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + DLC, DNA ~ YDNA5),
  # TLC associated only with DNA transposons
  F3_TLC_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF, DNA ~ YDNA5 + TLC),
  # DLC associated only with DNA transposons
  F4_DLC_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF, DNA ~ YDNA5 + DLC),
  # Matched structure suggested by current results
  F5_T_LTR_D_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC),
  # Crossed alternative
  F6_D_LTR_T_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + DLC, DNA ~ YDNA5 + TLC),
  # TLC associated with both TE classes
  F7_TLC_both = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + TLC),
  # DLC associated with both TE classes
  F8_DLC_both = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + DLC, DNA ~ YDNA5 + DLC),
  # Full environmental model
  F9_full = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC + DLC, DNA ~ YDNA5 + TLC + DLC),
  .common = c(TLC ~ DLC, RM ~ TLC, YDNA5 ~ RM))

models_DNA_LTR_DNDS <- define_model_set(
  # Supported environmental class backbone
  G0_bb = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC),
  # DNDS adds to LTR abundance
  G1_LTR = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC + DNDS, DNA ~ YDNA5 + DLC),
  # DNDS adds to DNA-transposon abundance
  G2_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC + DNDS),
  # DNDS adds to both TE classes
  G3_both = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC + DNDS, DNA ~ YDNA5 + DLC + DNDS),
  .common = c(TLC  ~ DLC, RM ~ TLC, DNDS ~ DLC + RM + TLC, YDNA5 ~ RM))
# .common = c(TLC  ~ DLC, RM ~ TLC))

models_DNA_LTR_HET <- define_model_set(
  # Supported environmental class backbone
  H0_bb = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC),
  # DNDS adds to LTR abundance
  H1_LTR = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC + HET, DNA ~ YDNA5 + DLC),
  # DNDS adds to DNA-transposon abundance
  H2_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC + HET),
  # DNDS adds to both TE classes
  H3_both = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC + HET, DNA ~ YDNA5 + DLC + HET),
  .common = c(TLC ~ DLC, RM ~ TLC))

# check YDNA5
models_young_DNA_env <- define_model_set(
  I0_none = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC),
  I1_RM = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC, YDNA5 ~ RM),
  I2_TLC = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC, YDNA5 ~ TLC),
  I3_DLC = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC, YDNA5 ~ DLC),
  I4_TLC_RM = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC, YDNA5 ~ TLC + RM),
  I5_DLC_RM = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC, YDNA5 ~ DLC + RM),
  I6_DLC_TLC = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC, YDNA5 ~ DLC + TLC),
  I7_full = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ YDNA5 + DLC, YDNA5 ~ DLC + TLC + RM),
  .common = c(TLC ~ DLC, RM ~ TLC))

# remove YDNA5
models_DNA_LTR_remove_YDNA5 <- define_model_set(
  J0_turnover = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF),
  J1_TLC_LTR = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC),
  J2_DLC_LTR = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + DLC),
  J3_TLC_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF, DNA ~ TLC),
  J4_DLC_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF, DNA ~ DLC),
  J5_T_LTR_D_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ DLC),
  J6_D_LTR_T_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + DLC, DNA ~ TLC),
  J7_TLC_both = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ TLC),
  J8_DLC_both = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + DLC, DNA ~ DLC),
  J9_full = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC + DLC, DNA ~ TLC + DLC),
  .common = c(TLC ~ DLC, RM ~ TLC))

models_DNA_LTR_DNDS_remove_YDNA5 <- define_model_set(
  K0_bb = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ DLC),
  K1_LTR = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC + DNDS, DNA ~ DLC),
  K2_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ DLC + DNDS),
  K3_both = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC + DNDS, DNA ~ DLC + DNDS),
  .common = c(TLC  ~ DLC, RM ~ TLC, DNDS ~ DLC + RM + TLC))

models_DNA_LTR_HET_remove_YDNA5 <- define_model_set(
  L0_bb = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ DLC),
  L1_LTR = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC + HET, DNA ~ DLC),
  L2_DNA = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC, DNA ~ DLC + HET),
  L3_both = c(GS ~ LTR + DNA + RM + DLC + TLC, LTR ~ YLTR5 + SLF + TLC + HET, DNA ~ DLC + HET),
  .common = c(TLC ~ DLC, RM ~ TLC, HET ~ RM))

# 10. Run phylogenetic path analyses
model_list <- list(
  RM_upstream = models_RM_upstream,
  env_GS = models_env_GS,
  env_TE = models_env_TE,
  DNDS = models_DNDS,
  HET = models_HET,
  DNA_LTR = models_DNA_LTR,
  DNA_LTR_DNDS = models_DNA_LTR_DNDS,
  DNA_LTR_HET = models_DNA_LTR_HET,
  young_DNA_env = models_young_DNA_env
)

results <- list()

for (name in names(model_list)) {

  cat(sprintf("\n========== Running models_%s ==========\n", name))

  results[[name]] <- phylo_path(
    model_list[[name]],
    data = data_full_scaled,
    tree = my_tree_clean,
    model = "lambda", 			# "lambda", "BM"
    method = "logistic_MPLE"	# "logistic_MPLE", "logistic_IG10", "logistic_MLE", "poisson_GEE"
  )
}

order_class <- c("DLC", "TLC", "RM", "SLF", "YLTR5", "LTR", "DNA", "GS")
results[["DNA_LTR_remove_YDNA5"]] <- phylo_path(models_DNA_LTR_remove_YDNA5, data = data_full_scaled, tree = my_tree_clean, model = "lambda", method = "logistic_MPLE", order = order_class)
order_class <- c("DLC", "TLC", "RM", "DNDS", "SLF", "YLTR5", "LTR", "DNA", "GS")
results[["DNA_LTR_DNDS_remove_YDNA5"]] <- phylo_path(models_DNA_LTR_DNDS_remove_YDNA5, data = data_full_scaled, tree = my_tree_clean, model = "lambda", method = "logistic_MPLE", order = order_class)
order_class <- c("DLC", "TLC", "RM", "HET", "SLF", "YLTR5", "LTR", "DNA", "GS")
results[["DNA_LTR_HET_remove_YDNA5"]] <- phylo_path(models_DNA_LTR_HET_remove_YDNA5, data = data_full_scaled, tree = my_tree_clean, model = "lambda", method = "logistic_MPLE", order = order_class)

# 11. D-sep diagnostics
# Purpose: Identify which conditional independence claims are violated by data.
# Interpretation:
#   - Significant tests (p < 0.05) indicate missing paths in the model.
#   - TLC ~ DLC will always be significant (physical collinearity, r = -0.747).
#   - Fewer significant tests = better model structure.
#   - p values are diagnostic tools, NOT pass/fail criteria for model acceptance.
#     Model selection is based on CICc.

sink(file.path(out_dir, "D-sep.diagnostics.txt"))
diagnose_dsep <- function(result, label) {
  cat(sprintf("\n===== D-sep Diagnostics: %s =====\n", label))
  all_sig <- list()
  for (m in names(result$d_sep)) {
    dsep_tbl <- result$d_sep[[m]]
    sig_mask <- !is.na(dsep_tbl$p) & dsep_tbl$p < 0.05
    n_sig <- sum(sig_mask)

    if (n_sig > 0) {
      cat(sprintf(
        "  [%s] %d significant d-sep tests p < 0.05:\n",
        m, n_sig
      ))

      sig_tbl <- dsep_tbl[sig_mask, c("d_sep", "p")]

      for (i in seq_len(nrow(sig_tbl))) {
        cat(sprintf(
          "    %-55s p = %.4f\n",
          sig_tbl$d_sep[i],
          sig_tbl$p[i]
        ))
      }

      all_sig[[m]] <- sig_tbl
    } else {
      cat(sprintf("  [%s] All d-sep p >= 0.05\n", m))
    }
  }

  invisible(all_sig)
}

dsep_results <- Map(
  function(x, n) diagnose_dsep(x, paste0("models_", n)),
  results,
  names(results)
)

sink()

# 12. Model comparison summary
# Output columns:
#   k    = number of d-sep tests (degrees of freedom)
#   q    = number of estimated parameters
#   C    = Fisher's C statistic (goodness of fit; lower = better fit)
#   p    = p-value for C statistic (p > 0.05 = model not rejected)
#   CICc = corrected C-based information criterion (primary selection criterion)
#   delta_CICc = CICc difference from best model
#   w    = Akaike weight (model probability; sums to 1 across model set)
#
# Selection rule: lowest CICc = best balance of fit and parsimony.
# Models within delta_CICc < avg_cutoff are candidates for model averaging.

sink(file.path(out_dir, "model.summary.txt"))
for (name in names(results)) {
  cat(sprintf("\n=== %s ===\n", name))
  cat(sprintf("Model averaging cutoff: delta_CICc < %d\n", avg_cutoff))

  s <- summary(results[[name]])
  print(s, n = Inf)
}
sink()

# 13. get coefficients #
# fitted_DAG objects contain coefficient matrices.
# Convert them to long-format tables with one row per path for readability.
# Columns: from (response) → to (predictor), estimate, SE, lower CI, upper CI
extract_path_table <- function(fitted_dag, model_name = "") {

  coef_mat  <- fitted_dag$coef
  se_mat    <- fitted_dag$se
  lower_mat <- fitted_dag$lower
  upper_mat <- fitted_dag$upper

  path_idx <- which(coef_mat != 0, arr.ind = TRUE)

  if (nrow(path_idx) == 0) {
    warning(sprintf("[%s] No non-zero path coefficients found.", model_name))
    return(data.frame())
  }

  path_table <- data.frame(
    model    = model_name,
    from     = rownames(coef_mat)[path_idx[, 1]], 
    to       = colnames(coef_mat)[path_idx[, 2]],
    estimate = coef_mat[path_idx],
    se       = if (!is.null(se_mat))    se_mat[path_idx]    else NA_real_,
    lower    = if (!is.null(lower_mat)) lower_mat[path_idx] else NA_real_,
    upper    = if (!is.null(upper_mat)) upper_mat[path_idx] else NA_real_,
    stringsAsFactors = FALSE
  )

  # Sort by response variable, then predictor
  path_table <- path_table[order(path_table$to, path_table$from), ]
  rownames(path_table) <- NULL

  return(path_table)
}

if (0) {
cat("Get Model-best coefficients")
best_models <- list()
tbl_best <- list()

plan(multisession, workers = parall_works_num)
res_list <- future_lapply(names(results), function(name) {
  tryCatch({
    model <- best(results[[name]], boot = bootstrap_set)
    tbl <- extract_path_table(
      model,
      model_name = paste0("best_", name)
    )
    
    outfile <- file.path(
      out_dir,
      paste0("coef_best_model_", name, ".csv")
    )
    
    write.csv(tbl, outfile, row.names = FALSE)
    
    list(
      name = name,
      model = model,
      table = tbl
    )
    
  }, error = function(e) {
    
    list(
      name = name,
      model = NULL,
      table = NULL
    )
  })
})

names(res_list) <- names(results)

best_models <- lapply(res_list, `[[`, "model")
tbl_best <- lapply(res_list, `[[`, "table")
}

cat("Get Model-averaged coefficients")
avg_models <- list()
tbl_avg <- list()

plan(multisession, workers = parall_works_num)
res_list <- future_lapply(names(results), function(name) {
  tryCatch({
    model <- average(results[[name]], 
					 cut_off = avg_cutoff,
                     avg_method = "conditional",		
	                 boot = bootstrap_set)
    tbl <- extract_path_table(
      model,
      model_name = paste0("averaged_", name)
    )

    outfile <- file.path(
      out_dir,
      paste0("coef_averaged_", name, ".csv")
    )

    write.csv(tbl, outfile, row.names = FALSE)

    list(
      name = name,
      model = model,
      table = tbl
    )

  }, error = function(e) {

    list(
      name = name,
      model = NULL,
      table = NULL
    )
  })
})

names(res_list) <- names(results)

avg_models <- lapply(res_list, `[[`, "model")
tbl_avg <- lapply(res_list, `[[`, "table")


# 14. Visualization
# Safer plotting helpers to avoid ggraph arc errors
safe_print_plot <- function(p) {
  tryCatch({
    if (inherits(p, "gg") || inherits(p, "ggplot")) {
      print(p)
    } else {
      # Many phylopath plot() calls draw directly; try printing anyway
      try(print(p), silent = TRUE)
    }
  }, error = function(e) {
    message("Plot skipped due to error: ", e$message)
    plot.new(); title(main = "Plot skipped due to plotting error")
  })
}

pdf( file.path(out_dir, "phylo.pdf"), width = 12, height = 8 )

positions_total <- data.frame(
  name = c(
          "TLC", 
    "DLC", "RM",  "TE", "GS"
  ),
  x = c(
          1.9,
    1.0,  1.7,  3, 5
 ),
  y = c(
          1.52,
    1.5,  1.495,  1.506, 1.5
  )
)

safe_print_plot( try(plot(avg_models[["env_TE"]], 
  manual_layout = positions_total,
  type          = "width",
  colors        = c("#3B4CC0", "#B40426"),
  text_size     = 4.2,
  box_x         = 18,
  box_y         = 11,
  edge_width    = 1.1,
  curvature     = 0.001,
  show.legend   = TRUE
), silent = TRUE) )

positions_class <- data.frame(
  name = c(
    "DLC", "TLC", "RM",
    "YLTR5", "SLF",
    "LTR", "DNA",
    "GS"
  ),
  x = c(
    1.0, 2.2, 2.0,
    1.2, 1.2, 
    3.3, 3.3, 
    4.7
  ),
  y = c(
    5, 6.0, 4.8,
    5.9, 6.1, 
    6.3, 4.5, 
    5.1
  )
)

safe_print_plot( try(plot(avg_models[["DNA_LTR_remove_YDNA5"]], 
  manual_layout = positions_class,
  type          = "width",
  colors        = c("#3B4CC0", "#B40426"),
  text_size     = 3.8,
  box_x         = 19,
  box_y         = 11,
  edge_width    = 1.0,
  curvature     = 0.08,
  show.legend   = TRUE
), silent = TRUE))

positions_class <- data.frame(
  name = c(
    "DLC", "TLC", "RM",
    "YLTR5", "SLF", "YDNA5",
    "LTR", "DNA",
    "GS"
  ),
  x = c(
    1.0, 2.2, 2.0,
    1.2, 1.2, 2.6,
    3.3, 3.3,
    4.7
  ),
  y = c(
    5, 6.0, 4.8,
    5.9, 6.1, 4.6,
    6.3, 4.5,
    5.1
  )
)

safe_print_plot( try(plot(avg_models[["DNA_LTR"]],
  manual_layout = positions_class,
  type          = "width",
  colors        = c("#3B4CC0", "#B40426"),
  text_size     = 3.8,
  box_x         = 19,
  box_y         = 11,
  edge_width    = 1.0,
  curvature     = 0.08,
  show.legend   = TRUE
), silent = TRUE))

# Coefficient plot (ranked by effect size)
if (0) {
print(
  coef_plot(best_models[["env_TE"]], order_by = "strength") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = sprintf(
      "Best model coefficients: env_TE (delta_CICc < %d)", avg_cutoff
    ))
)

print(
  coef_plot(best_models[["DNA_LTR_remove_YDNA5"]], order_by = "strength") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = sprintf(
      "Best model coefficients: DNA_LTR_remove_YDNA5 (delta_CICc < %d)", avg_cutoff
    ))
)

print(
  coef_plot(best_models[["DNA_LTR"]], order_by = "strength") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = sprintf(
      "Best model coefficients: DNA_LTR (delta_CICc < %d)", avg_cutoff
    ))
)
}

print(
  coef_plot(avg_models[["env_TE"]], order_by = "strength") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = sprintf(
      "Averaged model coefficients: env_TE (delta_CICc < %d)", avg_cutoff
    ))
)

print(
  coef_plot(avg_models[["DNA_LTR_remove_YDNA5"]], order_by = "strength") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = sprintf(
      "Averaged model coefficients: DNA_LTR_remove_YDNA5 (delta_CICc < %d)", avg_cutoff
    ))
)

print(
  coef_plot(avg_models[["DNA_LTR"]], order_by = "strength") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = sprintf(
      "Averaged model coefficients: DNA_LTR (delta_CICc < %d)", avg_cutoff
    ))
)

dev.off()

cat("Analysis complete.\n")
