#!/usr/bin/env Rscript
# =============================================================================
# Empirical Bayes disproportionality analysis of FAERS drug-drug-event pairs
# =============================================================================
# Purpose:
#   Estimate empirical Bayes geometric means (EBGM) and posterior credibility
#   bounds (EB05 and EB95) using the openEBGM R package.
#
# Input: CSV with columns caseid, product, outcome. Run the script separately
#   for (1) parent/composite outcomes (MACE, OCE) and (2) individual component
#   outcomes. Do not combine the two analysis levels in one input file.
#   Each dataset compares one ARPI–cardiovascular medication combination
#   against the cardiovascular medication alone. The product column includes
#   both exposure labels. Each row represents one case–exposure–outcome. Input data
#   must already be cleaned, deduplicated, and prepared for the analysis.
#
# Usage:
#   Rscript ebgm_analysis.R input/analysis_data.csv output/ebgm_results.csv
#
# Dependencies: openEBGM (install.packages("openEBGM"))
#
# IMPORTANT:
#   This script computes EBGM, EB05, EB95 and observed/expected ratio (N/E).
#   N/E is NOT a reporting odds ratio (ROR). It does not estimate causal risk.
#   EB05/EB95 are the 5th/95th posterior percentiles (90% interval).
# =============================================================================

suppressPackageStartupMessages(library(openEBGM))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript ebgm_analysis.R <input.csv> <output.csv>", call. = FALSE)
}
input_file <- args[[1L]]
output_file <- args[[2L]]
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

# 1. Import and validate required fields -------------------------------------
raw <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("caseid", "product", "outcome")
missing_columns <- setdiff(required, names(raw))
if (length(missing_columns)) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

# Retain only columns needed for the EBGM analysis. This avoids dropping
# otherwise eligible rows because of missing values in unrelated variables.
dat <- raw[, required, drop = FALSE]
names(dat) <- c("id", "var1", "var2")
for (nm in names(dat)) {
  dat[[nm]] <- trimws(as.character(dat[[nm]]))
  dat[[nm]][dat[[nm]] == ""] <- NA_character_
}
# Match original approach: blank outcome values become "none".
dat$var2[is.na(dat$var2)] <- "none"
dat <- dat[complete.cases(dat), , drop = FALSE]
if (nrow(dat) == 0L) stop("No eligible records after cleaning.")

# A case–exposure–event triple should occur at most once. Report duplicates
# rather than silently allowing them to inflate observed counts.
if (anyDuplicated(dat[c("id", "var1", "var2")])) {
  stop("Duplicate case-product-outcome rows detected. Deduplicate upstream before analysis.")
}
cat("Unique cases:", length(unique(dat$id)), "\n")
cat("Case–exposure–event records:", nrow(dat), "\n")

# 2. Generate observed and expected counts -----------------------------------
processed <- processRaw(dat)
stopifnot(all(is.finite(processed$N)), all(is.finite(processed$E)),
          all(processed$N >= 0), all(processed$E >= 0))
if (any(processed$E <= 0)) {
  stop("Expected counts must be positive for posterior estimation.")
}

# 3. Squash counts for computational efficiency ------------------------------
squashed <- autoSquash(processed)

# 4. Negative log-likelihood on an unconstrained parameter scale -------------
# Gamma shape/rate parameters are positive (exp transform); mixture
# probability p lies between 0 and 1 (logistic transform).
safe_negLL <- function(z, ni, ei, wi) {
  if (any(!is.finite(z))) return(1e100)
  pars <- c(pmax(exp(z[1:4]), 1e-8),
            pmin(pmax(plogis(z[5]), 1e-8), 1 - 1e-8))
  result <- tryCatch(
    suppressWarnings(negLLsquash(theta = pars, ni = ni, ei = ei,
                                 wi = wi, N_star = 1)),
    error = function(e) NA_real_
  )
  if (length(result) != 1L || !is.finite(result)) 1e100 else result
}

# 5. Optimize from multiple starting points ----------------------------------
starts <- list(
  c(log(0.2), log(0.1), log(2),   log(4), qlogis(1/3)),
  c(log(0.5), log(0.5), log(2),   log(2), qlogis(0.5)),
  c(log(1),   log(1),   log(3),   log(3), qlogis(0.25)),
  c(log(0.1), log(0.2), log(1.5), log(3), qlogis(0.75))
)
fits <- lapply(starts, function(st) {
  optim(par = st, fn = safe_negLL,
        ni = squashed$N, ei = squashed$E, wi = squashed$weight,
        method = "L-BFGS-B",
        lower = rep(-10, 5), upper = rep(10, 5),
        control = list(maxit = 5000, factr = 1e7, pgtol = 1e-8))
})
fit_table <- data.frame(
  start = seq_along(fits),
  objective = vapply(fits, function(x) x$value, numeric(1)),
  convergence = vapply(fits, function(x) x$convergence, integer(1))
)
print(fit_table)
valid <- which(is.finite(fit_table$objective) &
               fit_table$convergence == 0 &
               fit_table$objective < 1e99)
if (!length(valid)) stop("No optimizer run converged successfully.")
best <- fits[[valid[which.min(fit_table$objective[valid])]]]

# 6. Back-transform and validate mixture parameters --------------------------
z <- best$par
theta_hat <- c(alpha1 = exp(z[1]), beta1 = exp(z[2]),
               alpha2 = exp(z[3]), beta2 = exp(z[4]), p = plogis(z[5]))
if (any(!is.finite(theta_hat)) || any(theta_hat <= 0) ||
    theta_hat["p"] >= 1) stop("Invalid fitted mixture parameters.")
print(theta_hat)

# 7. Calculate posterior weights and signal estimates ------------------------
processed$Qn <- Qn(theta_hat = theta_hat, N = processed$N, E = processed$E)
processed$EBGM <- ebgm(theta_hat = theta_hat, N = processed$N,
                       E = processed$E, qn = processed$Qn)
processed$EB05 <- quantBisect(percent = 5, theta_hat = theta_hat,
                              N = processed$N, E = processed$E,
                              qn = processed$Qn)
processed$EB95 <- quantBisect(percent = 95, theta_hat = theta_hat,
                              N = processed$N, E = processed$E,
                              qn = processed$Qn)
processed$OE_ratio <- processed$N / processed$E

# 8. Export all exposure–event results --------------------------------------------
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write.csv(processed, output_file, row.names = FALSE)
cat("Results saved to:", output_file, "\n")

# 9. Optional: filter either exposure group and outcomes -----------------------
# Two separate analysis runs are used, with identical column names:
#   Parent analysis: outcome contains only composite labels.
#   Component analysis: outcome contains only individual component labels.
# In the parent-level input, outcome contains:
#   "Major adverse cardiovascular event" (MACE)
#   "Other cardiovascular event" (OCE)
# In the component-level input, MACE components are:
#   "Cardiovascular death", "Nonfatal myocardial infarction", "Nonfatal stroke"
# In the component-level input, OCE components are:
#   "General cardiac disorders", "Bleeding", "Thrombosis"
# IMPORTANT: Parent and component analyses are run separately. The script
# does not derive composites from components or pool the two datasets.
# Construct parent outcomes upstream using documented case-level rules.
# The same case may appear in both runs; do not sum their case counts.
# Uncomment and supply a combination or comparator label in product.
# target_product <- "<exposure_label>"
# target_outcomes <- c("Major adverse cardiovascular event",
#                      "Other cardiovascular event") # parent run
# # OR for the separate component run:
# target_outcomes <- c("Cardiovascular death",
#                      "Nonfatal myocardial infarction",
#                      "Nonfatal stroke",
#                      "General cardiac disorders",
#                      "Bleeding", "Thrombosis")
# selected <- subset(processed, var1 == target_product &
#                    var2 %in% target_outcomes)
# print(selected)

# Reproducibility: record R and package versions with your publication.
print(sessionInfo())
