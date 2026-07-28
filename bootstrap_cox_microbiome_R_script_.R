
# ============================================================
# Internal Cox bootstrap - HNSCC bacterial signature
# DSS, OS y DFS
# ============================================================

# Required packages
# install.packages(c("haven", "readxl", "survival", "dplyr", "openxlsx"))

library(haven)
library(readxl)
library(survival)
library(dplyr)
library(openxlsx)

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

# Option A: SPSS .sav file
# Change the filename if your file has a different name
file_sav <- "Variables microbioma Definitivos.sav"

# Option B: Excel .xlsx file
file_xlsx <- "Variables_en_binario_Microbioma.xlsx"

if (file.exists(file_sav)) {
  data <- read_sav(file_sav)
  cat("File loaded:", file_sav, "\n")
} else if (file.exists(file_xlsx)) {
  data <- read_excel(file_xlsx)
  cat("File loaded:", file_xlsx, "\n")
} else {
  stop("Neither the .sav nor the .xlsx file was found in the current working directory. Check getwd() and list.files().")
}

cat("\nAvailable variables:\n")
print(names(data))

# ------------------------------------------------------------
# 2. Convert required variables to numeric
# ------------------------------------------------------------

data <- data %>%
  mutate(
    Months = as.numeric(Months),
    DiseaseSpecificDeath = as.numeric(DiseaseSpecificDeath),
    OverallDeath = as.numeric(OverallDeath),
    Recurrence = as.numeric(Recurrence),
    MonthsToRecurrence = as.numeric(MonthsToRecurrence),
    BacterialSignature = as.numeric(BacterialSignature),
    TStage = as.numeric(TStage),
    NStage = as.numeric(NStage),
    TumorLocation = as.numeric(TumorLocation),
    Age = as.numeric(Age)
  )

# ------------------------------------------------------------
# 3. Define Cox models
# ------------------------------------------------------------
# DSS: time = Months; event = DiseaseSpecificDeath
# OS: time = Months; event = OverallDeath
# DFS: time = MonthsToRecurrence; event = Recurrence
#
# Bacterial signature = BacterialSignature

models <- list(
  DSS = list(
    time = "Months",
    event = "DiseaseSpecificDeath",
    covariates = c("TStage", "NStage", "TumorLocation", "BacterialSignature", "Age"),
    manuscript_hr = 6.215
  ),
  OS = list(
    time = "Months",
    event = "OverallDeath",
    covariates = c("TStage", "NStage", "TumorLocation", "BacterialSignature", "Age"),
    manuscript_hr = 5.417
  ),
  DFS = list(
    time = "MonthsToRecurrence",
    event = "Recurrence",
    covariates = c("TumorLocation", "BacterialSignature", "Age"),
    manuscript_hr = 6.842
  )
)

# ------------------------------------------------------------
# 4. Function to calculate the C-index in a specified dataset
# ------------------------------------------------------------
# reverse = TRUE because, in a Cox model, a higher linear predictor indicates higher risk
# and therefore a shorter time to the event.

calculate_cindex <- function(model, data_eval, time, event) {
  lp <- predict(model, newdata = data_eval, type = "lp")
  out <- concordance(
    as.formula(paste0("Surv(", time, ", ", event, ") ~ lp")),
    data = data_eval,
    reverse = TRUE
  )
  return(as.numeric(out$concordance))
}

# ------------------------------------------------------------
# 5. Main bootstrap function
# ------------------------------------------------------------

bootstrap_cox <- function(data, time, event, covariates,
                          endpoint, manuscript_hr,
                          n_boot = 1000, seed = 123) {

  set.seed(seed)

  vars <- c(time, event, covariates)
  model_data <- data[, vars] %>%
    na.omit()

  cox_formula <- as.formula(
    paste0(
      "Surv(", time, ", ", event, ") ~ ",
      paste(covariates, collapse = " + ")
    )
  )

  # Original model
  original_model <- coxph(cox_formula, data = model_data, x = TRUE)

  original_hr <- exp(coef(original_model)["BacterialSignature"])
  apparent_cindex <- summary(original_model)$concordance[1]

  n <- nrow(model_data)
  n_events <- sum(model_data[[event]] == 1, na.rm = TRUE)

  results <- data.frame(
    Endpoint = character(),
    Iteration = integer(),
    HR_BacterialSignature = numeric(),
    C_index_bootstrap = numeric(),
    C_index_original_eval = numeric(),
    Optimism = numeric(),
    Convergence = character()
  )

  for (i in 1:n_boot) {

    # Resampling with replacement
    indices <- sample(1:n, size = n, replace = TRUE)
    data_boot <- model_data[indices, ]

    fit_result <- tryCatch(
      {
        bootstrap_model <- coxph(cox_formula, data = data_boot, x = TRUE)

        hr_boot <- exp(coef(bootstrap_model)["BacterialSignature"])

        c_boot <- summary(bootstrap_model)$concordance[1]
        c_orig_eval <- calculate_cindex(bootstrap_model, model_data, time, event)

        data.frame(
          Endpoint = endpoint,
          Iteration = i,
          HR_BacterialSignature = hr_boot,
          C_index_bootstrap = c_boot,
          C_index_original_eval = c_orig_eval,
          Optimism = c_boot - c_orig_eval,
          Convergence = "OK"
        )
      },
      error = function(e) {
        data.frame(
          Endpoint = endpoint,
          Iteration = i,
          HR_BacterialSignature = NA,
          C_index_bootstrap = NA,
          C_index_original_eval = NA,
          Optimism = NA,
          Convergence = "Failed"
        )
      }
    )

    results <- rbind(results, fit_result)
  }

  # Retain valid iterations for the summary
  valid_results <- results %>%
    filter(!is.na(HR_BacterialSignature),
           is.finite(HR_BacterialSignature),
           !is.na(C_index_bootstrap),
           !is.na(C_index_original_eval))

  mean_optimism <- mean(valid_results$Optimism, na.rm = TRUE)
  corrected_cindex <- apparent_cindex - mean_optimism

  summary_table <- data.frame(
    Endpoint = endpoint,
    Time_variable = time,
    Event_variable = event,
    N_patients = n,
    N_events = n_events,
    Manuscript_HR = manuscript_hr,
    Recalculated_HR = as.numeric(original_hr),
    Bootstrap_HR_median = median(valid_results$HR_BacterialSignature, na.rm = TRUE),
    Bootstrap_HR_mean = mean(valid_results$HR_BacterialSignature, na.rm = TRUE),
    Bootstrap_95CI_lower = quantile(valid_results$HR_BacterialSignature, 0.025, na.rm = TRUE),
    Bootstrap_95CI_upper = quantile(valid_results$HR_BacterialSignature, 0.975, na.rm = TRUE),
    Apparent_C_index = as.numeric(apparent_cindex),
    Optimism_medio = mean_optimism,
    Optimism_corrected_C_index = corrected_cindex,
    Valid_iterations = nrow(valid_results),
    Failed_or_infinite_iterations = n_boot - nrow(valid_results)
  )

  return(list(
    original_model = original_model,
    summary_table = summary_table,
    iterations = results
  ))
}

# ------------------------------------------------------------
# 6. Run the bootstrap for DSS, OS, and DFS
# ------------------------------------------------------------

boot_DSS <- bootstrap_cox(
  data = data,
  time = models$DSS$time,
  event = models$DSS$event,
  covariates = models$DSS$covariates,
  endpoint = "DSS",
  manuscript_hr = models$DSS$manuscript_hr,
  n_boot = 1000,
  seed = 123
)

boot_OS <- bootstrap_cox(
  data = data,
  time = models$OS$time,
  event = models$OS$event,
  covariates = models$OS$covariates,
  endpoint = "OS",
  manuscript_hr = models$OS$manuscript_hr,
  n_boot = 1000,
  seed = 124
)

boot_DFS <- bootstrap_cox(
  data = data,
  time = models$DFS$time,
  event = models$DFS$event,
  covariates = models$DFS$covariates,
  endpoint = "DFS",
  manuscript_hr = models$DFS$manuscript_hr,
  n_boot = 1000,
  seed = 125
)

# ------------------------------------------------------------
# 7. Combine results
# ------------------------------------------------------------

final_summary <- bind_rows(
  boot_DSS$summary_table,
  boot_OS$summary_table,
  boot_DFS$summary_table
)

final_iterations <- bind_rows(
  boot_DSS$iterations,
  boot_OS$iterations,
  boot_DFS$iterations
)

print(final_summary)

# ------------------------------------------------------------
# 8. Export the final Excel workbook
# ------------------------------------------------------------

wb <- createWorkbook()

addWorksheet(wb, "Read first")
writeData(wb, "Read first", data.frame(
  Explanation = c(
    "Internal bootstrap of multivariable Cox models.",
    "A total of 1000 bootstrap iterations are performed for each endpoint.",
    "In each iteration, the cohort is resampled with replacement.",
    "The same Cox model is refitted and the HR for BacterialSignature is extracted.",
    "The C-index is analogous to the AUC/ROC, but adapted for survival analysis.",
    "The optimism-corrected C-index estimates the expected performance in new data.",
    "DSS: Surv(Months, DiseaseSpecificDeath).",
    "OS: Surv(Months, OverallDeath).",
    "DFS: Surv(MonthsToRecurrence, Recurrence)."
  )
))

addWorksheet(wb, "Summary")
writeData(wb, "Summary", final_summary)

addWorksheet(wb, "DSS")
writeData(wb, "DSS", boot_DSS$iterations)

addWorksheet(wb, "OS")
writeData(wb, "OS", boot_OS$iterations)

addWorksheet(wb, "DFS")
writeData(wb, "DFS", boot_DFS$iterations)

addWorksheet(wb, "Models_used")
models_used <- data.frame(
  Endpoint = c("DSS", "OS", "DFS"),
  Model = c(
    "coxph(Surv(Months, DiseaseSpecificDeath) ~ TStage + NStage + TumorLocation + BacterialSignature + Age)",
    "coxph(Surv(Months, OverallDeath) ~ TStage + NStage + TumorLocation + BacterialSignature + Age)",
    "coxph(Surv(MonthsToRecurrence, Recurrence) ~ TumorLocation + BacterialSignature + Age)"
  )
)
writeData(wb, "Models_used", models_used)

# Basic formatting
for (s in names(wb)) {
  freezePane(wb, s, firstRow = TRUE)
  setColWidths(wb, s, cols = 1:30, widths = "auto")
}

saveWorkbook(
  wb,
  file = "Bootstrap_Cox_1000_iterations_DSS_OS_DFS_R.xlsx",
  overwrite = TRUE
)

cat("\nExcel file generated: Bootstrap_Cox_1000_iterations_DSS_OS_DFS_R.xlsx\n")

