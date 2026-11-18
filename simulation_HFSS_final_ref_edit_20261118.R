#rm(list=ls())
#gc()
#setwd("/mnt/storage_fast4/Salt_England/IMPACTncd_Engl")

#Updated 20261104 (hallfast): START install into a writable user library. In Docker, R runs as your own user, so the
#system library (/usr/local/lib/R/site-library) is read-only and a plain install_github() fails ("not writable").
#The user library (R_LIBS_USER, created later by global.R anyway) is created here first; if it is not writable, the
#session's temporary folder is used instead. hallfast is then loaded so the forked simulation workers inherit it.
if (!requireNamespace("hallfast", quietly = TRUE)) {
  hallfast_lib <- Sys.getenv("R_LIBS_USER")
  if (nzchar(hallfast_lib)) dir.create(hallfast_lib, recursive = TRUE, showWarnings = FALSE)
  if (!nzchar(hallfast_lib) || file.access(hallfast_lib, mode = 2) != 0) {
    hallfast_lib <- file.path(tempdir(), "R_libs")
    dir.create(hallfast_lib, recursive = TRUE, showWarnings = FALSE)
  }
  .libPaths(c(hallfast_lib, .libPaths()))
  remotes::install_github("ChristK/hallfast", lib = hallfast_lib, upgrade = "never", build_vignettes = FALSE)
  rm(hallfast_lib)
}
invisible(loadNamespace("hallfast"))
#Updated 20261104 (hallfast): END install into a writable user library

source("./global.R") #Updated_by_Edi

#remotes::install_version("data.table", version = "1.16.4", repos = "http://cran.r-project.org")

source("./additional_costs.R")
source("./additional_prev_change.R")

IMPACTncd <- Simulation$new("./sim_design_Edi_final.yaml")

#Updated 20261104: START timing of each simulation step
### time_step() runs one step (a scenario run or an export), prints how long it took, and saves all
### timings so far to run_times.csv in the output folder (written after every step, so the file is
### up to date even if the job stops part-way).
run_times <- data.table(step = character(), started = character(), finished = character(), minutes = numeric())
time_step <- function(step, expr) {
  t0 <- Sys.time()
  message(sprintf("[%s] START %s", format(t0, "%Y-%m-%d %H:%M:%S"), step))
  force(expr)                                      # the step itself runs here
  t1 <- Sys.time()
  mins <- as.numeric(difftime(t1, t0, units = "mins"))
  run_times <<- rbind(run_times, data.table(step = step,
                                            started = format(t0, "%Y-%m-%d %H:%M:%S"),
                                            finished = format(t1, "%Y-%m-%d %H:%M:%S"),
                                            minutes = round(mins, 2)))
  message(sprintf("[%s] DONE  %s in %.1f min (%.2f h) | total so far %.2f h",
                  format(t1, "%Y-%m-%d %H:%M:%S"), step, mins, mins / 60, sum(run_times$minutes) / 60))
  out_dir <- IMPACTncd$design$sim_prm$output_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(run_times, file.path(out_dir, "run_times.csv"))
  invisible(NULL)
}
#Updated 20261104: END timing of each simulation step
#after running the line above, I got Warning message:
#file(s) not found in repo "ChristK/IMPACTncd_Engl"

start_runs <- 1
n_runs <- 100L
tax_rate <- 0.08 #change for sensitivity analysis here

Hall <- TRUE #if FALSE, will be using Christiansen & Garby energy model
Hall_PAL <- 1.6 #Updated 20261103 (hallfast): new parameter. physical activity level for the Hall model (hallfast default 1.6, SACN 2011; the old 'bw' package used 1.5)
Hall_pre_entry <- TRUE #Updated 20261104 (hallfast): new parameter. TRUE = people entering the simulated population after sc_year get the exposure they had before entry (adult years only); see hall_update_bmi()
Hall_within_year <- "mean" #Updated 20261104 (hallfast): new parameter. "mean" = annual BMI change is the mean of the start- and end-of-year changes; "end" = end-of-year change (20261103 behaviour)

#Updated 20261104 (hallfast): START revised helper function (pre-entry exposure; within-year averaging). Replaces the 20261103 version.
########## Hall model (hallfast) helper, used by all scenarios when Hall == TRUE
### Converts the policy-attributable change in energy intake (kcal/day, relative to baseline)
### into a change in BMI using the dynamic adult body weight model of Hall et al. (2011),
### via hallfast::hall_bw_long().
###
### - Every person is followed year by year from sc_year onwards: body composition (fat, lean,
###   glycogen, extracellular fluid, adaptive thermogenesis) carries over from one year to the
###   next, so the BMI effect accumulates over time, and the energy change may differ by year
###   (e.g. the compensation term, or waning).
### - Baseline body weight = BMI in the person's first modelled year x assumed height^2
###   (men 1.76 m, women 1.62 m; Health Survey for England 2021 means), as in the Christiansen
###   model.
### - Pre-entry exposure (pre_entry = TRUE, set by Hall_pre_entry). People who enter the
###   simulated population after sc_year (e.g. by ageing into the modelled age range) were
###   already exposed to the policy before they entered. For each of them,
###   n_pre = min(entry year - sc_year, floor(entry age - adult_age)) unobserved years are added
###   before their first row, with their entry-year weight and entry-year energy change, so they
###   enter with the effect accumulated since sc_year. Only adult years are added, because the
###   Hall et al. (2011) model is an adult model. Using the entry-year energy change is
###   conservative when the net reduction was larger in earlier years (growing grain/dairy
###   compensation term; waning scenarios). In-migrants were not exposed before arrival: if the
###   model includes them, exclude them from this step (or set Hall_pre_entry <- FALSE).
### - Within-year timing (within_year = "mean", set by Hall_within_year). The BMI change for
###   year y is the mean of the changes at the start and end of year y, a trapezoid
###   approximation of the average over the year (within ~1% of the exact person-time average
###   cumulatively over 10 years). "end" uses the change at the end of year y (the 20261103
###   behaviour), which front-loads the effect (about 1.7-fold in year 1, about 5% cumulatively
###   over 10 years).
### - The BMI change is added to that year's no-policy bmi_curr_xps. This additive
###   approximation ignores the small dependence of the response on later changes in no-policy
###   BMI (error within ~2% for drifts of -2 to +3 BMI units over 20 years).
### - No appetite feedback (appetite_k = 0), as the energy change is already a net change
###   (it includes compensation).
### - pre_entry = FALSE with within_year = "end" reproduces the 20261103 helper exactly.
###
### Arguments
###   pop         data.table with pid, year, age, sex ("men"/"women"), bmi_curr_xps and
###               energy_col; modified by reference.
###   energy_col  name of the column with the policy-attributable change in energy intake
###               (kcal/day, negative = reduction).
###   sc_year     first policy year (same short-year format as pop$year, e.g. 26L).
###   PAL         physical activity level passed to hallfast (default Hall_PAL).
###   pre_entry   add pre-entry exposure for people whose first modelled year is after sc_year
###               (default Hall_pre_entry).
###   adult_age   youngest age from which pre-entry exposure is added (default 18).
###   within_year "mean" or "end" (default Hall_within_year; see above).
###   solver      function with the interface of hallfast::hall_bw_long(); exposed only so the
###               helper can be tested with another solver.
### Returns pop invisibly (modified by reference).
hall_update_bmi <- function(pop, energy_col, sc_year, PAL = Hall_PAL,
                            pre_entry = Hall_pre_entry, adult_age = 18L,
                            within_year = Hall_within_year,
                            solver = hallfast::hall_bw_long) {
  within_year <- match.arg(within_year, c("mean", "end"))
  stopifnot(is.logical(pre_entry), length(pre_entry) == 1L, !is.na(pre_entry))
  idx <- pop[year >= sc_year, which = TRUE]
  if (length(idx) == 0L) return(invisible(pop))
  hd <- data.table(
    row  = idx,
    pid  = pop$pid[idx],
    year = pop$year[idx],
    age  = pop$age[idx],
    sex  = as.character(pop$sex[idx]),
    bmi  = pop$bmi_curr_xps[idx],
    dEI  = as.numeric(pop[[energy_col]][idx])
  )
  hd[is.na(dEI) | !is.finite(dEI), dEI := 0]
  setorder(hd, pid, year) # rows of each person contiguous and in chronological order
  hd[, ht := fifelse(sex == "men", 1.76, 1.62)]
  hd[, bw := bmi * ht^2]

  ## Pre-entry exposure: add unobserved adult years before the first row of late entrants
  if (pre_entry) {
    first <- hd[hd[, .I[1L], by = pid]$V1]          # first modelled row of each person
    first[, n_pre := pmax(0L, pmin(as.integer(year - sc_year),
                                   as.integer(floor(age - adult_age))))]
    first <- first[!is.na(n_pre) & n_pre > 0L]
    if (nrow(first) > 0L) {
      pre <- first[rep(seq_len(.N), n_pre)]           # one row per pre-entry year
      pre[, k := seq_len(.N), by = pid]               # 1 .. n_pre
      pre[, `:=`(row  = NA_integer_,                  # NA marks rows that are not in pop
                 year = year - n_pre + k - 1L,        # the years between sc_year and entry
                 age  = age  - n_pre + k - 1L)]       # correspondingly younger
      # bw, bmi, ht, sex and dEI are carried from the entry row (weight assumed stable
      # before entry; entry-year energy change applied in each pre-entry year)
      hd <- rbindlist(list(pre[, names(hd), with = FALSE], hd), use.names = TRUE)
      setorder(hd, pid, year)
    }
  }

  ## Baseline weight = first row of each person (an added pre-entry row for late entrants)
  hd[, bw0 := bw[match(pid, pid)]]
  ok <- hd[, !is.na(bw0) & is.finite(bw0) & bw0 > 0]
  hd[, bw_end := NA_real_]
  hd[ok, bw_end := solver(pid, bw, ht, age, sex, dEI, PAL = PAL)]

  ## Annual BMI change
  if (within_year == "mean") {
    # Weight at the start of each year = weight at the end of the previous year;
    # for each person's first row it is the baseline weight.
    hd[, bw_start := shift(bw_end)]
    hd[c(TRUE, pid[-1L] != pid[-.N]), bw_start := bw0]
    hd[, dbmi := ((bw_start + bw_end) / 2 - bw0) / ht^2]
  } else {
    hd[, dbmi := (bw_end - bw0) / ht^2]
  }
  hd[is.na(dbmi) | !is.finite(dbmi), dbmi := 0]

  hd <- hd[!is.na(row)]                               # drop the added pre-entry rows
  pop[hd$row, bmi_curr_xps := bmi_curr_xps + hd$dbmi]
  invisible(pop)
}
#Updated 20261104 (hallfast): END revised helper function

#Updated 20261103: START per-parameter seeds (common random numbers across scenarios)
### Each scenario defines draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }
### and every random parameter is drawn with its own fixed k, so for the same mc the same parameter
### gets the same random value in every scenario (regardless of which other parameters a scenario draws).
### Differences between scenarios then reflect the policies, not different parameter draws.
### k:  1 elasticity_HFSS and r_new (tax effect by income; same k so both scenarios use the same random number)   2 (not used)   3 sbp_effect (BMI->SBP)   4 tchol_effect (BMI->tchol)
###     5 r_new_sbp (salt->SBP)   6 elasticity_fruit   7 elasticity_veg
###     8 lbl_eff_energy     9 lbl_eff_salt    10 lbl_eff_energy1   11 lbl_eff_energy2   12 lbl_eff_salt1   13 lbl_eff_salt2
###    14 ref_eff_energy    15 ref_eff_salt    16 tax_eff_fruit (subsidy)   17 tax_eff_veg (subsidy)
### Keep k < 100 and use a new k for any new parameter.
#Updated 20261103: END per-parameter seeds

#Updated 20261103: START cap HFSS intake at total intake
### Set cap_hfss_to_total <- FALSE to use the HFSS values as they are (the previous behaviour).
cap_hfss_to_total <- TRUE
hfss_capped <- function(hfss, total) {
  if (!cap_hfss_to_total) return(hfss)
  fifelse(is.na(total), hfss, pmin(hfss, total)) # if total is missing, keep the HFSS value
}
#Updated 20261103: END cap HFSS intake at total intake

########## BASELINE SCENARIO
  #del_logs()$
  #del_outputs()$
time_step("sc0", IMPACTncd$
  run(mc = start_runs:n_runs, multicore = TRUE, scenario_nam = "sc0")) #Updated 20261104: timing
  #export_summaries(multicore = TRUE)





########## PATHWAYS THROUGH ENERGY AND SALT

########## Translating changes in energy kcal into body weight/BMI
########## Hall Method was used (Hall et al. 2011, https://doi.org/10.1016/S0140-6736(11)60812-X), solved with the hallfast package https://github.com/ChristK/hallfast
########## Dynamic year-by-year body weight trajectory per person from sc_year onwards (see hall_update_bmi() above) #Updated 20261103 (was: https://github.com/INSP-RH/bw lookup table) #Updated 20261104: with pre-entry exposure for late entrants and within-year averaging

########## Christiansen & Garby https://doi.org/10.1046/j.1365-2362.2002.01036.x
### men: ∆BWss = 17·7 ∆(EI/PAL)
### women: ∆BWss = 20·7 ∆(EI/PAL)
### EI= energy intake MJ/day --> 1 kcal = 4.2/1000 MJ
### PAL= physical activity level ---> kept constant at 1.5 https://openknowledge.fao.org/server/api/core/bitstreams/62ae7aeb-9536-4e43-b2d0-55120e662824/content


########## Translating changes in salt g/d into SBP
########## Formula Mozaffarian et al. https://www.nejm.org/doi/full/10.1056/NEJMoa1304127#APPNEJMoa1304127SUP (in the supp. text)
### "Y = α+ B1X1 + B2X2 + B3X3, where Y = the BP change in each trial in mm Hg,
### standardized to a sodium reduction of 100 mmol/d; B1 = the mean age of subjects in each trial in years, centered
### at age 50; B2 = whether the trial was conducted in normotensive or hypertensive subjects; B3 = whether the trial
### was conducted in non-Blacks or Blacks; and α = the effect of a 100 mmol/d sodium reduction on systolic BP
### when mean subject age=50, subjects were normotensive, and subjects were not Black. Each trial was weighted
### by the inverse variance of its BP change, standardized to a sodium reduction of 100 mmol/d.The final coefficients (SE) were:
### α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 ).
### In other words, for a white, normotensive population at age 50, each 100 mmol/d sodium reduction lowered
### systolic BP by 3.74 mm Hg (95%CI: 2.29, 5.18). For every year above or below age 50, there was 0.105 mm
### Hg (95%CI: 0.047, 0.164) larger or smaller BP reduction, respectively: e.g., for normotensives at age 25, a 100
### mmol/d reduction lowered BP by 1.12 mm Hg (3.74 – 25 yrs*0.105/yr); and
### at age 70, by 5.84 mm Hg (3.74 + 20 yrs*0.105/yr). Additionally, there was 2.49 mm Hg greater reduction among Blacks (95%CI: 0.13, 4.85),
### and 1.87 mm Hg greater reduction among hypertensives (95%CI: 0.12, 3.63)."
### ----> 100 mmol/d (2.3 g/d sodium, 5.85 g/d salt)





######## CONSUMER RESPONSE ONLY ##################################################################################################################################################
######## TAXATION ##################################################################################################################################################
### Based on a study in the Scottish context, an increase 10% price of HFSS products will reduce the consumption by 8% (an average of 6%+10%/2))
### https://doi.org/https://doi.org/10.1016/j.jclepro.2023.138564

################# Scenario 1: Taxation 8% (Consumer response only) [cons_sc1_tax] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- draw(1, runif(1, min = 0.6, max = 1.0)) #per 1% increase in price #based on discussion with Chris 20250805  #Updated 20261103: draw()
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~5.8 kcal for energy compensation in 2019 (10% increase in 10 years)
    sp$pop[year>= sc_year, energy_tax:= energy_diet*tax_eff_energy + 5.8 * (1 + 0.10)^((year - 19L) / 10L)]
    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_tax", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_tax") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_tax*(4.2/1000))/1.5),
                                   20.7 * ((energy_tax*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_tax", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
           ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)] #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~0.015 g for salt compensation in 2019 (18% reduction in 10 years)
    sp$pop[year>= sc_year, salt_tax:= salt_diet_g*tax_eff_salt + 0.015 * (1 - 0.18)^((year - 19L) / 10L)]
    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_tax >0, sbp_change_tax :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*0.60, sbp_change_tax :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_tax <0, sbp_change_tax := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_tax) & (sbp_change_tax >= 0),
                                   sbp_curr_xps - sbp_change_tax,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_tax", "salt_diet_g", "sbp_change_tax") := NULL]

 ## Reduction in fruits and vegetables due to taxation
    #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
    #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   137.4731 549.0456 #this is alpha, and beta
    elasticity_fruit <- draw(6, rbeta(1, 137.4731, 549.0456))*(-1)  #Updated 20261103: draw()

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)))] #tax_eff_fruit is negative #Updated 20261103


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   84.65086 84.65087 #this is alpha, and beta
    elasticity_veg <- draw(7, rbeta(1, 84.65086, 84.65087))*(-1)  #Updated 20261103: draw()

    tax_eff_veg <- tax_rate*elasticity_veg #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)))] #tax_eff_veg is negative #Updated 20261103

    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_tax", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


################# Scenario 2: Taxation 8% (Consumer response only, differences in the effects by INCOME) [cons_sc1_tax_effINC] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    #https://doi.org/10.1371/journal.pmed.1002057
    #Based on evidence from Mexico, low SES households purchased 10.2% less taxed foods than expected (p < 0.05),
    #whereas medium SES households purchased 5.8% less taxed foods (p < 0.05), and high-income households’ purchases did not change (2.3%, p> 0.05)
    # as used in this study https://doi.org/10.1186/s12889-025-21745-0

    #Based on the same study https://doi.org/10.1371/journal.pmed.1002057, the absolute reductions are as follows:
    #low SES: -44 (-72, -16), with baseline 434, --> equivalent to -10.2% (-16.6%, -3.7%)
    #medium SES: -28 (-46, -11), with baseline 491, --> equivalent to -5.8% (-9.4%, -2.2%)
    #high SES: -13 (-39, 12), with baseline 569, --> equivalent to -2.3% (-6.9%, 2.1%)

    #p <- c(0.025, 0.5, 0.975)
    #q1 <- c(0.037, 0.102, 0.166)
    #q2 <- c(0.022, 0.058, 0.094)
    #p0 <- c(0, 0.5, 0.975)
    #q3 <- c(0, 0.023, 0.069)
    #
    #params <- get.beta.par(q = q1, p = p)
    #$par
    #[1] 11.46264 98.38558
    #params <- get.beta.par(q = q2, p = p)
    #$par
    #[1]  12.46944 197.57645
    #params <- get.beta.par(q = q3, p = p0)
    #$par
    #[1]  2.307398 84.614135
    r_new <- draw(1, runif(1, min = 0, max = 1)) # to make it consistent  #Updated 20261103: draw(); k = 1 (same random number as elasticity_HFSS), so the income-specific tax effects sit at the same point of their uncertainty range as the main tax scenario

    tax_eff_energy1 <- (qbeta(r_new, 11.46264, 98.38558))*(-1)
    if (tax_eff_energy1 > 0) tax_eff_energy1 <- 0 #negative

    tax_eff_energy3 <- (qbeta(r_new, 12.46944, 197.57645))*(-1)
    if (tax_eff_energy3 > 0) tax_eff_energy3 <- 0 #negative

    tax_eff_energy5 <- (qbeta(r_new, 2.307398, 84.614135))*(-1)
    if (tax_eff_energy5 > 0) tax_eff_energy5 <- 0 #negative

    tax_eff_energy2 <- (tax_eff_energy1 + tax_eff_energy3)/2
    tax_eff_energy4 <- (tax_eff_energy3 + tax_eff_energy5)/2

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    # income levels: "1 Lowest", "2", "3", "4", "5 Highest" #Updated 20261103
    # 379042    349334    362763    369401    352700
    # sp$pop[year>= sc_year, energy_tax:= energy_diet*tax_eff_energy]
    #Updated 20261104: stop if the income labels differ from those used below; a mismatch would silently give no effect to the unmatched quintiles (and swap others)
    stopifnot("income labels must be '1 Lowest', '2', '3', '4', '5 Highest'" =
                setequal(unique(as.character(na.omit(sp$pop$income))), c("1 Lowest", "2", "3", "4", "5 Highest")))
    if (sp$pop[year >= sc_year, anyNA(income)]) warning("Some person-years have missing income; they get no tax effect in cons_sc1_tax_effINC.") #Updated 20261104
    sp$pop[year>= sc_year & income == "1 Lowest", energy_tax:= energy_diet*tax_eff_energy1] #income-specific effect (Mexico evidence) #Updated 20261103
    sp$pop[year>= sc_year & income == "2", energy_tax:= energy_diet*tax_eff_energy2]
    sp$pop[year>= sc_year & income == "3", energy_tax:= energy_diet*tax_eff_energy3]
    sp$pop[year>= sc_year & income == "4", energy_tax:= energy_diet*tax_eff_energy4]
    sp$pop[year>= sc_year & income == "5 Highest", energy_tax:= energy_diet*tax_eff_energy5]

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~5.8 kcal for energy compensation in 2019 (10% increase in 10 years)
    sp$pop[year >= sc_year,
           energy_tax := energy_tax + 5.8 * (1 + 0.10)^((year - 19L) / 10L)]

    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_tax", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_tax") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_tax*(4.2/1000))/1.5),
                                   20.7 * ((energy_tax*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_tax", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
    )]

    tax_eff_salt1 <- tax_eff_energy1 # the effect was assumed the same (the proportional change)
    tax_eff_salt2 <- tax_eff_energy2
    tax_eff_salt3 <- tax_eff_energy3
    tax_eff_salt4 <- tax_eff_energy4
    tax_eff_salt5 <- tax_eff_energy5

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)] #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #sp$pop[year>= sc_year, salt_tax:= salt_diet_g*tax_eff_salt]
    sp$pop[year>= sc_year & income == "1 Lowest", salt_tax:= salt_diet_g*tax_eff_salt1] #income-specific effect (Mexico evidence) #Updated 20261103
    sp$pop[year>= sc_year & income == "2", salt_tax:= salt_diet_g*tax_eff_salt2]
    sp$pop[year>= sc_year & income == "3", salt_tax:= salt_diet_g*tax_eff_salt3]
    sp$pop[year>= sc_year & income == "4", salt_tax:= salt_diet_g*tax_eff_salt4]
    sp$pop[year>= sc_year & income == "5 Highest", salt_tax:= salt_diet_g*tax_eff_salt5]

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~0.015 g for salt compensation in 2019 (18% reduction in 10 years)
    sp$pop[year >= sc_year,
           salt_tax := salt_tax + 0.015 * (1 - 0.18)^((year - 19L) / 10L)]

    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_tax >0, sbp_change_tax :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*0.60, sbp_change_tax :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_tax <0, sbp_change_tax := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_tax) & (sbp_change_tax >= 0),
                                   sbp_curr_xps - sbp_change_tax,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_tax", "salt_diet_g", "sbp_change_tax") := NULL]

 ## Reduction in fruits and vegetables due to taxation
    #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
    #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   137.4731 549.0456 #this is alpha, and beta
    elasticity_fruit <- draw(6, rbeta(1, 137.4731, 549.0456))*(-1)  #Updated 20261103: draw()

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)))] #tax_eff_fruit is negative #Updated 20261103


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   84.65086 84.65087 #this is alpha, and beta
    elasticity_veg <- draw(7, rbeta(1, 84.65086, 84.65087))*(-1)  #Updated 20261103: draw()

    tax_eff_veg <- tax_rate*elasticity_veg #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)))] #tax_eff_veg is negative #Updated 20261103

    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_tax_effINC", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_effINC")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


################# Scenario 3: Taxation 8% (Consumer response only, without reduction in fruits and veggies) [cons_sc1_tax_nfvg] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- draw(1, runif(1, min = 0.6, max = 1.0)) #per 1% increase in price #based on discussion with Chris 20250805  #Updated 20261103: draw()
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~5.8 kcal for energy compensation in 2019 (10% increase in 10 years)
    sp$pop[year>= sc_year, energy_tax:= energy_diet*tax_eff_energy + 5.8 * (1 + 0.10)^((year - 19L) / 10L)]
    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_tax", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_tax") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_tax*(4.2/1000))/1.5),
                                   20.7 * ((energy_tax*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_tax", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)] #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~0.015 g for salt compensation in 2019 (18% reduction in 10 years)
    sp$pop[year>= sc_year, salt_tax:= salt_diet_g*tax_eff_salt + 0.015 * (1 - 0.18)^((year - 19L) / 10L)]
    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_tax >0, sbp_change_tax :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*0.60, sbp_change_tax :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_tax <0, sbp_change_tax := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_tax) & (sbp_change_tax >= 0),
                                   sbp_curr_xps - sbp_change_tax,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_tax", "salt_diet_g", "sbp_change_tax") := NULL]

 # ## Reduction in fruits and vegetables due to taxation
 #    #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
 #    #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
 #    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
 #    # install.packages("rriskDistributions")
 #    # library(rriskDistributions)
 #    # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
 #    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
 #    #
 #    # params <- get.beta.par(q = q, p = p)
 #    # params
 #    #
 #    # $par
 #    # [1]   137.4731 549.0456 #this is alpha, and beta
 #    elasticity_fruit <- rbeta(1, 137.4731, 549.0456)*(-1)
 #
 #    tax_eff_fruit <- tax_rate*elasticity_fruit #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
 #    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative
 #
 #    sp$pop[year >= sc_year,
 #           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)))] #tax_eff_fruit is negative #Updated 20261103
 #
 #
 #    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
 #    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
 #    # install.packages("rriskDistributions")
 #    # library(rriskDistributions)
 #    # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
 #    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
 #    #
 #    # params <- get.beta.par(q = q, p = p)
 #    # params
 #    #
 #    # $par
 #    # [1]   84.65086 84.65087 #this is alpha, and beta
 #    elasticity_veg <- rbeta(1, 84.65086, 84.65087)*(-1)
 #
 #    tax_eff_veg <- tax_rate*elasticity_veg #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
 #    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative
 #
 #    sp$pop[year >= sc_year,
 #           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)))] #tax_eff_veg is negative #Updated 20261103

    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_tax_nfvg", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_nfvg")) #Updated 20261104: timing
#
# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)





######## NUTRITION WARNING LABELLING ##################################################################################################################################################
### MAIN
### Song et al. NW labelling reduced energy purchased by 3.8% (95% CI: [0.9%, 6.8%]) per 100 grams of foods,
### but did not significantly differ when compared to MLT labelling in reducing total energy intake per 100 grams: 0.8% (95% CI: [-2.0%, 3.7%]).
### No significant effects were reported for NW labelling in reducing salt per 100 grams of food when compared to either control (2.1%; 95% CI: [-4.0%; 8.2%])
### or MLT labelling (0.8%; 95% CI: [-4.9%; 6.5%]). Therefore, NW labelling may not outperform MLT labelling in reducing energy intake per 100-gram unit.

### SENSITIVITY
### Song et al. NW labelling reduced energy intake by 12.9% (95% CI: [8%; 17.9%]) and outperformed MLT labelling by 6.4% (95% CI: [0.4%; 12.5%]) in reducing energy intake
### Non-significant effects for i) NW labelling (vs. control) (7.8%; 95% CI: [-1.2%; 16.8%]) and ii) NW (vs. MLT labelling) (1.3%; 95% CI: [-7.7%; 10.3%]) in reducing total salt (grams).
### https://doi.org/10.1371/journal.pmed.1003765

### Croker et al. NW labelling reduced energy intake by 4.43 kcal (95% CI: [0.12, 8.74]) per 100 gram of food
### and sodium by 33.78 mg (95% CI: [8.16, 59.40]) per 100 gram of food (will be divided by 393.4 to translate into salt (grams))
### https://doi.org/10.1111/jhn.12758

################# Scenario 4: Song et al. (Consumer response only) [cons_sc1_lbl_sn] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.009, 0.038, 0.068)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   8.893301 217.128711 #this is alpha, and beta
    lbl_eff_energy <- draw(8, rbeta(1, 8.893301, 217.128711))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
    #prop_packaged <- 0.80
    prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels

    #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the MLT become mandatory
    #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
    current_coverage_MLT <- 0.75
    extra_coverage_MLT <- 1- current_coverage_MLT

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)*prop_packaged] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS energy explicitly as food weight (g) x energy density (kcal per gram of HFSS food).
    #Song et al. effect is a % change in energy per 100 g (i.e. in energy density), so: change = weight x density x % effect (= energy_diet x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not kcal.
    sp$pop[, energy_density := fifelse(hfss_food_weight > 0, energy_diet / hfss_food_weight, 0)] #kcal per gram of HFSS food
    sp$pop[year>= sc_year, energy_lbl:= (hfss_food_weight*energy_density)*(extra_coverage_MLT*lbl_eff_energy)]
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative
    sp$pop[, energy_density := NULL] #Updated 20261103

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_lbl", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_lbl") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_lbl*(4.2/1000))/1.5),
                                   20.7 * ((energy_lbl*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_lbl", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.021, 0.082)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   1.454906 53.419114 #this is alpha, and beta
    lbl_eff_salt <- draw(9, rbeta(1, 1.454906, 53.419114))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)*prop_packaged] #see above #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS salt explicitly as food weight (g) x salt density (g salt per gram of HFSS food).
    #Song et al. effect is a % change in salt per 100 g (i.e. in salt density), so: change = weight x density x % effect (= salt_diet_g x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not grams of salt.
    sp$pop[, salt_density := fifelse(hfss_food_weight > 0, salt_diet_g / hfss_food_weight, 0)] #g salt per gram of HFSS food
    sp$pop[year>= sc_year, salt_lbl:= (hfss_food_weight*salt_density)*(extra_coverage_MLT*lbl_eff_salt)]
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative
    sp$pop[, salt_density := NULL] #Updated 20261103

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_lbl >0, sbp_change_lbl :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*prop_packaged*0.60, sbp_change_lbl :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_lbl <0, sbp_change_lbl := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_lbl) & (sbp_change_lbl >= 0),
                                   sbp_curr_xps - sbp_change_lbl,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_lbl", "salt_diet_g", "sbp_change_lbl") := NULL]
    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_lbl_sn", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_lbl_sn")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


################# Scenario 5: Song et al. (Consumer response only, sensitivity reduction in total intakes) [cons_sc1_lbl_sn_sensi] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    #lbl_eff_energy1 <- rnorm(1, mean = -0.129, sd = (-0.08 - (-0.179))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.08, 0.129, 0.179)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   24.06169 160.75193 #this is alpha, and beta
    lbl_eff_energy1 <- draw(10, rbeta(1, 24.06169, 160.75193))*(-1)  #Updated 20261103: draw()
    if (lbl_eff_energy1 > 0) lbl_eff_energy1 <- 0 #negative

    #lbl_eff_energy2 <- rnorm(1, mean = -0.064, sd = (-0.004 - (-0.125))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.004, 0.064, 0.125)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   6.360689 88.532728 #this is alpha, and beta #Updated 20261103 (was copied from lbl_eff_energy1)
    lbl_eff_energy2 <- draw(11, rbeta(1, 6.360689, 88.532728))*(-1)  #Updated 20261103: draw()
    if (lbl_eff_energy2 > 0) lbl_eff_energy2 <- 0 #negative

    #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
    #prop_packaged <- 0.80
    prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels

    #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the MLT become mandatory
    #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
    current_coverage_MLT <- 0.75
    extra_coverage_MLT <- 1- current_coverage_MLT

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)*prop_packaged] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    #Based on different effect sizes 12.9% vs 6.4%, we will apply 6.4% to 75% of the products currently featuring MLT as the effect of NW outperform NL
    sp$pop[year>= sc_year, energy_lbl:= energy_diet*((current_coverage_MLT*lbl_eff_energy2) + (extra_coverage_MLT*lbl_eff_energy1))]
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_lbl", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_lbl") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_lbl*(4.2/1000))/1.5),
                                   20.7 * ((energy_lbl*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_lbl", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    #lbl_eff_salt1 <- rnorm(1, mean = -0.078, sd = (0.168 - (-0.012))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.078, 0.168)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   4.584322 50.638770 #this is alpha, and beta
    lbl_eff_salt1 <- draw(12, rbeta(1, 4.584322, 50.638770))*(-1)  #Updated 20261103: draw()
    if (lbl_eff_salt1 > 0) lbl_eff_salt1 <- 0 #negative

    #lbl_eff_salt2 <- rnorm(1, mean = -0.013, sd = (0.103 - (-0.077))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.013, 0.103)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   0.6301193 26.3996677  #this is alpha, and beta
    lbl_eff_salt2 <- draw(13, rbeta(1, 0.6301193, 26.3996677))*(-1)  #Updated 20261103: draw()
    if (lbl_eff_salt2 > 0) lbl_eff_salt2 <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)*prop_packaged] #see above #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_lbl:= salt_diet_g*((current_coverage_MLT*lbl_eff_salt2) + (extra_coverage_MLT*lbl_eff_salt1))]
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_lbl >0, sbp_change_lbl :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*prop_packaged*0.60, sbp_change_lbl :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_lbl <0, sbp_change_lbl := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_lbl) & (sbp_change_lbl >= 0),
                                   sbp_curr_xps - sbp_change_lbl,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_lbl", "salt_diet_g", "sbp_change_lbl") := NULL]
    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_lbl_sn_sensi", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_lbl_sn_sensi")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


# ################# Scenario 2: Croker et al. (Consumer response only)
# IMPACTncd$update_primary_prevention_scn(
#   function(sp) {
#     #setDT(sp$pop)
#     # Ensure that .Random.seed exists (it is created when a random number is generated)
#     if (!exists(".Random.seed", envir = .GlobalEnv)) {
#       runif(1)
#     }
#
#     # Save the current random seed state
#     saved_seed <- .Random.seed
#
#     # Change the seed (this will change the state of the RNG)
#     set.seed(93684 + sp$mc_aggr)
#
#  ## Reduction in energy
#     lbl_eff_energy <- rnorm(1, mean = -4.43, sd = (8.74 - 0.12)/3.92) #abs reduction per 100 gram of food
#     if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative
#
#     #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
#     #prop_packaged <- 0.80
#     prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels
#
#     #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
#     #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
#     #so, 25% of the products will be labeled if the MLT become mandatory
#     #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
#     current_coverage_MLT <- 0.75
#     extra_coverage_MLT <- 1- current_coverage_MLT
#
#     sp$pop[, energy_diet:= hfss_energy*prop_packaged]
#     sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive
#
#     #We will estimate the impact from 2026
#     sc_year <- 26L
#     #Based on the effect per 100 grams of foods, we calculated the reduction based on consumption of HFSS
#     sp$pop[year>= sc_year, energy_lbl:= (hfss_food_weight/100*prop_packaged*extra_coverage_MLT)*lbl_eff_energy] #multiply with total grams consumed from HFSS
#     #test <- sp$pop[energy_lbl*(-1) > energy_diet, c("energy_lbl", "energy_diet")] # due to energy_lbl was calculated based on food_grams
#     sp$pop[energy_lbl*(-1) > energy_diet*extra_coverage_MLT, energy_lbl := energy_diet*extra_coverage_MLT*(-1)] #make sure the reduction is not > energy from HFSS
#     sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative
#
#     sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation
#
#     if (Hall) {
#     #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
#     #Translating change in energy intake (kcal/day) into change in body weight using the
#     #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
#     #Each person's body weight trajectory is simulated year by year from sc_year onwards
#     #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
#     hall_update_bmi(sp$pop, energy_col = "energy_lbl", sc_year = sc_year)
#     sp$pop[, c("energy_diet", "energy_lbl") := NULL]
#     #Updated 20261103 (hallfast): END
#
#     } else{
#       #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
#       sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_lbl*(4.2/1000))/1.5),
#                                    20.7 * ((energy_lbl*(4.2/1000))/1.5))]
#
#       ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
#       sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
#       sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4
#
#       sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
#       sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)
#
#       sp$pop[year >= sc_year,
#              bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
#                                      bmi_curr_xps* relative_change_bw,
#                                      bmi_curr_xps)]
#       #sp$pop[, summary(bmi_curr_xps)]
#       sp$pop[, c("energy_diet", "energy_lbl", "change_bw", "weight_assum", "relative_change_bw") := NULL]
#       # Restore the original seed state
#     }
#
#     sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
#     sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]
#
#     #adding mediated effect through SBP and Tchol #20250901
#     #https://doi.org/10.1111/jch.14661
#     #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
#     #Only applying effect to those with BMI >= 25
#
#     #https://doi.org/10.1210/clinem/dgaa673
#     #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
#     #which is -0.229mmol/L (-0.305, -0.154)
#
#     sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
#     sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
#     sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed
#
#     tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
#     tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
#     tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed
#
#     sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
#            `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
#                  tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
#     ]
#
#     sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]
#
#  ## Reduction in salt
#     sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
#                   black = as.integer(ethnicity=="black african")
#                   )]
#
#     lbl_eff_salt <- (rnorm(1, mean = -33.78, sd = (59.40 - 8.16)/3.92))/393.4 # this is sodium in mg and translate into salt (grams) #abs reduction per 100 gram of food
#     if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative
#
#     sp$pop[, salt_diet_g:= hfss_salt*prop_packaged]
#     sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive
#
#     sp$pop[year>= sc_year, salt_lbl:= (hfss_food_weight/100*prop_packaged*extra_coverage_MLT)*lbl_eff_salt] #multiply with total grams consumed from HFSS
#     sp$pop[salt_lbl*(-1) > salt_diet_g*extra_coverage_MLT, salt_lbl := salt_diet_g*extra_coverage_MLT*(-1)] #make sure the reduction is not > salt from HFSS
#     sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative
#
#     #Translating reduction in salt (gram) into reduction in SBP (mmHG)
#     #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
#     r_new_sbp <- runif(1, min = 0, max = 1)
#     a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
#     b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
#     b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
#     b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)
#
#     sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
#     sp$pop[salt_lbl >0, sbp_change_lbl :=  0] #no change in SBP if no reduction in salt
#     #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
#     sp$pop[salt_diet_g <=3.75*prop_packaged*0.60, sbp_change_lbl :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
#     sp$pop[sbp_change_lbl <0, sbp_change_lbl := 0]
#
#     sp$pop[year >= sc_year,
#            sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_lbl) & (sbp_change_lbl >= 0),
#                                    sbp_curr_xps - sbp_change_lbl,
#                                    sbp_curr_xps)]
#     #sp$pop[, summary(sbp_curr_xps)]
#     sp$pop[, c("hypertensive", "black", "salt_lbl", "salt_diet_g", "sbp_change_lbl") := NULL]
#     # Restore the original seed state
#     .Random.seed <- saved_seed
#   }
# )
#
# IMPACTncd$
#   run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc2_lbl_cr")
# #
# # print("SUMMARY EXPORTS!!!")
# #
# # IMPACTncd$export_summaries(multicore = TRUE)


######## COMBINED EFFECTS OF TAXATION AND NUTRITION WARNING LABELLING ##################################################################################################################################################
################# Scenario 6: Taxation 8% and Song et al. (Consumer response only) [cons_sc1_tax_lbl_sn] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- draw(1, runif(1, min = 0.6, max = 1.0)) #per 1% increase in price #based on discussion with Chris 20250805  #Updated 20261103: draw()
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    # label Song et al.
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.009, 0.038, 0.068)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   8.893301 217.128711 #this is alpha, and beta
    lbl_eff_energy <- draw(8, rbeta(1, 8.893301, 217.128711))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
    #prop_packaged <- 0.80
    prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels

    #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the MLT become mandatory
    #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
    current_coverage_MLT <- 0.75
    extra_coverage_MLT <- 1- current_coverage_MLT

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~5.8 kcal for energy compensation in 2019 (10% increase in 10 years)
    # The combined effect can be calculated as
    # 1 – [(1+PRi)*(1+PRj)] PR= percent reduction < 0
    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS energy explicitly as food weight (g) x energy density (kcal per gram of HFSS food).
    #Song et al. effect is a % change in energy per 100 g (i.e. in energy density), so: change = weight x density x % effect (= energy_diet x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not kcal.
    sp$pop[, energy_density := fifelse(hfss_food_weight > 0, energy_diet / hfss_food_weight, 0)] #kcal per gram of HFSS food
    sp$pop[year>= sc_year, energy_combined:= ((hfss_food_weight*energy_density)*(1+tax_eff_energy)*(1+(extra_coverage_MLT*lbl_eff_energy)) - (hfss_food_weight*energy_density) + 5.8 * (1 + 0.10)^((year - 19L) / 10L))] #Updated 20261104: labelling factor no longer applied to the grain/dairy substitution term (substitution is due to the tax only; grains/dairy are not HFSS and not labelled)
    #sp$pop[year>= sc_year, energy_combined:= energy_diet*(((1+tax_eff_energy)*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy)))-1)] #negative value = reduction # the same as above #Updated 20261103
    sp$pop[is.na(energy_combined) | is.infinite(energy_combined) | energy_combined > 0, energy_combined := 0] #negative
    sp$pop[, energy_density := NULL] #Updated 20261103

    # sp$pop[year>= sc_year, energy_tax:= energy_diet*tax_eff_energy]
    # sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative
    #
    #
    # sp$pop[year >= sc_year,
    #        energy_tax := energy_tax + 5.8 * (1 + 0.10)^((year - 19L) / 10L)]
    #
    # sp$pop[year>= sc_year, energy_lbl:= (energy_diet+energy_tax)*(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy)]
    # sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative
    #
    # sp$pop[year>= sc_year, energy_combined:= energy_tax + energy_lbl] #negative
    # sp$pop[is.na(energy_combined) | is.infinite(energy_combined) | energy_combined > 0, energy_combined := 0] #negative


    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_combined", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_combined") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_combined*(4.2/1000))/1.5),
                                   20.7 * ((energy_combined*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_combined", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
    )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.021, 0.082)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   1.454906 53.419114 #this is alpha, and beta
    lbl_eff_salt <- draw(9, rbeta(1, 1.454906, 53.419114))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)] #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~0.015 g for salt compensation in 2019 (18% reduction in 10 years)
    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS salt explicitly as food weight (g) x salt density (g salt per gram of HFSS food).
    #Song et al. effect is a % change in salt per 100 g (i.e. in salt density), so: change = weight x density x % effect (= salt_diet_g x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not grams of salt.
    sp$pop[, salt_density := fifelse(hfss_food_weight > 0, salt_diet_g / hfss_food_weight, 0)] #g salt per gram of HFSS food
    sp$pop[year>= sc_year, salt_combined:= ((hfss_food_weight*salt_density)*(1+tax_eff_salt)*(1+(extra_coverage_MLT*lbl_eff_salt)) - (hfss_food_weight*salt_density) + 0.015 * (1 - 0.18)^((year - 19L) / 10L))] #Updated 20261104: labelling factor no longer applied to the grain/dairy substitution term (substitution is due to the tax only; grains/dairy are not HFSS and not labelled)
    #sp$pop[year>= sc_year, salt_combined:= salt_diet_g*(((1+tax_eff_salt)*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_salt)))-1)] #negative value = reduction #Updated 20261103
    sp$pop[is.na(salt_combined) | is.infinite(salt_combined) | salt_combined > 0, salt_combined := 0] #negative
    sp$pop[, salt_density := NULL] #Updated 20261103

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change :=  (salt_combined*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_combined >0, sbp_change :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*0.60, sbp_change :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change <0, sbp_change := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change) & (sbp_change >= 0),
                                   sbp_curr_xps - sbp_change,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_combined", "salt_diet_g", "sbp_change") := NULL]

 ## Reduction in fruits and vegetables due to taxation
    #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
    #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   137.4731 549.0456 #this is alpha, and beta
    elasticity_fruit <- draw(6, rbeta(1, 137.4731, 549.0456))*(-1)  #Updated 20261103: draw()

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)))] #tax_eff_fruit is negative #Updated 20261103


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   84.65086 84.65087 #this is alpha, and beta
    elasticity_veg <- draw(7, rbeta(1, 84.65086, 84.65087))*(-1)  #Updated 20261103: draw()

    tax_eff_veg <- tax_rate*elasticity_veg #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)))] #tax_eff_veg is negative #Updated 20261103

    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_tax_lbl_sn", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_lbl_sn")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


# ################# Scenario 1: Taxation 8% and Song et al. (Consumer response only, sensitivity reduction in total intakes)
# IMPACTncd$update_primary_prevention_scn(
#   function(sp) {
#     #setDT(sp$pop)
#     # Ensure that .Random.seed exists (it is created when a random number is generated)
#     if (!exists(".Random.seed", envir = .GlobalEnv)) {
#       runif(1)
#     }
#
#     # Save the current random seed state
#     saved_seed <- .Random.seed
#
#     # Change the seed (this will change the state of the RNG)
#     set.seed(93684 + sp$mc_aggr)
#
#  ## Reduction in energy
#     # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
#     #tax_rate <- 0.08
#     #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
#     elasticity_HFSS <- runif(1, min = 0.6, max = 1.0) #per 1% increase in price #based on discussion with Chris 20250805
#     tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
#     if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative
#
#     # label Song et al.
#     #lbl_eff_energy1 <- rnorm(1, mean = -0.129, sd = (-0.08 - (-0.179))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.08, 0.129, 0.179)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   24.06169 160.75193 #this is alpha, and beta
#     lbl_eff_energy1 <- rbeta(1, 24.06169, 160.75193)*(-1)
#     if (lbl_eff_energy1 > 0) lbl_eff_energy1 <- 0 #negative
#
#     #lbl_eff_energy2 <- rnorm(1, mean = -0.064, sd = (-0.004 - (-0.125))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.004, 0.064, 0.125)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   24.06169 160.75193 #this is alpha, and beta
#     lbl_eff_energy2 <- rbeta(1, 6.360689, 88.532728)*(-1)
#     if (lbl_eff_energy2 > 0) lbl_eff_energy2 <- 0 #negative
#
#     #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
#     #prop_packaged <- 0.80
#     prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels
#
#     #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
#     #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
#     #so, 25% of the products will be labeled if the MLT become mandatory
#     #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
#     current_coverage_MLT <- 0.75
#     extra_coverage_MLT <- 1- current_coverage_MLT
#
#     sp$pop[, energy_diet:= hfss_energy]
#     sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive
#
#     #We will estimate the impact from 2026
#     sc_year <- 26L
#
#     # The combined effect can be calculated as
#     # 1 – [(1+PRi)*(1+PRj)] PR= percent reduction < 0
#     #Based on different effect sizes 12.9% vs 6.4%, we will apply 6.4% to 75% of the products currently featuring MLT as the effect of NW outperform NL
#     sp$pop[year>= sc_year, energy_combined:= energy_diet*(((1+tax_eff_energy)*(1+((current_coverage_MLT*lbl_eff_energy2) + (extra_coverage_MLT*lbl_eff_energy1))))-1)] # -1 to make it negative
#     sp$pop[is.na(energy_combined) | is.infinite(energy_combined) | energy_combined > 0, energy_combined := 0] #negative
#
#     sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation
#
#     if (Hall) {
#     #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
#     #Translating change in energy intake (kcal/day) into change in body weight using the
#     #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
#     #Each person's body weight trajectory is simulated year by year from sc_year onwards
#     #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
#     hall_update_bmi(sp$pop, energy_col = "energy_combined", sc_year = sc_year)
#     sp$pop[, c("energy_diet", "energy_combined") := NULL]
#     #Updated 20261103 (hallfast): END
#
#     } else{
#       #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
#       sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_combined*(4.2/1000))/1.5),
#                                    20.7 * ((energy_combined*(4.2/1000))/1.5))]
#
#       ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
#       sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
#       sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4
#
#       sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
#       sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)
#
#       sp$pop[year >= sc_year,
#              bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
#                                      bmi_curr_xps* relative_change_bw,
#                                      bmi_curr_xps)]
#       #sp$pop[, summary(bmi_curr_xps)]
#       sp$pop[, c("energy_diet", "energy_combined", "change_bw", "weight_assum", "relative_change_bw") := NULL]
#       # Restore the original seed state
#     }
#
#     sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
#     sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]
#
#     #adding mediated effect through SBP and Tchol #20250901
#     #https://doi.org/10.1111/jch.14661
#     #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
#     #Only applying effect to those with BMI >= 25
#
#     #https://doi.org/10.1210/clinem/dgaa673
#     #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
#     #which is -0.229mmol/L (-0.305, -0.154)
#
#     sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
#     sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
#     sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed
#
#     tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
#     tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
#     tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed
#
#     sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
#            `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
#                  tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
#     ]
#
#     sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]
#
#  ## Reduction in salt
#     sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
#                   black = as.integer(ethnicity=="black african")
#     )]
#
#     tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)
#
#     sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
#                   black = as.integer(ethnicity=="black african")
#     )]
#
#     #lbl_eff_salt1 <- rnorm(1, mean = -0.078, sd = (0.168 - (-0.012))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0, 0.078, 0.168)   # lower CI, mean, upper CI #the opposite effect was converted to 0
#     # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   4.584322 50.638770 #this is alpha, and beta
#     lbl_eff_salt1 <- rbeta(1, 4.584322, 50.638770)*(-1)
#     if (lbl_eff_salt1 > 0) lbl_eff_salt1 <- 0 #negative
#
#     #lbl_eff_salt2 <- rnorm(1, mean = -0.013, sd = (0.103 - (-0.077))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0, 0.013, 0.103)   # lower CI, mean, upper CI #the opposite effect was converted to 0
#     # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   0.6301193 26.3996677  #this is alpha, and beta
#     lbl_eff_salt2 <- rbeta(1, 0.6301193, 26.3996677)*(-1)
#     if (lbl_eff_salt2 > 0) lbl_eff_salt2 <- 0 #negative
#
#
#     sp$pop[, salt_diet_g:= hfss_salt]
#     sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive
#
#     sp$pop[year>= sc_year, salt_combined:= salt_diet_g*(((1+tax_eff_salt)*(1+((current_coverage_MLT*lbl_eff_salt2) + (extra_coverage_MLT*lbl_eff_salt1))))-1)] # -1 to make it negative
#     sp$pop[is.na(salt_combined) | is.infinite(salt_combined) | salt_combined > 0, salt_combined := 0] #negative
#
#     #Translating reduction in salt (gram) into reduction in SBP (mmHG)
#     #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
#     r_new_sbp <- runif(1, min = 0, max = 1)
#     a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
#     b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
#     b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
#     b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)
#
#     sp$pop[, sbp_change :=  (salt_combined*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
#     sp$pop[salt_combined >0, sbp_change :=  0] #no change in SBP if no reduction in salt
#     #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
#     sp$pop[salt_diet_g <=3.75*0.60, sbp_change :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
#     sp$pop[sbp_change <0, sbp_change := 0]
#
#     sp$pop[year >= sc_year,
#            sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change) & (sbp_change >= 0),
#                                    sbp_curr_xps - sbp_change,
#                                    sbp_curr_xps)]
#     #sp$pop[, summary(sbp_curr_xps)]
#     sp$pop[, c("hypertensive", "black", "salt_combined", "salt_diet_g", "sbp_change") := NULL]
#
#  ## Reduction in fruits and vegetables due to taxation
#     #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
#     #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   137.4731 549.0456 #this is alpha, and beta
#     elasticity_fruit <- rbeta(1, 137.4731, 549.0456)*(-1)
#
#     tax_eff_fruit <- tax_rate*elasticity_fruit #tax rate in %
#     if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative
#
#     sp$pop[year >= sc_year,
#            fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric
#
#
#     #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   84.65086 84.65087 #this is alpha, and beta
#     elasticity_veg <- rbeta(1, 84.65086, 84.65087)*(-1)
#
#     tax_eff_veg <- tax_rate*elasticity_veg #tax rate in %
#     if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative
#
#     sp$pop[year >= sc_year,
#            veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative
#
#     # Restore the original seed state
#     .Random.seed <- saved_seed
#   }
# )
#
# IMPACTncd$
#   run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_lbl_sn_sensi")
#
# # print("SUMMARY EXPORTS!!!")
# #
# # IMPACTncd$export_summaries(multicore = TRUE)
#
#
# ################# Scenario 2: Taxation 8% and Croker et al. (Consumer response only)
# IMPACTncd$update_primary_prevention_scn(
#   function(sp) {
#     #setDT(sp$pop)
#     # Ensure that .Random.seed exists (it is created when a random number is generated)
#     if (!exists(".Random.seed", envir = .GlobalEnv)) {
#       runif(1)
#     }
#
#     # Save the current random seed state
#     saved_seed <- .Random.seed
#
#     # Change the seed (this will change the state of the RNG)
#     set.seed(93684 + sp$mc_aggr)
#
#  ## Reduction in energy
#     # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
#     #tax_rate <- 0.08
#     #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
#     elasticity_HFSS <- runif(1, min = 0.6, max = 1.0) #per 1% increase in price #based on discussion with Chris 20250805
#     tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
#     if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative
#
#     # label Croker et al.
#     lbl_eff_energy <- rnorm(1, mean = -4.43, sd = (8.74 - 0.12)/3.92) #per 100 gram of food
#     if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative
#
#     #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
#     #prop_packaged <- 0.80
#     prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels
#
#     #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
#     #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
#     #so, 25% of the products will be labeled if the MLT become mandatory
#     #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
#     current_coverage_MLT <- 0.75
#     extra_coverage_MLT <- 1- current_coverage_MLT
#
#     sp$pop[, energy_diet:= hfss_energy]
#     sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive
#
#     #We will estimate the impact from 2026
#     sc_year <- 26L
#
#     # The combined effect can be calculated as
#     # 1 – [(1+PRi)*(1+PRj)] PR= percent reduction < 0
#     sp$pop[year>= sc_year, energy_tax:= energy_diet*tax_eff_energy]
#     # The second policy (NW) was implemented after reduction in energy intake or food grams, post energy = (1+tax_negative)*energy_diet -> the proportion = (1+tax_negative)
#     # We assumed that % reduction in energy intake is equal to % reduction in food grams
#     sp$pop[year>= sc_year, energy_lbl:= ((1+tax_eff_energy)*hfss_food_weight/100*extra_coverage_MLT)*lbl_eff_energy] #multiply with total grams consumed from HFSS
#     sp$pop[energy_lbl*(-1) > (1+tax_eff_energy)*energy_diet*extra_coverage_MLT, energy_lbl := (1+tax_eff_energy)*energy_diet*extra_coverage_MLT*(-1)] #make sure the reduction is not > energy from HFSS
#     sp$pop[year>= sc_year, energy_combined:= energy_tax + energy_lbl]
#     sp$pop[is.na(energy_combined) | is.infinite(energy_combined) | energy_combined > 0, energy_combined := 0] #negative
#
#     sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation
#
#     if (Hall) {
#     #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
#     #Translating change in energy intake (kcal/day) into change in body weight using the
#     #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
#     #Each person's body weight trajectory is simulated year by year from sc_year onwards
#     #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
#     hall_update_bmi(sp$pop, energy_col = "energy_combined", sc_year = sc_year)
#     sp$pop[, c("energy_diet", "energy_combined") := NULL]
#     #Updated 20261103 (hallfast): END
#
#     } else{
#       #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
#       sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_combined*(4.2/1000))/1.5),
#                                    20.7 * ((energy_combined*(4.2/1000))/1.5))]
#
#       ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
#       sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
#       sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4
#
#       sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
#       sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)
#
#       sp$pop[year >= sc_year,
#              bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
#                                      bmi_curr_xps* relative_change_bw,
#                                      bmi_curr_xps)]
#       #sp$pop[, summary(bmi_curr_xps)]
#       sp$pop[, c("energy_diet", "energy_tax", "energy_lbl", "energy_combined", "change_bw", "weight_assum", "relative_change_bw") := NULL]
#       # Restore the original seed state
#     }
#
#     sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
#     sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]
#
#     #adding mediated effect through SBP and Tchol #20250901
#     #https://doi.org/10.1111/jch.14661
#     #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
#     #Only applying effect to those with BMI >= 25
#
#     #https://doi.org/10.1210/clinem/dgaa673
#     #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
#     #which is -0.229mmol/L (-0.305, -0.154)
#
#     sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
#     sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
#     sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed
#
#     tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
#     tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
#     tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed
#
#     sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
#            `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
#                  tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
#     ]
#
#     sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]
#
#  ## Reduction in salt
#     sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
#                   black = as.integer(ethnicity=="black african")
#     )]
#
#     # tax effect for 8% tax rate is 6.4% (based on a study in the Scottish context, an increase 10% price of HFSS products will reduce the consumption by 8% (an average of 6%+10%/2))
#     tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)
#
#     lbl_eff_salt <- (rnorm(1, mean = -33.78, sd = (59.40 - 8.16)/3.92))/393.4 # this is sodium in mg and translate into salt (grams)
#     if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative
#
#     sp$pop[, salt_diet_g:= hfss_salt]
#     sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive
#
#     sp$pop[year>= sc_year, salt_tax:= salt_diet_g*tax_eff_salt]
#
#     sp$pop[year>= sc_year, salt_lbl:= ((1+tax_eff_salt)*hfss_food_weight/100*extra_coverage_MLT)*lbl_eff_salt] #multiply with total grams consumed from HFSS
#     sp$pop[salt_lbl*(-1) > (1+tax_eff_salt)*salt_diet_g*extra_coverage_MLT, salt_lbl := (1+tax_eff_salt)*salt_diet_g*extra_coverage_MLT*(-1)] #make sure the reduction is not > salt from HFSS
#     sp$pop[year>= sc_year, salt_combined:= salt_tax + salt_lbl]
#     sp$pop[is.na(salt_combined) | is.infinite(salt_combined) | salt_combined > 0, salt_combined := 0] #negative
#
#     #Translating reduction in salt (gram) into reduction in SBP (mmHG)
#     #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
#     r_new_sbp <- runif(1, min = 0, max = 1)
#     a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
#     b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
#     b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
#     b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)
#
#     sp$pop[, sbp_change :=  (salt_combined*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
#     sp$pop[salt_combined >0, sbp_change :=  0] #no change in SBP if no reduction in salt
#     #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
#     sp$pop[salt_diet_g <=3.75*0.60, sbp_change :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
#     sp$pop[sbp_change <0, sbp_change := 0]
#
#     sp$pop[year >= sc_year,
#            sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change) & (sbp_change >= 0),
#                                    sbp_curr_xps - sbp_change,
#                                    sbp_curr_xps)]
#     #sp$pop[, summary(sbp_curr_xps)]
#     sp$pop[, c("hypertensive", "black", "salt_tax", "salt_lbl", "salt_combined", "salt_diet_g", "sbp_change") := NULL]
#
#  ## Reduction in fruits and vegetables due to taxation
#     #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
#     #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   137.4731 549.0456 #this is alpha, and beta
#     elasticity_fruit <- rbeta(1, 137.4731, 549.0456)*(-1)
#
#     tax_eff_fruit <- tax_rate*elasticity_fruit #tax rate in %
#     if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative
#
#     sp$pop[year >= sc_year,
#            fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric
#
#
#     #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   84.65086 84.65087 #this is alpha, and beta
#     elasticity_veg <- rbeta(1, 84.65086, 84.65087)*(-1)
#
#     tax_eff_veg <- tax_rate*elasticity_veg #tax rate in %
#     if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative
#
#     sp$pop[year >= sc_year,
#            veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative
#
#     # Restore the original seed state
#     .Random.seed <- saved_seed
#   }
# )
#
# IMPACTncd$
#   run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc2_tax_lbl_cr")
#
# # print("SUMMARY EXPORTS!!!")
# #
# # IMPACTncd$export_summaries(multicore = TRUE)






### NW labelling reformulation was updated 20260317

######## CONSUMER RESPONSE AND REFORMULATION ##################################################################################################################################################
######## TAXATION ##################################################################################################################################################
################# Scenario 7: Taxation 8% (Consumer response and reformulation) [cons_ref_sc1_tax] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- draw(1, runif(1, min = 0.6, max = 1.0)) #per 1% increase in price #based on discussion with Chris 20250805  #Updated 20261103: draw()
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    # Reformulation effect
    #ref_eff_energy <- rnorm(1, mean = -0.06, sd = ((-0.06*0.85) - (-0.06*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.051, 0.06, 0.069)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   162.5123 2541.8669  #this is alpha, and beta
    ref_eff_energy <- draw(14, rbeta(1, 162.5123, 2541.8669))*(-1)  #Updated 20261103: draw()
    if (ref_eff_energy > 0) ref_eff_energy <- 0 #negative

    # We assume that 40% of the manufactures reformulate the products to avoid taxation
    prop_ref <- 0.40
    prop_tax <- 1 - prop_ref #the proportion of products with taxation

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~5.8 kcal for energy compensation in 2019 (10% increase in 10 years)
    sp$pop[year>= sc_year, energy_tax:= energy_diet*((prop_tax*tax_eff_energy) + (prop_ref*ref_eff_energy)) + (prop_tax* (5.8 * (1 + 0.10)^((year - 19L) / 10L)))]
    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_tax", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_tax") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_tax*(4.2/1000))/1.5),
                                   20.7 * ((energy_tax*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_tax", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    # Reformulation effect
    #ref_eff_salt <- rnorm(1, mean = -0.15, sd = ((-0.15*0.85) - (-0.15*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.1275, 0.15, 0.1725)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   146.4041 828.3033  #this is alpha, and beta
    ref_eff_salt <- draw(15, rbeta(1, 146.4041, 828.3033))*(-1)  #Updated 20261103: draw()
    if (ref_eff_salt > 0) ref_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)] #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~0.015 g for salt compensation in 2019 (18% reduction in 10 years)
    sp$pop[year>= sc_year, salt_tax:= salt_diet_g*((prop_tax*tax_eff_salt) + (prop_ref*ref_eff_salt)) + prop_tax*(0.015 * (1 - 0.18)^((year - 19L) / 10L))]
    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_tax >0, sbp_change_tax :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*0.60, sbp_change_tax :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_tax <0, sbp_change_tax := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_tax) & (sbp_change_tax >= 0),
                                   sbp_curr_xps - sbp_change_tax,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_tax", "salt_diet_g", "sbp_change_tax") := NULL]

 ## Reduction in fruits and vegetables due to taxation
    #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
    #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   137.4731 549.0456 #this is alpha, and beta
    elasticity_fruit <- draw(6, rbeta(1, 137.4731, 549.0456))*(-1)  #Updated 20261103: draw()

    tax_eff_fruit <- prop_tax*tax_rate*elasticity_fruit #tax_rate is a proportion (0.08 = 8%) #Updated 20261103: scaled by prop_tax, as only the non-reformulated (taxed) share of HFSS products is taxed
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)))] #tax_eff_fruit is negative #Updated 20261103


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   84.65086 84.65087 #this is alpha, and beta
    elasticity_veg <- draw(7, rbeta(1, 84.65086, 84.65087))*(-1)  #Updated 20261103: draw()

    tax_eff_veg <- prop_tax*tax_rate*elasticity_veg #tax_rate is a proportion (0.08 = 8%) #Updated 20261103: scaled by prop_tax, as only the non-reformulated (taxed) share of HFSS products is taxed
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)))] #tax_eff_veg is negative #Updated 20261103

    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_ref_sc1_tax", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_ref_sc1_tax")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


######## NUTRITION WARNING LABELLING ##################################################################################################################################################
################# Scenario 8: Song et al. (Consumer response and reformulation) [cons_ref_sc1_lbl_sn] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.009, 0.038, 0.068)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   8.893301 217.128711 #this is alpha, and beta
    lbl_eff_energy <- draw(8, rbeta(1, 8.893301, 217.128711))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    # Reformulation effect 0.9% (95% CI: [–3.1, 4.9])
    #ref_eff_energy <- rnorm(1, mean = -0.009, sd = ((0.049- (-0.031))/3.92))
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.009, 0.049)   # lower CI, mean, upper CI
    # p <- c(0, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   0.9497959 71.3291126 #this is alpha, and beta
    # ref_eff_energy <- rbeta(1, 0.9497959, 71.3291126)*(-1) #per 100 grams
    # if (ref_eff_energy > 0) ref_eff_energy <- 0 #negative

    #Updated 20260317 #the effect above 0.9% was positive (increase energy content, not negative)
    # Reformulation effect −3.9% (p=0.388), based on this study https://www.mdpi.com/2072-6643/12/8/2371,
    # Used in this study: https://www.thelancet.com/journals/lanepe/article/PIIS2666-7762(25)00298-4/fulltext
    # Calculating potential CIs
    # library(neatStats)
    # ci_from_p(diff = 0.039, p_value = 0.388)
    # ci_lower    ci_upper
    # -0.04830619  0.12630619
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.039, 0.126)   # lower CI, mean, upper CI
    # p <- c(0, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   1.932238 40.036484 #this is alpha, and beta
    ref_eff_energy <- draw(14, rbeta(1, 1.932238, 40.036484))*(-1)  #Updated 20261103: draw()
    if (ref_eff_energy > 0) ref_eff_energy <- 0 #negative

    # We assume that 20% of the manufactures reformulate the products to avoid labelling
    prop_ref <- 0.20
    prop_lbl <- 1 - prop_ref #the proportion of products with labelling

    #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
    #prop_packaged <- 0.80
    prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels

    #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the MLT become mandatory
    #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
    current_coverage_MLT <- 0.75
    extra_coverage_MLT <- 1- current_coverage_MLT

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)*prop_packaged] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS energy explicitly as food weight (g) x energy density (kcal per gram of HFSS food).
    #Song et al. effect is a % change in energy per 100 g (i.e. in energy density), so: change = weight x density x % effect (= energy_diet x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not kcal.
    sp$pop[, energy_density := fifelse(hfss_food_weight > 0, energy_diet / hfss_food_weight, 0)] #kcal per gram of HFSS food
    sp$pop[year>= sc_year, energy_lbl:= (hfss_food_weight*energy_density)*(((extra_coverage_MLT*lbl_eff_energy)*prop_lbl) + (ref_eff_energy*prop_ref))]
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative
    sp$pop[, energy_density := NULL] #Updated 20261103

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_lbl", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_lbl") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_lbl*(4.2/1000))/1.5),
                                   20.7 * ((energy_lbl*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_lbl", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.021, 0.082)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   1.454906 53.419114 #this is alpha, and beta
    lbl_eff_salt <- draw(9, rbeta(1, 1.454906, 53.419114))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative


    # Reformulation effect 8.9% (95% CI: [–17.3%, –0.6%]
    #ref_eff_salt <- rnorm(1, mean = -0.089, sd = ((-0.006- (-0.173))/3.92))
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.006, 0.089, 0.173)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   6.184802 60.264501 #this is alpha, and beta
    # ref_eff_salt <- rbeta(1, 6.184802, 60.264501)*(-1)
    # if (ref_eff_salt > 0) ref_eff_salt <- 0 #negative

    #Updated 20260317
    # Reformulation effect −9.2% (p = 0.118), based on this study https://www.mdpi.com/2072-6643/12/8/2371,
    # Calculating potential CIs
    # library(neatStats)
    # ci_from_p(diff = 0.092, p_value = 0.118)
    # ci_lower    ci_upper
    # -0.02336739  0.20736739
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.092, 0.207)   # lower CI, mean, upper CI
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   3.940064 35.985975 #this is alpha, and beta
    ref_eff_salt <- draw(15, rbeta(1, 3.940064, 35.985975))*(-1)  #Updated 20261103: draw()
    if (ref_eff_salt > 0) ref_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)*prop_packaged] #see above #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS salt explicitly as food weight (g) x salt density (g salt per gram of HFSS food).
    #Song et al. effect is a % change in salt per 100 g (i.e. in salt density), so: change = weight x density x % effect (= salt_diet_g x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not grams of salt.
    sp$pop[, salt_density := fifelse(hfss_food_weight > 0, salt_diet_g / hfss_food_weight, 0)] #g salt per gram of HFSS food
    sp$pop[year>= sc_year, salt_lbl:= (hfss_food_weight*salt_density)*(((extra_coverage_MLT*lbl_eff_salt)*prop_lbl) + (ref_eff_salt*prop_ref))]
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative
    sp$pop[, salt_density := NULL] #Updated 20261103

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_lbl >0, sbp_change_lbl :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*prop_packaged*0.60, sbp_change_lbl :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_lbl <0, sbp_change_lbl := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_lbl) & (sbp_change_lbl >= 0),
                                   sbp_curr_xps - sbp_change_lbl,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_lbl", "salt_diet_g", "sbp_change_lbl") := NULL]
    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_ref_sc1_lbl_sn", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_ref_sc1_lbl_sn")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


# ################# Scenario 1: Song et al. (Consumer response and reformulation, sensitivity reduction in total intakes)
# IMPACTncd$update_primary_prevention_scn(
#   function(sp) {
#     #setDT(sp$pop)
#     # Ensure that .Random.seed exists (it is created when a random number is generated)
#     if (!exists(".Random.seed", envir = .GlobalEnv)) {
#       runif(1)
#     }
#
#     # Save the current random seed state
#     saved_seed <- .Random.seed
#
#     # Change the seed (this will change the state of the RNG)
#     set.seed(93684 + sp$mc_aggr)
#
#  ## Reduction in energy
#     #lbl_eff_energy1 <- rnorm(1, mean = -0.129, sd = (-0.08 - (-0.179))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.08, 0.129, 0.179)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   24.06169 160.75193 #this is alpha, and beta
#     lbl_eff_energy1 <- rbeta(1, 24.06169, 160.75193)*(-1)
#     if (lbl_eff_energy1 > 0) lbl_eff_energy1 <- 0 #negative
#
#     #lbl_eff_energy2 <- rnorm(1, mean = -0.064, sd = (-0.004 - (-0.125))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.004, 0.064, 0.125)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   24.06169 160.75193 #this is alpha, and beta
#     lbl_eff_energy2 <- rbeta(1, 6.360689, 88.532728)*(-1)
#     if (lbl_eff_energy2 > 0) lbl_eff_energy2 <- 0 #negative
#
#     # Reformulation effect 0.9% (95% CI: [–3.1, 4.9])
#     #ref_eff_energy <- rnorm(1, mean = -0.009, sd = ((0.049- (-0.031))/3.92))
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0, 0.009, 0.049)   # lower CI, mean, upper CI
#     # p <- c(0, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   0.9497959 71.3291126 #this is alpha, and beta
#     # ref_eff_energy <- rbeta(1, 0.9497959, 71.3291126)*(-1) #per 100 grams
#     # if (ref_eff_energy > 0) ref_eff_energy <- 0 #negative
#
#     #Updated 20260317 #the effect above 0.9% was positive (increase energy content, not negative)
#     # Reformulation effect −3.9% (p=0.388), based on this study https://www.mdpi.com/2072-6643/12/8/2371,
#     # Used in this study: https://www.thelancet.com/journals/lanepe/article/PIIS2666-7762(25)00298-4/fulltext
#     # Calculating potential CIs
#     # library(neatStats)
#     # ci_from_p(diff = 0.039, p_value = 0.388)
#     # ci_lower    ci_upper
#     # -0.04830619  0.12630619
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0, 0.039, 0.126)   # lower CI, mean, upper CI
#     # p <- c(0, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   1.932238 40.036484 #this is alpha, and beta
#     ref_eff_energy <- rbeta(1, 1.932238, 40.036484)*(-1)
#     if (ref_eff_energy > 0) ref_eff_energy <- 0 #negative
#
#     # We assume that 20% of the manufactures reformulate the products to avoid labelling
#     prop_ref <- 0.20
#     prop_lbl <- 1 - prop_ref #the proportion of products with labelling
#
#     #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
#     #prop_packaged <- 0.80
#     prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels
#
#     #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
#     #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
#     #so, 25% of the products will be labeled if the MLT become mandatory
#     #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
#     current_coverage_MLT <- 0.75
#     extra_coverage_MLT <- 1- current_coverage_MLT
#
#     sp$pop[, energy_diet:= hfss_energy*prop_packaged]
#     sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive
#
#     #We will estimate the impact from 2026
#     sc_year <- 26L
#     #Based on different effect sizes 12.9% vs 6.4%, we will apply 6.4% to 75% of the products currently featuring MLT as the effect of NW outperform NL
#     sp$pop[year>= sc_year, energy_lbl:= energy_diet*((((current_coverage_MLT*lbl_eff_energy2) + (extra_coverage_MLT*lbl_eff_energy1))*prop_lbl) + (ref_eff_energy*prop_ref))]
#     sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative
#
#     sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation
#
#     if (Hall) {
#     #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
#     #Translating change in energy intake (kcal/day) into change in body weight using the
#     #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
#     #Each person's body weight trajectory is simulated year by year from sc_year onwards
#     #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
#     hall_update_bmi(sp$pop, energy_col = "energy_lbl", sc_year = sc_year)
#     sp$pop[, c("energy_diet", "energy_lbl") := NULL]
#     #Updated 20261103 (hallfast): END
#
#     } else{
#       #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
#       sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_lbl*(4.2/1000))/1.5),
#                                    20.7 * ((energy_lbl*(4.2/1000))/1.5))]
#
#       ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
#       sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
#       sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4
#
#       sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
#       sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)
#
#       sp$pop[year >= sc_year,
#              bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
#                                      bmi_curr_xps* relative_change_bw,
#                                      bmi_curr_xps)]
#       #sp$pop[, summary(bmi_curr_xps)]
#       sp$pop[, c("energy_diet", "energy_lbl", "change_bw", "weight_assum", "relative_change_bw") := NULL]
#       # Restore the original seed state
#     }
#
#     sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
#     sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]
#
#     #adding mediated effect through SBP and Tchol #20250901
#     #https://doi.org/10.1111/jch.14661
#     #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
#     #Only applying effect to those with BMI >= 25
#
#     #https://doi.org/10.1210/clinem/dgaa673
#     #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
#     #which is -0.229mmol/L (-0.305, -0.154)
#
#     sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
#     sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
#     sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed
#
#     tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
#     tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
#     tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed
#
#     sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
#            `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
#                  tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
#     ]
#
#     sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]
#
#  ## Reduction in salt
#     sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
#                   black = as.integer(ethnicity=="black african")
#                   )]
#
#     #lbl_eff_salt1 <- rnorm(1, mean = -0.078, sd = (0.168 - (-0.012))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0, 0.078, 0.168)   # lower CI, mean, upper CI #the opposite effect was converted to 0
#     # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   4.584322 50.638770 #this is alpha, and beta
#     lbl_eff_salt1 <- rbeta(1, 4.584322, 50.638770)*(-1)
#     if (lbl_eff_salt1 > 0) lbl_eff_salt1 <- 0 #negative
#
#     #lbl_eff_salt2 <- rnorm(1, mean = -0.013, sd = (0.103 - (-0.077))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0, 0.013, 0.103)   # lower CI, mean, upper CI #the opposite effect was converted to 0
#     # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   0.6301193 26.3996677  #this is alpha, and beta
#     lbl_eff_salt2 <- rbeta(1, 0.6301193, 26.3996677)*(-1)
#     if (lbl_eff_salt2 > 0) lbl_eff_salt2 <- 0 #negative
#
#     # Reformulation effect 8.9% (95% CI: [–17.3%, –0.6%]
#     #ref_eff_salt <- rnorm(1, mean = -0.089, sd = ((-0.006- (-0.173))/3.92))
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.006, 0.089, 0.173)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   6.184802 60.264501 #this is alpha, and beta
#     # ref_eff_salt <- rbeta(1, 6.184802, 60.264501)*(-1)
#     # if (ref_eff_salt > 0) ref_eff_salt <- 0 #negative
#
#     #Updated 20260317
#     # Reformulation effect −9.2% (p = 0.118), based on this study https://www.mdpi.com/2072-6643/12/8/2371,
#     # Calculating potential CIs
#     # library(neatStats)
#     # ci_from_p(diff = 0.092, p_value = 0.118)
#     # ci_lower    ci_upper
#     # -0.02336739  0.20736739
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0, 0.092, 0.207)   # lower CI, mean, upper CI
#     # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   3.940064 35.985975 #this is alpha, and beta
#     ref_eff_salt <- rbeta(1, 3.940064, 35.985975)*(-1)
#     if (ref_eff_salt > 0) ref_eff_salt <- 0 #negative
#
#     sp$pop[, salt_diet_g:= hfss_salt*prop_packaged] #see above
#     sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive
#
#     sp$pop[year>= sc_year, salt_lbl:= salt_diet_g*((((current_coverage_MLT*lbl_eff_salt2) + (extra_coverage_MLT*lbl_eff_salt1))*prop_lbl) + (ref_eff_salt*prop_ref))]
#     sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative
#
#     #Translating reduction in salt (gram) into reduction in SBP (mmHG)
#     #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
#     r_new_sbp <- runif(1, min = 0, max = 1)
#     a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
#     b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
#     b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
#     b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)
#
#     sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
#     sp$pop[salt_lbl >0, sbp_change_lbl :=  0] #no change in SBP if no reduction in salt
#     #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
#     sp$pop[salt_diet_g <=3.75*prop_packaged*0.60, sbp_change_lbl :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
#     sp$pop[sbp_change_lbl <0, sbp_change_lbl := 0]
#
#     sp$pop[year >= sc_year,
#            sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_lbl) & (sbp_change_lbl >= 0),
#                                    sbp_curr_xps - sbp_change_lbl,
#                                    sbp_curr_xps)]
#     #sp$pop[, summary(sbp_curr_xps)]
#     sp$pop[, c("hypertensive", "black", "salt_lbl", "salt_diet_g", "sbp_change_lbl") := NULL]
#     # Restore the original seed state
#     .Random.seed <- saved_seed
#   }
# )
#
# IMPACTncd$
#   run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_ref_sc1_lbl_sn_sensi")
#
# # print("SUMMARY EXPORTS!!!")
# #
# # IMPACTncd$export_summaries(multicore = TRUE)


######## COMBINED EFFECTS OF TAXATION AND NUTRITION WARNING LABELLING ##################################################################################################################################################
################# Scenario 9: Taxation 8% and Song et al. (Consumer response and reformulation) [cons_ref_sc1_tax_lbl_sn] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- draw(1, runif(1, min = 0.6, max = 1.0)) #per 1% increase in price #based on discussion with Chris 20250805  #Updated 20261103: draw()
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    # label Song et al.
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.009, 0.038, 0.068)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   8.893301 217.128711 #this is alpha, and beta
    lbl_eff_energy <- draw(8, rbeta(1, 8.893301, 217.128711))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    # Reformulation effect (We used from the taxation only as the coverage is higher)
    #ref_eff_energy <- rnorm(1, mean = -0.06, sd = ((-0.06*0.85) - (-0.06*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.051, 0.06, 0.069)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   162.5123 2541.8669  #this is alpha, and beta
    ref_eff_energy <- draw(14, rbeta(1, 162.5123, 2541.8669))*(-1)  #Updated 20261103: draw()
    if (ref_eff_energy > 0) ref_eff_energy <- 0 #negative

    # We assume that 40% of the manufactures reformulate the products to avoid taxation
    prop_ref <- 0.40
    prop_tax <- 1 - prop_ref #the proportion of products with taxation

    #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
    #prop_packaged <- 0.80
    prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels

    #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the MLT become mandatory
    #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
    current_coverage_MLT <- 0.75
    extra_coverage_MLT <- 1- current_coverage_MLT

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~5.8 kcal for energy compensation in 2019 (10% increase in 10 years)
    # The combined effect can be calculated as
    # 1 – [(1+PRi)*(1+PRj)] PR= percent reduction < 0
    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS energy explicitly as food weight (g) x energy density (kcal per gram of HFSS food).
    #Song et al. effect is a % change in energy per 100 g (i.e. in energy density), so: change = weight x density x % effect (= energy_diet x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not kcal.
    sp$pop[, energy_density := fifelse(hfss_food_weight > 0, energy_diet / hfss_food_weight, 0)] #kcal per gram of HFSS food
    sp$pop[year>= sc_year, energy_combined:= prop_tax*(((hfss_food_weight*energy_density)*(1+tax_eff_energy)*(1+(extra_coverage_MLT*lbl_eff_energy)) - (hfss_food_weight*energy_density) + 5.8 * (1 + 0.10)^((year - 19L) / 10L))) #Updated 20261104: labelling factor no longer applied to the grain/dairy substitution term (substitution is due to the tax only; grains/dairy are not HFSS and not labelled)
                                             + ((hfss_food_weight*energy_density)*prop_ref*ref_eff_energy)] #negative value = reduction #Updated 20261103
    sp$pop[is.na(energy_combined) | is.infinite(energy_combined) | energy_combined > 0, energy_combined := 0] #negative
    sp$pop[, energy_density := NULL] #Updated 20261103

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_combined", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_combined") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_combined*(4.2/1000))/1.5),
                                   20.7 * ((energy_combined*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_combined", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
    )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.021, 0.082)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   1.454906 53.419114 #this is alpha, and beta
    lbl_eff_salt <- draw(9, rbeta(1, 1.454906, 53.419114))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    # Reformulation effect
    #ref_eff_salt <- rnorm(1, mean = -0.15, sd = ((-0.15*0.85) - (-0.15*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.1275, 0.15, 0.1725)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   146.4041 828.3033  #this is alpha, and beta
    ref_eff_salt <- draw(15, rbeta(1, 146.4041, 828.3033))*(-1)  #Updated 20261103: draw()
    if (ref_eff_salt > 0) ref_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)] #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~0.015 g for salt compensation in 2019 (18% reduction in 10 years)
    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS salt explicitly as food weight (g) x salt density (g salt per gram of HFSS food).
    #Song et al. effect is a % change in salt per 100 g (i.e. in salt density), so: change = weight x density x % effect (= salt_diet_g x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not grams of salt.
    sp$pop[, salt_density := fifelse(hfss_food_weight > 0, salt_diet_g / hfss_food_weight, 0)] #g salt per gram of HFSS food
    sp$pop[year>= sc_year, salt_combined:= prop_tax*(((hfss_food_weight*salt_density)*(1+tax_eff_salt)*(1+(extra_coverage_MLT*lbl_eff_salt)) - (hfss_food_weight*salt_density) + 0.015 * (1 - 0.18)^((year - 19L) / 10L))) #Updated 20261104: labelling factor no longer applied to the grain/dairy substitution term (substitution is due to the tax only; grains/dairy are not HFSS and not labelled)
                                           + ((hfss_food_weight*salt_density)*prop_ref*ref_eff_salt)] #negative value = reduction #Updated 20261103

    sp$pop[is.na(salt_combined) | is.infinite(salt_combined) | salt_combined > 0, salt_combined := 0] #negative
    sp$pop[, salt_density := NULL] #Updated 20261103

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change :=  (salt_combined*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_combined >0, sbp_change :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*0.60, sbp_change :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change <0, sbp_change := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change) & (sbp_change >= 0),
                                   sbp_curr_xps - sbp_change,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_combined", "salt_diet_g", "sbp_change") := NULL]

 ## Reduction in fruits and vegetables due to taxation
    #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
    #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   137.4731 549.0456 #this is alpha, and beta
    elasticity_fruit <- draw(6, rbeta(1, 137.4731, 549.0456))*(-1)  #Updated 20261103: draw()

    tax_eff_fruit <- prop_tax*tax_rate*elasticity_fruit #tax_rate is a proportion (0.08 = 8%) #Updated 20261103: scaled by prop_tax, as only the non-reformulated (taxed) share of HFSS products is taxed
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)))] #tax_eff_fruit is negative #Updated 20261103


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   84.65086 84.65087 #this is alpha, and beta
    elasticity_veg <- draw(7, rbeta(1, 84.65086, 84.65087))*(-1)  #Updated 20261103: draw()

    tax_eff_veg <- prop_tax*tax_rate*elasticity_veg #tax_rate is a proportion (0.08 = 8%) #Updated 20261103: scaled by prop_tax, as only the non-reformulated (taxed) share of HFSS products is taxed
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)))] #tax_eff_veg is negative #Updated 20261103

    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_ref_sc1_tax_lbl_sn", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_ref_sc1_tax_lbl_sn")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


# ################# Scenario 1: Taxation 8% and Song et al. (Consumer response only, sensitivity reduction in total intakes)
# IMPACTncd$update_primary_prevention_scn(
#   function(sp) {
#     #setDT(sp$pop)
#     # Ensure that .Random.seed exists (it is created when a random number is generated)
#     if (!exists(".Random.seed", envir = .GlobalEnv)) {
#       runif(1)
#     }
#
#     # Save the current random seed state
#     saved_seed <- .Random.seed
#
#     # Change the seed (this will change the state of the RNG)
#     set.seed(93684 + sp$mc_aggr)
#
#  ## Reduction in energy
#     # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
#     #tax_rate <- 0.08
#     #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
#     elasticity_HFSS <- runif(1, min = 0.6, max = 1.0) #per 1% increase in price #based on discussion with Chris 20250805
#     tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
#     if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative
#
#     # label Song et al.
#     #lbl_eff_energy1 <- rnorm(1, mean = -0.129, sd = (-0.08 - (-0.179))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.08, 0.129, 0.179)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   24.06169 160.75193 #this is alpha, and beta
#     lbl_eff_energy1 <- rbeta(1, 24.06169, 160.75193)*(-1)
#     if (lbl_eff_energy1 > 0) lbl_eff_energy1 <- 0 #negative
#
#     #lbl_eff_energy2 <- rnorm(1, mean = -0.064, sd = (-0.004 - (-0.125))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.004, 0.064, 0.125)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   24.06169 160.75193 #this is alpha, and beta
#     lbl_eff_energy2 <- rbeta(1, 6.360689, 88.532728)*(-1)
#     if (lbl_eff_energy2 > 0) lbl_eff_energy2 <- 0 #negative
#
#     # Reformulation effect (We used from the taxation only as the coverage is higher)
#     #ref_eff_energy <- rnorm(1, mean = -0.06, sd = ((-0.06*0.85) - (-0.06*1.15))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.051, 0.06, 0.069)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   162.5123 2541.8669  #this is alpha, and beta
#     ref_eff_energy <- rbeta(1, 162.5123, 2541.8669)*(-1)
#     if (ref_eff_energy > 0) ref_eff_energy <- 0 #negative
#
#     # We assume that 40% of the manufactures reformulate the products to avoid taxation
#     prop_ref <- 0.40
#     prop_tax <- 1 - prop_ref #the proportion of products with taxation
#
#     #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
#     #prop_packaged <- 0.80
#     prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels
#
#     #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
#     #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
#     #so, 25% of the products will be labeled if the MLT become mandatory
#     #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
#     current_coverage_MLT <- 0.75
#     extra_coverage_MLT <- 1- current_coverage_MLT
#
#     sp$pop[, energy_diet:= hfss_energy]
#     sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive
#
#     #We will estimate the impact from 2026
#     sc_year <- 26L
#
#     # The combined effect can be calculated as
#     # 1 – [(1+PRi)*(1+PRj)] PR= percent reduction < 0
#     #Based on different effect sizes 12.9% vs 6.4%, we will apply 6.4% to 75% of the products currently featuring MLT as the effect of NW outperform NL
#     sp$pop[year>= sc_year, energy_combined:= energy_diet*prop_tax*(((1+tax_eff_energy)*(1+((current_coverage_MLT*lbl_eff_energy2) + (extra_coverage_MLT*lbl_eff_energy1))))-1)
#                                              + (energy_diet*prop_ref*ref_eff_energy)] # -1 to make it negative
#     sp$pop[is.na(energy_combined) | is.infinite(energy_combined) | energy_combined > 0, energy_combined := 0] #negative
#
#     sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation
#
#     if (Hall) {
#     #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
#     #Translating change in energy intake (kcal/day) into change in body weight using the
#     #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
#     #Each person's body weight trajectory is simulated year by year from sc_year onwards
#     #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
#     hall_update_bmi(sp$pop, energy_col = "energy_combined", sc_year = sc_year)
#     sp$pop[, c("energy_diet", "energy_combined") := NULL]
#     #Updated 20261103 (hallfast): END
#
#     } else{
#       #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
#       sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_combined*(4.2/1000))/1.5),
#                                    20.7 * ((energy_combined*(4.2/1000))/1.5))]
#
#       ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
#       sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
#       sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4
#
#       sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
#       sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)
#
#       sp$pop[year >= sc_year,
#              bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
#                                      bmi_curr_xps* relative_change_bw,
#                                      bmi_curr_xps)]
#       #sp$pop[, summary(bmi_curr_xps)]
#       sp$pop[, c("energy_diet", "energy_combined", "change_bw", "weight_assum", "relative_change_bw") := NULL]
#       # Restore the original seed state
#     }
#
#     sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
#     sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]
#
#     #adding mediated effect through SBP and Tchol #20250901
#     #https://doi.org/10.1111/jch.14661
#     #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
#     #Only applying effect to those with BMI >= 25
#
#     #https://doi.org/10.1210/clinem/dgaa673
#     #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
#     #which is -0.229mmol/L (-0.305, -0.154)
#
#     sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
#     sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
#     sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed
#
#     tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
#     tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
#     tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed
#
#     sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
#            `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
#                  tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
#     ]
#
#     sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]
#
#  ## Reduction in salt
#     sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
#                   black = as.integer(ethnicity=="black african")
#     )]
#
#     tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)
#
#     #lbl_eff_salt1 <- rnorm(1, mean = -0.078, sd = (0.168 - (-0.012))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0, 0.078, 0.168)   # lower CI, mean, upper CI #the opposite effect was converted to 0
#     # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   4.584322 50.638770 #this is alpha, and beta
#     lbl_eff_salt1 <- rbeta(1, 4.584322, 50.638770)*(-1)
#     if (lbl_eff_salt1 > 0) lbl_eff_salt1 <- 0 #negative
#
#     #lbl_eff_salt2 <- rnorm(1, mean = -0.013, sd = (0.103 - (-0.077))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0, 0.013, 0.103)   # lower CI, mean, upper CI #the opposite effect was converted to 0
#     # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   0.6301193 26.3996677  #this is alpha, and beta
#     lbl_eff_salt2 <- rbeta(1, 0.6301193, 26.3996677)*(-1)
#     if (lbl_eff_salt2 > 0) lbl_eff_salt2 <- 0 #negative
#
#     # Reformulation effect
#     #ref_eff_salt <- rnorm(1, mean = -0.15, sd = ((-0.15*0.85) - (-0.15*1.15))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.1275, 0.15, 0.1725)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   146.4041 828.3033  #this is alpha, and beta
#     ref_eff_salt <- rbeta(1, 146.4041, 828.3033)*(-1)
#     if (ref_eff_salt > 0) ref_eff_salt <- 0 #negative
#
#     sp$pop[, salt_diet_g:= hfss_salt]
#     sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive
#
#     sp$pop[year>= sc_year, salt_combined:= salt_diet_g*prop_tax*(((1+tax_eff_salt)*(1+((current_coverage_MLT*lbl_eff_salt2) + (extra_coverage_MLT*lbl_eff_salt1))))-1)
#                                            + (salt_diet_g*prop_ref*ref_eff_salt)] # -1 to make it negative
#     sp$pop[is.na(salt_combined) | is.infinite(salt_combined) | salt_combined > 0, salt_combined := 0] #negative
#
#     #Translating reduction in salt (gram) into reduction in SBP (mmHG)
#     #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
#     r_new_sbp <- runif(1, min = 0, max = 1)
#     a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
#     b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
#     b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
#     b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)
#
#     sp$pop[, sbp_change :=  (salt_combined*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
#     sp$pop[salt_combined >0, sbp_change :=  0] #no change in SBP if no reduction in salt
#     #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
#     sp$pop[salt_diet_g <=3.75*0.60, sbp_change :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
#     sp$pop[sbp_change <0, sbp_change := 0]
#
#     sp$pop[year >= sc_year,
#            sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change) & (sbp_change >= 0),
#                                    sbp_curr_xps - sbp_change,
#                                    sbp_curr_xps)]
#     #sp$pop[, summary(sbp_curr_xps)]
#     sp$pop[, c("hypertensive", "black", "salt_combined", "salt_diet_g", "sbp_change") := NULL]
#
#  ## Reduction in fruits and vegetables due to taxation
#     #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
#     #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   137.4731 549.0456 #this is alpha, and beta
#     elasticity_fruit <- rbeta(1, 137.4731, 549.0456)*(-1)
#
#     tax_eff_fruit <- tax_rate*elasticity_fruit #tax rate in %
#     if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative
#
#     sp$pop[year >= sc_year,
#            fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric
#
#
#     #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
#     # using Betta distribution for percentage, as suggested by Chris:
#     # install.packages("rriskDistributions")
#     # library(rriskDistributions)
#     # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
#     # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
#     #
#     # params <- get.beta.par(q = q, p = p)
#     # params
#     #
#     # $par
#     # [1]   84.65086 84.65087 #this is alpha, and beta
#     elasticity_veg <- rbeta(1, 84.65086, 84.65087)*(-1)
#
#     tax_eff_veg <- tax_rate*elasticity_veg #tax rate in %
#     if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative
#
#     sp$pop[year >= sc_year,
#            veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative
#
#     # Restore the original seed state
#     .Random.seed <- saved_seed
#   }
# )
#
# IMPACTncd$
#   run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_ref_sc1_tax_lbl_sn_sensi")
#
# # print("SUMMARY EXPORTS!!!")
# #
# # IMPACTncd$export_summaries(multicore = TRUE)





######## SENSITIVITY ANALYSES WITH COMPENSATION (FOR CONSUMER RESPONSE ONLY!!!) ##################################################################################################################################################
######## 11% COMPENSATION ON ENERGY AND SALT!!! #Updated 20261103: header corrected (was "ON ENERGY ONLY"), the compensation is applied to both energy and salt https://doi.org/10.1186/s12966-022-01287-z

comps <- 0.11

######## TAXATION ##################################################################################################################################################
################# Scenario 10: Taxation 8% (Consumer response only) [cons_sc1_tax_comps] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- draw(1, runif(1, min = 0.6, max = 1.0)) #per 1% increase in price #based on discussion with Chris 20250805  #Updated 20261103: draw()
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    sp$pop[year>= sc_year, energy_tax:= energy_diet*tax_eff_energy*(1-comps) + 5.8 * (1 + 0.10)^((year - 19L) / 10L)] #the effect was reduced for compensation #Updated 20261103: (1-comps) applied to the policy effect only, not to the grain/dairy compensation term
    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_tax", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_tax") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_tax*(4.2/1000))/1.5),
                                   20.7 * ((energy_tax*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_tax", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
           ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)] #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_tax:= salt_diet_g*tax_eff_salt*(1-comps) + 0.015 * (1 - 0.18)^((year - 19L) / 10L)] #the effect was reduced for compensation (energy and salt) #Updated 20261103: (1-comps) applied to the policy effect only, not to the grain/dairy compensation term
    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_tax >0, sbp_change_tax :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*0.60, sbp_change_tax :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_tax <0, sbp_change_tax := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_tax) & (sbp_change_tax >= 0),
                                   sbp_curr_xps - sbp_change_tax,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_tax", "salt_diet_g", "sbp_change_tax") := NULL]

 ## Reduction in fruits and vegetables due to taxation
    #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
    #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   137.4731 549.0456 #this is alpha, and beta
    elasticity_fruit <- draw(6, rbeta(1, 137.4731, 549.0456))*(-1)  #Updated 20261103: draw()

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)))] #tax_eff_fruit is negative #Updated 20261103


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   84.65086 84.65087 #this is alpha, and beta
    elasticity_veg <- draw(7, rbeta(1, 84.65086, 84.65087))*(-1)  #Updated 20261103: draw()

    tax_eff_veg <- tax_rate*elasticity_veg #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)))] #tax_eff_veg is negative #Updated 20261103

    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_tax_comps", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_comps")) #Updated 20261104: timing
#
# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


######## NUTRITION WARNING LABELLING ##################################################################################################################################################
################# Scenario 11: Song et al. (Consumer response only) [cons_sc1_lbl_sn_comps] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.009, 0.038, 0.068)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   8.893301 217.128711 #this is alpha, and beta
    lbl_eff_energy <- draw(8, rbeta(1, 8.893301, 217.128711))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
    #prop_packaged <- 0.80
    prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels

    #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the MLT become mandatory
    #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
    current_coverage_MLT <- 0.75
    extra_coverage_MLT <- 1- current_coverage_MLT

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)*prop_packaged] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS energy explicitly as food weight (g) x energy density (kcal per gram of HFSS food).
    #Song et al. effect is a % change in energy per 100 g (i.e. in energy density), so: change = weight x density x % effect (= energy_diet x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not kcal.
    sp$pop[, energy_density := fifelse(hfss_food_weight > 0, energy_diet / hfss_food_weight, 0)] #kcal per gram of HFSS food
    sp$pop[year>= sc_year, energy_lbl:= (hfss_food_weight*energy_density)*(extra_coverage_MLT*lbl_eff_energy)*(1-comps)] #the effect was reduced for compensation
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative
    sp$pop[, energy_density := NULL] #Updated 20261103

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_lbl", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_lbl") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_lbl*(4.2/1000))/1.5),
                                   20.7 * ((energy_lbl*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_lbl", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.021, 0.082)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   1.454906 53.419114 #this is alpha, and beta
    lbl_eff_salt <- draw(9, rbeta(1, 1.454906, 53.419114))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)*prop_packaged] #see above #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS salt explicitly as food weight (g) x salt density (g salt per gram of HFSS food).
    #Song et al. effect is a % change in salt per 100 g (i.e. in salt density), so: change = weight x density x % effect (= salt_diet_g x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not grams of salt.
    sp$pop[, salt_density := fifelse(hfss_food_weight > 0, salt_diet_g / hfss_food_weight, 0)] #g salt per gram of HFSS food
    sp$pop[year>= sc_year, salt_lbl:= (hfss_food_weight*salt_density)*(extra_coverage_MLT*lbl_eff_salt)*(1-comps)] #the effect was reduced for compensation (energy and salt)
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative
    sp$pop[, salt_density := NULL] #Updated 20261103

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_lbl >0, sbp_change_lbl :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*prop_packaged*0.60, sbp_change_lbl :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_lbl <0, sbp_change_lbl := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_lbl) & (sbp_change_lbl >= 0),
                                   sbp_curr_xps - sbp_change_lbl,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_lbl", "salt_diet_g", "sbp_change_lbl") := NULL]
    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_lbl_sn_comps", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_lbl_sn_comps")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


######## COMBINED EFFECTS OF TAXATION AND NUTRITION WARNING LABELLING ##################################################################################################################################################
################# Scenario 12: Taxation 8% and Song et al. (Consumer response only) [cons_sc1_tax_lbl_sn_comps] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- draw(1, runif(1, min = 0.6, max = 1.0)) #per 1% increase in price #based on discussion with Chris 20250805  #Updated 20261103: draw()
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    # label Song et al.
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.009, 0.038, 0.068)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   8.893301 217.128711 #this is alpha, and beta
    lbl_eff_energy <- draw(8, rbeta(1, 8.893301, 217.128711))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
    #prop_packaged <- 0.80
    prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels

    #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the MLT become mandatory
    #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
    current_coverage_MLT <- 0.75
    extra_coverage_MLT <- 1- current_coverage_MLT

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    # The combined effect can be calculated as
    # 1 – [(1+PRi)*(1+PRj)] PR= percent reduction < 0
    #sp$pop[year>= sc_year, energy_combined:= (energy_diet*(1+tax_eff_energy)*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy))) - energy_diet]
    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS energy explicitly as food weight (g) x energy density (kcal per gram of HFSS food).
    #Song et al. effect is a % change in energy per 100 g (i.e. in energy density), so: change = weight x density x % effect (= energy_diet x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not kcal.
    sp$pop[, energy_density := fifelse(hfss_food_weight > 0, energy_diet / hfss_food_weight, 0)] #kcal per gram of HFSS food
    sp$pop[year>= sc_year, energy_combined:= (((hfss_food_weight*energy_density)*(1+tax_eff_energy)*(1+(extra_coverage_MLT*lbl_eff_energy)))-(hfss_food_weight*energy_density))*(1-comps) + 5.8 * (1 + 0.10)^((year - 19L) / 10L)] #negative value = reduction # the same as above #the effect was reduced for compensation #Updated 20261103: (1-comps) applied to the policy effect only, not to the grain/dairy compensation term #Updated 20261104: labelling factor no longer applied to the grain/dairy substitution term (substitution is due to the tax only; grains/dairy are not HFSS and not labelled)
    sp$pop[is.na(energy_combined) | is.infinite(energy_combined) | energy_combined > 0, energy_combined := 0] #negative
    sp$pop[, energy_density := NULL] #Updated 20261103

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_combined", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_combined") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_combined*(4.2/1000))/1.5),
                                   20.7 * ((energy_combined*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_combined", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
    )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.021, 0.082)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   1.454906 53.419114 #this is alpha, and beta
    lbl_eff_salt <- draw(9, rbeta(1, 1.454906, 53.419114))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)] #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS salt explicitly as food weight (g) x salt density (g salt per gram of HFSS food).
    #Song et al. effect is a % change in salt per 100 g (i.e. in salt density), so: change = weight x density x % effect (= salt_diet_g x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not grams of salt.
    sp$pop[, salt_density := fifelse(hfss_food_weight > 0, salt_diet_g / hfss_food_weight, 0)] #g salt per gram of HFSS food
    sp$pop[year>= sc_year, salt_combined:= (((hfss_food_weight*salt_density)*(1+tax_eff_salt)*(1+(extra_coverage_MLT*lbl_eff_salt)))-(hfss_food_weight*salt_density))*(1-comps) + 0.015 * (1 - 0.18)^((year - 19L) / 10L)] #negative value = reduction # the same as above #the effect was reduced for compensation (energy and salt) #Updated 20261103: (1-comps) applied to the policy effect only, not to the grain/dairy compensation term #Updated 20261104: labelling factor no longer applied to the grain/dairy substitution term (substitution is due to the tax only; grains/dairy are not HFSS and not labelled)
    sp$pop[is.na(salt_combined) | is.infinite(salt_combined) | salt_combined > 0, salt_combined := 0] #negative
    sp$pop[, salt_density := NULL] #Updated 20261103

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change :=  (salt_combined*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_combined >0, sbp_change :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*0.60, sbp_change :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change <0, sbp_change := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change) & (sbp_change >= 0),
                                   sbp_curr_xps - sbp_change,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_combined", "salt_diet_g", "sbp_change") := NULL]

 ## Reduction in fruits and vegetables due to taxation
    #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
    #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   137.4731 549.0456 #this is alpha, and beta
    elasticity_fruit <- draw(6, rbeta(1, 137.4731, 549.0456))*(-1)  #Updated 20261103: draw()

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)))] #tax_eff_fruit is negative #Updated 20261103


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   84.65086 84.65087 #this is alpha, and beta
    elasticity_veg <- draw(7, rbeta(1, 84.65086, 84.65087))*(-1)  #Updated 20261103: draw()

    tax_eff_veg <- tax_rate*elasticity_veg #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)))] #tax_eff_veg is negative #Updated 20261103

    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_tax_lbl_sn_comps", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_lbl_sn_comps")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)





######## SENSITIVITY ANALYSES WITH DIFFERENT COVERAGE FOR NW ##################################################################################################################################################
#For nutrient warning labelling, the authors apply effects
#only to the 25% of HFSS products assumed not to feature MTL labelling,
#based on overall packaged product data.
#However, they note that HFSS products may actually have lower FOPNL coverage
#(46% vs 59% for non-HFSS).
#Given this uncertainty could substantially affect results,
#could the authors provide a sensitivity analysis varying this coverage assumption?

################# Scenario 13: Song et al. (Consumer response only) [cons_sc1_lbl_sn_cover] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.009, 0.038, 0.068)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   8.893301 217.128711 #this is alpha, and beta
    lbl_eff_energy <- draw(8, rbeta(1, 8.893301, 217.128711))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
    #prop_packaged <- 0.80
    prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels

    #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the MLT become mandatory
    #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
    current_coverage_MLT <- 0.46 #updated for a sensitivity analysis 20260317, 46% of the HFSS products featuring FOPNL or multiple traffic light labelling, and therefore https://doi.org/10.1016/j.cdnut.2024.104509
    extra_coverage_MLT <- 1- current_coverage_MLT

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)*prop_packaged] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS energy explicitly as food weight (g) x energy density (kcal per gram of HFSS food).
    #Song et al. effect is a % change in energy per 100 g (i.e. in energy density), so: change = weight x density x % effect (= energy_diet x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not kcal.
    sp$pop[, energy_density := fifelse(hfss_food_weight > 0, energy_diet / hfss_food_weight, 0)] #kcal per gram of HFSS food
    sp$pop[year>= sc_year, energy_lbl:= (hfss_food_weight*energy_density)*(extra_coverage_MLT*lbl_eff_energy)]
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative
    sp$pop[, energy_density := NULL] #Updated 20261103

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_lbl", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_lbl") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_lbl*(4.2/1000))/1.5),
                                   20.7 * ((energy_lbl*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_lbl", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.021, 0.082)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   1.454906 53.419114 #this is alpha, and beta
    lbl_eff_salt <- draw(9, rbeta(1, 1.454906, 53.419114))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)*prop_packaged] #see above #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS salt explicitly as food weight (g) x salt density (g salt per gram of HFSS food).
    #Song et al. effect is a % change in salt per 100 g (i.e. in salt density), so: change = weight x density x % effect (= salt_diet_g x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not grams of salt.
    sp$pop[, salt_density := fifelse(hfss_food_weight > 0, salt_diet_g / hfss_food_weight, 0)] #g salt per gram of HFSS food
    sp$pop[year>= sc_year, salt_lbl:= (hfss_food_weight*salt_density)*(extra_coverage_MLT*lbl_eff_salt)]
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative
    sp$pop[, salt_density := NULL] #Updated 20261103

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_lbl >0, sbp_change_lbl :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*prop_packaged*0.60, sbp_change_lbl :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_lbl <0, sbp_change_lbl := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_lbl) & (sbp_change_lbl >= 0),
                                   sbp_curr_xps - sbp_change_lbl,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_lbl", "salt_diet_g", "sbp_change_lbl") := NULL]
    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_lbl_sn_cover", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_lbl_sn_cover")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)




######## SENSITIVITY ANALYSES WITH DIFFERENT EFFECTS FOR NW ##################################################################################################################################################
#We conducted a sensitivity analysis based on real-world evidence from Chile. https://doi.org/10.1371/journal.pmed.1004463
#We used evidence from the Phase 1 (less strict nutrient thresholds for classifying high-in products compared to Phase 2) in Chile,
#and conservatively assumed that energy and salt in purchases of HFSS products to be reduced by 17% and 15%, respectively,
#applied to 25% HFSS products without FOPNL in the UK to avoid overestimation.
#These effects are greater than the effects from the meta-analysis we used to inform the main modelling analysis.

################# Scenario 14: Evidence from Chile Phase 1 (Consumer response only) [cons_sc1_lbl_sensi_eff] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    # Effect = 17.1% (15.3% - 18.9%) https://doi.org/10.1371/journal.pmed.1004463
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.153, 0.171, 0.189)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   288.4971 1397.5370 #this is alpha, and beta
    lbl_eff_energy <- draw(8, rbeta(1, 288.4971, 1397.5370))*(-1)  #Updated 20261103: draw()
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
    #prop_packaged <- 0.80
    prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels

    #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the MLT become mandatory
    #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
    current_coverage_MLT <- 0.75
    extra_coverage_MLT <- 1- current_coverage_MLT

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)*prop_packaged] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    sp$pop[year>= sc_year, energy_lbl:= energy_diet*(extra_coverage_MLT*lbl_eff_energy)]
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_lbl", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_lbl") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_lbl*(4.2/1000))/1.5),
                                   20.7 * ((energy_lbl*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_lbl", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    # Effect = 15.2% (12.9% - 17.6%) https://doi.org/10.1371/journal.pmed.1004463
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.129, 0.152, 0.176)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   1.454906 53.419114 #this is alpha, and beta
    lbl_eff_salt <- draw(9, rbeta(1, 136.3641, 759.3622))*(-1)  #Updated 20261103: draw()
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)*prop_packaged] #see above #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_lbl:= salt_diet_g*(extra_coverage_MLT*lbl_eff_salt)]
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_lbl >0, sbp_change_lbl :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*prop_packaged*0.60, sbp_change_lbl :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_lbl <0, sbp_change_lbl := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_lbl) & (sbp_change_lbl >= 0),
                                   sbp_curr_xps - sbp_change_lbl,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_lbl", "salt_diet_g", "sbp_change_lbl") := NULL]
    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_lbl_sensi_eff", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_lbl_sensi_eff")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


################# Sensitivity Scenario 15: Taxation 10% (Consumer response only) and Subsidies 15% on FNVs [cons_sc1_10tax_15sub] #Updated 20261103
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- draw(1, runif(1, min = 0.5, max = 0.9)) #per 1% increase in price #based on discussion with Chris 20250805 see scenario 2 in https://doi.org/10.1016/j.jclepro.2023.138564  #Updated 20261103: draw()
    tax_eff_energy <- 0.10*elasticity_HFSS*(-1) #sensitivity 10% tax
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    sp$pop[year>= sc_year, energy_tax:= energy_diet*tax_eff_energy + 7.3 * (1 + 0.10)^((year - 19L) / 10L)] #2% effect of compensation
    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_tax", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_tax") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_tax*(4.2/1000))/1.5),
                                   20.7 * ((energy_tax*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_tax", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
           ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)] #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_tax:= salt_diet_g*tax_eff_salt + 0.018 * (1 - 0.18)^((year - 19L) / 10L)] #2% effect of compensation
    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_tax >0, sbp_change_tax :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*0.60, sbp_change_tax :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_tax <0, sbp_change_tax := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_tax) & (sbp_change_tax >= 0),
                                   sbp_curr_xps - sbp_change_tax,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_tax", "salt_diet_g", "sbp_change_tax") := NULL]

 ## Increased in fruits and vegetables due to subsidies
    #A 15% subsidy on fruits and vegetables resulted in an increase in the quantities demanded by 11% and 7% respectively (15% uncertainties were added)
    # mean = 0.11, sd = ((0.11*1.15) - (0.11*0.85))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.0935, 0.11, 0.1265)   # lower CI, mean, upper CI #here the effect directly, not elasticity
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]  153.5565 1240.4041 #this is alpha, and beta
    sub_eff_fruit <- draw(16, rbeta(1, 153.5565, 1240.4041))*(1) #+ due to subsidies  #Updated 20261103: draw()

    # sub_eff_fruit <- 0.15*elasticity_fruit #tax_rate is a proportion (0.08 = 8%) #15% subsidies #Updated 20261103
    # if (sub_eff_fruit > 0) sub_eff_fruit <- 0 #negative #Updated 20261103

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+sub_eff_fruit)))] #sub_eff_fruit is positive (subsidy) #Updated 20261103


    # mean = 0.07, sd = ((0.07*1.15) - (0.07*0.85))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.0595, 0.07, 0.0805)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   158.5147 2102.5560 #this is alpha, and beta
    sub_eff_veg <- draw(17, rbeta(1, 158.5147, 2102.5560))*(1) #+ due to subsidies  #Updated 20261103: draw()

    # sub_eff_veg <- 0.15*elasticity_veg #tax_rate is a proportion (0.08 = 8%) #15% subsidies #Updated 20261103
    # if (sub_eff_veg > 0) sub_eff_veg <- 0 #negative #Updated 20261103

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+sub_eff_veg)))] #sub_eff_veg is positive (subsidy) #Updated 20261103

    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_10tax_15sub", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_10tax_15sub")) #Updated 20261104: timing

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)





######## SENSITIVITY ANALYSES WITH WANING POLICY EFFECT (FOR CONSUMER RESPONSE ONLY!!!) #################################################################### #Updated 20261103: waning effect
### Same as the 11% compensation scenarios above (Scenarios 10-12, incl. grain/dairy substitution and, for the tax, the fruit/veg changes),
### but the whole policy effect wanes by 5% per year (compounding): effect x (1 - wane_rate)^(year - sc_year),
### i.e. 100% in 2026, 95% in 2027, 90.25% in 2028 ... about 42% in 2043. The grain/dairy substitution and fruit/veg changes wane with it.
### The per-parameter seeds (draw(k, ...)) are the same as in Scenarios 10-12, so differences from those scenarios reflect the waning only.
wane_rate <- 0.05 #Updated 20261103: waning effect

################# Scenario 16: Taxation 8%, 11% compensation, policy effect waning 5%/year (Consumer response only) [cons_sc1_tax_comps_wane] #Updated 20261103: waning effect
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- draw(1, runif(1, min = 0.6, max = 1.0)) #per 1% increase in price #based on discussion with Chris 20250805  #Updated 20261103: draw()
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    sp$pop[year>= sc_year, energy_tax:= (energy_diet*tax_eff_energy*(1-comps) + 5.8 * (1 + 0.10)^((year - 19L) / 10L))*(1 - wane_rate)^(year - sc_year)] #the effect was reduced for compensation #Updated 20261103: (1-comps) applied to the policy effect only, not to the grain/dairy compensation term #Updated 20261103: waning effect (whole effect incl. substitution x (1 - wane_rate)^(year - sc_year))
    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_tax", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_tax") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_tax*(4.2/1000))/1.5),
                                   20.7 * ((energy_tax*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_tax", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
           ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)] #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_tax:= (salt_diet_g*tax_eff_salt*(1-comps) + 0.015 * (1 - 0.18)^((year - 19L) / 10L))*(1 - wane_rate)^(year - sc_year)] #the effect was reduced for compensation (energy and salt) #Updated 20261103: (1-comps) applied to the policy effect only, not to the grain/dairy compensation term #Updated 20261103: waning effect (whole effect incl. substitution x (1 - wane_rate)^(year - sc_year))
    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_tax >0, sbp_change_tax :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*0.60, sbp_change_tax :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_tax <0, sbp_change_tax := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_tax) & (sbp_change_tax >= 0),
                                   sbp_curr_xps - sbp_change_tax,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_tax", "salt_diet_g", "sbp_change_tax") := NULL]

 ## Reduction in fruits and vegetables due to taxation
    #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
    #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   137.4731 549.0456 #this is alpha, and beta
    elasticity_fruit <- draw(6, rbeta(1, 137.4731, 549.0456))*(-1)  #Updated 20261103: draw()

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit*(1 - wane_rate)^(year - sc_year))))] #tax_eff_fruit is negative #Updated 20261103 #Updated 20261103: waning effect


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   84.65086 84.65087 #this is alpha, and beta
    elasticity_veg <- draw(7, rbeta(1, 84.65086, 84.65087))*(-1)  #Updated 20261103: draw()

    tax_eff_veg <- tax_rate*elasticity_veg #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg*(1 - wane_rate)^(year - sc_year))))] #tax_eff_veg is negative #Updated 20261103 #Updated 20261103: waning effect

    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_tax_comps_wane", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_comps_wane")) #Updated 20261104: timing


################# Scenario 17: Song et al. NW labelling, 11% compensation, policy effect waning 5%/year (Consumer response only) [cons_sc1_lbl_sn_comps_wane] #Updated 20261103: waning effect
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.009, 0.038, 0.068)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   8.893301 217.128711 #this is alpha, and beta
    lbl_eff_energy <- draw(8, rbeta(1, 8.893301, 217.128711))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
    #prop_packaged <- 0.80
    prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels

    #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the MLT become mandatory
    #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
    current_coverage_MLT <- 0.75
    extra_coverage_MLT <- 1- current_coverage_MLT

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)*prop_packaged] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS energy explicitly as food weight (g) x energy density (kcal per gram of HFSS food).
    #Song et al. effect is a % change in energy per 100 g (i.e. in energy density), so: change = weight x density x % effect (= energy_diet x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not kcal.
    sp$pop[, energy_density := fifelse(hfss_food_weight > 0, energy_diet / hfss_food_weight, 0)] #kcal per gram of HFSS food
    sp$pop[year>= sc_year, energy_lbl:= ((hfss_food_weight*energy_density)*(extra_coverage_MLT*lbl_eff_energy)*(1-comps))*(1 - wane_rate)^(year - sc_year)] #the effect was reduced for compensation #Updated 20261103: waning effect (whole effect incl. substitution x (1 - wane_rate)^(year - sc_year))
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative
    sp$pop[, energy_density := NULL] #Updated 20261103

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_lbl", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_lbl") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_lbl*(4.2/1000))/1.5),
                                   20.7 * ((energy_lbl*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_lbl", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
                  )]

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.021, 0.082)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   1.454906 53.419114 #this is alpha, and beta
    lbl_eff_salt <- draw(9, rbeta(1, 1.454906, 53.419114))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)*prop_packaged] #see above #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS salt explicitly as food weight (g) x salt density (g salt per gram of HFSS food).
    #Song et al. effect is a % change in salt per 100 g (i.e. in salt density), so: change = weight x density x % effect (= salt_diet_g x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not grams of salt.
    sp$pop[, salt_density := fifelse(hfss_food_weight > 0, salt_diet_g / hfss_food_weight, 0)] #g salt per gram of HFSS food
    sp$pop[year>= sc_year, salt_lbl:= ((hfss_food_weight*salt_density)*(extra_coverage_MLT*lbl_eff_salt)*(1-comps))*(1 - wane_rate)^(year - sc_year)] #the effect was reduced for compensation (energy and salt) #Updated 20261103: waning effect (whole effect incl. substitution x (1 - wane_rate)^(year - sc_year))
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative
    sp$pop[, salt_density := NULL] #Updated 20261103

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_lbl >0, sbp_change_lbl :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*prop_packaged*0.60, sbp_change_lbl :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change_lbl <0, sbp_change_lbl := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change_lbl) & (sbp_change_lbl >= 0),
                                   sbp_curr_xps - sbp_change_lbl,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_lbl", "salt_diet_g", "sbp_change_lbl") := NULL]
    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_lbl_sn_comps_wane", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_lbl_sn_comps_wane")) #Updated 20261104: timing


################# Scenario 18: Taxation 8% and Song et al. NW labelling, 11% compensation, policy effect waning 5%/year (Consumer response only) [cons_sc1_tax_lbl_sn_comps_wane] #Updated 20261103: waning effect
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    #setDT(sp$pop)
    # Ensure that .Random.seed exists (it is created when a random number is generated)
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }

    # Save the current random seed state
    saved_seed <- .Random.seed

    # Change the seed (this will change the state of the RNG)
    set.seed(93684 + sp$mc_aggr)
    #Updated 20261103: per-parameter seeds, so each parameter gets the same value in every scenario for the same mc
    #(common random numbers; the k for each parameter is the same across all scenarios, see list at the top of the script)
    draw <- function(k, expr) { set.seed(93684 + sp$mc_aggr * 100 + k); expr }

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- draw(1, runif(1, min = 0.6, max = 1.0)) #per 1% increase in price #based on discussion with Chris 20250805  #Updated 20261103: draw()
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    # label Song et al.
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.009, 0.038, 0.068)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   8.893301 217.128711 #this is alpha, and beta
    lbl_eff_energy <- draw(8, rbeta(1, 8.893301, 217.128711))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    #We assume that 80% of the products purchased are packaged (vs. 20% fresh) https://www.food.gov.uk/research/consumer-responses-to-food-labelling-a-rapid-evidence-review
    #prop_packaged <- 0.80
    prop_packaged <- 1 #based on meeting with Eric on 20250602, we used 100%, as all non-OOH HFSS can be added labels

    #Proportion of products that currently feature a multiple traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the MLT become mandatory
    #so, 25% for the effect of NW labeling for findings from MA suggesting similar effect for MLT and NW
    current_coverage_MLT <- 0.75
    extra_coverage_MLT <- 1- current_coverage_MLT

    sp$pop[, energy_diet:= hfss_capped(hfss_energy, total_energy)] #Updated 20261103: HFSS energy capped at total energy (see cap_hfss_to_total)
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    # The combined effect can be calculated as
    # 1 – [(1+PRi)*(1+PRj)] PR= percent reduction < 0
    #sp$pop[year>= sc_year, energy_combined:= (energy_diet*(1+tax_eff_energy)*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy))) - energy_diet]
    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS energy explicitly as food weight (g) x energy density (kcal per gram of HFSS food).
    #Song et al. effect is a % change in energy per 100 g (i.e. in energy density), so: change = weight x density x % effect (= energy_diet x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not kcal.
    sp$pop[, energy_density := fifelse(hfss_food_weight > 0, energy_diet / hfss_food_weight, 0)] #kcal per gram of HFSS food
    sp$pop[year>= sc_year, energy_combined:= ((((hfss_food_weight*energy_density)*(1+tax_eff_energy)*(1+(extra_coverage_MLT*lbl_eff_energy)))-(hfss_food_weight*energy_density))*(1-comps) + 5.8 * (1 + 0.10)^((year - 19L) / 10L))*(1 - wane_rate)^(year - sc_year)] #negative value = reduction # the same as above #the effect was reduced for compensation #Updated 20261103: (1-comps) applied to the policy effect only, not to the grain/dairy compensation term #Updated 20261103: waning effect (whole effect incl. substitution x (1 - wane_rate)^(year - sc_year)) #Updated 20261104: labelling factor no longer applied to the grain/dairy substitution term (substitution is due to the tax only; grains/dairy are not HFSS and not labelled)
    sp$pop[is.na(energy_combined) | is.infinite(energy_combined) | energy_combined > 0, energy_combined := 0] #negative
    sp$pop[, energy_density := NULL] #Updated 20261103

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Updated 20261103 (hallfast): START replaced Hall lookup table (Data_Hall_bmi_13_94_kcal_1_200.fst) with hallfast
      #Translating change in energy intake (kcal/day) into change in body weight using the
      #dynamic adult model of Hall et al. (2011), solved with hallfast (https://github.com/ChristK/hallfast).
      #Each person's body weight trajectory is simulated year by year from sc_year onwards
      #(cumulative effect; energy change can differ by year), see hall_update_bmi() at the top.
      hall_update_bmi(sp$pop, energy_col = "energy_combined", sc_year = sc_year)
      sp$pop[, c("energy_diet", "energy_combined") := NULL]
      #Updated 20261103 (hallfast): END

    } else{
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Christiansen model
      sp$pop[, change_bw:= fifelse(sex=="men", 17.7 * ((energy_combined*(4.2/1000))/1.5),
                                   20.7 * ((energy_combined*(4.2/1000))/1.5))]

      ## As the formula was for an absolute change in body weight, we calculated % change in body weight using baseline body weight from BMI,
      sp$pop[sex=="men", weight_assum:= bmi_curr_xps*(1.76*1.76)] #175.9 https://digital.nhs.uk/data-and-information/publications/statistical/health-survey-for-england/2021/part-4-trends#:~:text=Part%204:%20Trends-,Height%20and%20weight,For%20more%20information:%20Table%207.
      sp$pop[sex=="women", weight_assum:= bmi_curr_xps*(1.62*1.62)] #162.4

      sp$pop[, relative_change_bw:= (weight_assum + change_bw)/weight_assum] #change_bw is negative
      sp$pop[is.na(relative_change_bw)| is.nan(relative_change_bw) | is.infinite(relative_change_bw) | relative_change_bw < 0, relative_change_bw := 1] #if there is any unreliable value, we assigned 1 (indicating no change)

      sp$pop[year >= sc_year,
             bmi_curr_xps := fifelse((bmi_curr_xps > relative_change_bw) & (relative_change_bw > 0), #as percentage, max relative_change_bw=1
                                     bmi_curr_xps* relative_change_bw,
                                     bmi_curr_xps)]
      #sp$pop[, summary(bmi_curr_xps)]
      sp$pop[, c("energy_diet", "energy_combined", "change_bw", "weight_assum", "relative_change_bw") := NULL]
    }

    sp$pop[year>= sc_year, bmi_chnge := bmi_curr_xps - bmi_curr_xps_old] #negative
    sp$pop[year>= sc_year & bmi_chnge>0, bmi_chnge:=0]

    #adding mediated effect through SBP and Tchol #20250901
    #https://doi.org/10.1111/jch.14661
    #The SBP reduced by 5.79 mmHg (95% CI, 3.54–8.05) after a mean BMI reduction of 2.27 kg/m2
    #Only applying effect to those with BMI >= 25

    #https://doi.org/10.1210/clinem/dgaa673
    #Using the pharmacotherapy effect for mean reduction in tchol from a 1 kg/m2 in BMI after 12 months: -8.87 (–11.79,–5.94)mg/DL (supp table 3)
    #which is -0.229mmol/L (-0.305, -0.154)

    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- draw(3, rnorm(1, mean = -5.79, sd= sbp_sd))/2.27 #change provided is for 2.27 decrease in BMI  #Updated 20261103: draw()
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- draw(4, rnorm(1, mean = -0.229, tchol_sd))  #Updated 20261103: draw()
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is not reversed #Updated 20261103

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #baseline BMI (before the scenario) >= 25 #Updated 20261103
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity %in% c("black african", "black caribbean")) #Updated 20261103: include black caribbean (Mozaffarian coefficient is for Black people)
    )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0, 0.021, 0.082)   # lower CI, mean, upper CI #the opposite effect was converted to 0
    # p <- c(0, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   1.454906 53.419114 #this is alpha, and beta
    lbl_eff_salt <- draw(9, rbeta(1, 1.454906, 53.419114))*(-1) #per 100 grams  #Updated 20261103: draw()
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_capped(hfss_salt, total_salt_fooddiary)] #Updated 20261103: HFSS salt capped at total salt (see cap_hfss_to_total)
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20261103: removed hfss_food_weight/100 multiplier and wrote HFSS salt explicitly as food weight (g) x salt density (g salt per gram of HFSS food).
    #Song et al. effect is a % change in salt per 100 g (i.e. in salt density), so: change = weight x density x % effect (= salt_diet_g x % effect).
    #Note: hfss_food_weight x % effect alone would be grams of food, not grams of salt.
    sp$pop[, salt_density := fifelse(hfss_food_weight > 0, salt_diet_g / hfss_food_weight, 0)] #g salt per gram of HFSS food
    sp$pop[year>= sc_year, salt_combined:= ((((hfss_food_weight*salt_density)*(1+tax_eff_salt)*(1+(extra_coverage_MLT*lbl_eff_salt)))-(hfss_food_weight*salt_density))*(1-comps) + 0.015 * (1 - 0.18)^((year - 19L) / 10L))*(1 - wane_rate)^(year - sc_year)] #negative value = reduction # the same as above #the effect was reduced for compensation (energy and salt) #Updated 20261103: (1-comps) applied to the policy effect only, not to the grain/dairy compensation term #Updated 20261103: waning effect (whole effect incl. substitution x (1 - wane_rate)^(year - sc_year)) #Updated 20261104: labelling factor no longer applied to the grain/dairy substitution term (substitution is due to the tax only; grains/dairy are not HFSS and not labelled)
    sp$pop[is.na(salt_combined) | is.infinite(salt_combined) | salt_combined > 0, salt_combined := 0] #negative
    sp$pop[, salt_density := NULL] #Updated 20261103

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- draw(5, runif(1, min = 0, max = 1))  #Updated 20261103: draw()
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # effect per 100 mmol/d sodium (= 5.85 g/d salt); divided by 5.85 below to get per g/d of salt #Updated 20261103
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change :=  (salt_combined*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value = reduction in SBP #Updated 20261103
    sp$pop[salt_combined >0, sbp_change :=  0] #no change in SBP if no reduction in salt
    #Based on Mozaffarian et al. (see above), mean sodium of 1500 mg or 3.75 g of salt (1500/1000*2.5) is associated with reduced SBP in RCTs.
    sp$pop[salt_diet_g <=3.75*0.60, sbp_change :=  0] #no effect in SBP change if baseline salt <= 3.75 grams/day #60% of salt from HFSS
    sp$pop[sbp_change <0, sbp_change := 0]

    sp$pop[year >= sc_year,
           sbp_curr_xps := fifelse((sbp_curr_xps > sbp_change) & (sbp_change >= 0),
                                   sbp_curr_xps - sbp_change,
                                   sbp_curr_xps)]
    #sp$pop[, summary(sbp_curr_xps)]
    sp$pop[, c("hypertensive", "black", "salt_combined", "salt_diet_g", "sbp_change") := NULL]

 ## Reduction in fruits and vegetables due to taxation
    #reductions in fruits and vegetables by 2% and 5%, respectively due to a 10% tax (15% uncertainties were added)
    #elasticity_fruit <- rnorm(1, mean = -0.2, sd = ((-0.2*0.85) - (-0.2*1.15))/3.92) #((-0.2-(-0.2*0.15))-(-0.2+(-0.2*0.15)))/3.92
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.17, 0.2, 0.23)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   137.4731 549.0456 #this is alpha, and beta
    elasticity_fruit <- draw(6, rbeta(1, 137.4731, 549.0456))*(-1)  #Updated 20261103: draw()

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit*(1 - wane_rate)^(year - sc_year))))] #tax_eff_fruit is negative #Updated 20261103 #Updated 20261103: waning effect


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Beta distribution for percentage, as suggested by Chris: #Updated 20261103
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.425, 0.5, 0.575)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975)     # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   84.65086 84.65087 #this is alpha, and beta
    elasticity_veg <- draw(7, rbeta(1, 84.65086, 84.65087))*(-1)  #Updated 20261103: draw()

    tax_eff_veg <- tax_rate*elasticity_veg #tax_rate is a proportion (0.08 = 8%) #Updated 20261103
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg*(1 - wane_rate)^(year - sc_year))))] #tax_eff_veg is negative #Updated 20261103 #Updated 20261103: waning effect

    # Restore the original seed state
    assign(".Random.seed", saved_seed, envir = .GlobalEnv) #Updated 20261103: restore the global seed (".Random.seed <- saved_seed" only created a local copy)
  }
)

time_step("cons_sc1_tax_lbl_sn_comps_wane", IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_lbl_sn_comps_wane")) #Updated 20261104: timing


print("SUMMARY EXPORTS NOW!!!")

time_step("export_summaries", IMPACTncd$export_summaries(multicore = TRUE)) #Updated 20261104: timing

print("CREATE TABLES NOW!!!")

time_step("export_tables", IMPACTncd$export_tables( #Updated 20261104: timing
  baseline_year_for_change_outputs = 2026L,
  two_agegrps = FALSE,
  wtp = c(25000, 35000, 70000),
  multicore = TRUE,
  strata = list(
    ons = list(
      "year",
      c("year", "sex"),
      c("year", "qimd"),
      c("year", "sex", "qimd"),
      c("year", "agegrp"),
      c("year", "agegrp", "sex"),
      c("year", "agegrp", "sex", "qimd")
    ),
    esp = list(
      "year",
      c("year", "sex"),
      c("year", "qimd"),
      c("year", "sex", "qimd")
    ),
    mrtl_ons = list(
      "year",
      c("year", "sex"),
      c("year", "qimd"),
      c("year", "agegrp"),
      c("year", "agegrp", "sex"),
      c("year", "agegrp", "sex", "qimd")
    ),
    mrtl_esp = list(
      "year",
      c("year", "sex"),
      c("year", "qimd"),
      c("year", "sex", "qimd")
    ),
    disease_char = list(
      "year",
      c("year", "sex"),
      c("year", "qimd"),
      c("year", "sex", "qimd")
    ),
    xps_ons = list(
      "year",
      c("year", "sex"),
      c("year", "agegrp20"),
      c("year", "agegrp20", "sex"),
      c("year", "agegrp20", "sex", "qimd")
    ),
    xps_esp = list(
      "year",
      c("year", "sex"),
      c("year", "qimd"),
      c("year", "sex", "qimd")
    )
  )
)) #Updated 20261104: timing


#Updated 20261104: timing summary
print(run_times)
message(sprintf("TOTAL time: %.1f min (%.2f h)", sum(run_times$minutes), sum(run_times$minutes) / 60))
