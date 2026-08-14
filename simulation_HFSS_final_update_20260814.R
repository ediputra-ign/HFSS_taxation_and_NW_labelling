#rm(list=ls())
#gc()
#setwd("/mnt/storage_fast4/Salt_England/IMPACTncd_Engl")
source("./global.R") #Updated_by_Edi

#remotes::install_version("data.table", version = "1.16.4", repos = "http://cran.r-project.org")

# source("./additional_costs.R")
# source("./additional_prev_change.R")

IMPACTncd <- Simulation$new("./sim_design_Edi_final.yaml")
#after running the line above, I got Warning message:
#file(s) not found in repo "ChristK/IMPACTncd_Engl"

start_runs <- 1
n_runs <- 100L
tax_rate <- 0.08 #change for sensitivity analysis here

Hall <- TRUE #if FALSE, will be using Christiansen & Garby energy model

########## BASELINE SCENARIO
IMPACTncd$
  #del_logs()$
  #del_outputs()$
  run(mc = start_runs:n_runs, multicore = TRUE, scenario_nam = "sc0")#$
  #export_summaries(multicore = TRUE)





########## PATHWAYS THROUGH ENERGY AND SALT

########## Translating changes in energy kcal into body weight/BMI
########## Hall Method was used. https://github.com/INSP-RH/bw

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

################# Scenario 1: Taxation 8% (Consumer response only)
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

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- runif(1, min = 0.6, max = 1.0) #per 1% increase in price #based on discussion with Chris 20250805
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    sp$pop[, energy_diet:= hfss_energy]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~5.8 kcal for energy compensation in 2019 (10% increase in 10 years)
    sp$pop[year>= sc_year, energy_tax:= energy_diet*tax_eff_energy + 5.8 * (1 + 0.10)^((year - 19L) / 10L)]
    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                 as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_tax))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_tax", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
           ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
                  )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    sp$pop[, salt_diet_g:= hfss_salt]
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~0.015 g for salt compensation in 2019 (18% reduction in 10 years)
    sp$pop[year>= sc_year, salt_tax:= salt_diet_g*tax_eff_salt + 0.015 * (1 - 0.18)^((year - 19L) / 10L)]
    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_fruit <- rbeta(1, 137.4731, 549.0456)*(-1)

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax rate in %
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_veg <- rbeta(1, 84.65086, 84.65087)*(-1)

    tax_eff_veg <- tax_rate*elasticity_veg #tax rate in %
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative

    # Restore the original seed state
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax")

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


################# Scenario 1: Taxation 8% (Consumer response only, differences in the effects by INCOME)
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
    r_new <- runif(1, min = 0, max = 1) # to make it consistent

    tax_eff_energy1 <- (qbeta(r_new, 11.46264, 98.38558))*(-1)
    if (tax_eff_energy1 > 0) tax_eff_energy1 <- 0 #negative

    tax_eff_energy3 <- (qbeta(r_new, 12.46944, 197.57645))*(-1)
    if (tax_eff_energy3 > 0) tax_eff_energy3 <- 0 #negative

    tax_eff_energy5 <- (qbeta(r_new, 2.307398, 84.614135))*(-1)
    if (tax_eff_energy5 > 0) tax_eff_energy5 <- 0 #negative

    tax_eff_energy2 <- (tax_eff_energy1 + tax_eff_energy3)/2
    tax_eff_energy4 <- (tax_eff_energy3 + tax_eff_energy5)/2

    sp$pop[, energy_diet:= hfss_energy]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    # 1 Highest         2         3         4  5 Lowest
    # 379042    349334    362763    369401    352700
    # sp$pop[year>= sc_year, energy_tax:= energy_diet*tax_eff_energy]
    sp$pop[year>= sc_year & income == "1 Lowest", energy_tax:= energy_diet*tax_eff_energy1] #the effect was reversed
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
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                                  as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_tax))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_tax", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
    )]

    tax_eff_salt1 <- tax_eff_energy1 # the effect was assumed the same (the proportional change)
    tax_eff_salt2 <- tax_eff_energy2
    tax_eff_salt3 <- tax_eff_energy3
    tax_eff_salt4 <- tax_eff_energy4
    tax_eff_salt5 <- tax_eff_energy5

    sp$pop[, salt_diet_g:= hfss_salt]
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #sp$pop[year>= sc_year, salt_tax:= salt_diet_g*tax_eff_salt]
    sp$pop[year>= sc_year & income == "1 Lowest", salt_tax:= salt_diet_g*tax_eff_salt1] #the effect was reversed
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
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_fruit <- rbeta(1, 137.4731, 549.0456)*(-1)

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax rate in %
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_veg <- rbeta(1, 84.65086, 84.65087)*(-1)

    tax_eff_veg <- tax_rate*elasticity_veg #tax rate in %
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative

    # Restore the original seed state
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_effINC")

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


################# Scenario 1: Taxation 8% (Consumer response only, without reduction in fruits and veggies)
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

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- runif(1, min = 0.6, max = 1.0) #per 1% increase in price #based on discussion with Chris 20250805
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    sp$pop[, energy_diet:= hfss_energy]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~5.8 kcal for energy compensation in 2019 (10% increase in 10 years)
    sp$pop[year>= sc_year, energy_tax:= energy_diet*tax_eff_energy + 5.8 * (1 + 0.10)^((year - 19L) / 10L)]
    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                 as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_tax))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_tax", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
                  )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    sp$pop[, salt_diet_g:= hfss_salt]
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~0.015 g for salt compensation in 2019 (18% reduction in 10 years)
    sp$pop[year>= sc_year, salt_tax:= salt_diet_g*tax_eff_salt + 0.015 * (1 - 0.18)^((year - 19L) / 10L)]
    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
 #    # using Betta distribution for percentage, as suggested by Chris:
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
 #    tax_eff_fruit <- tax_rate*elasticity_fruit #tax rate in %
 #    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative
 #
 #    sp$pop[year >= sc_year,
 #           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric
 #
 #
 #    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
 #    # using Betta distribution for percentage, as suggested by Chris:
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
 #    tax_eff_veg <- tax_rate*elasticity_veg #tax rate in %
 #    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative
 #
 #    sp$pop[year >= sc_year,
 #           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative

    # Restore the original seed state
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_nfvg")
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

################# Scenario 1: Song et al. (Consumer response only)
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

 ## Reduction in energy
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_energy <- rbeta(1, 8.893301, 217.128711)*(-1) #per 100 grams
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

    sp$pop[, energy_diet:= hfss_energy*prop_packaged]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    sp$pop[year>= sc_year, energy_lbl:= energy_diet*(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy)]
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                 as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_lbl))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_lbl", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
                  )]

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_salt <- rbeta(1, 1.454906, 53.419114)*(-1) #per 100 grams
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_salt*prop_packaged] #see above
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_lbl:= salt_diet_g*(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_salt)]
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_lbl_sn")

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


################# Scenario 1: Song et al. (Consumer response only, sensitivity reduction in total intakes)
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

 ## Reduction in energy
    #lbl_eff_energy1 <- rnorm(1, mean = -0.129, sd = (-0.08 - (-0.179))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_energy1 <- rbeta(1, 24.06169, 160.75193)*(-1)
    if (lbl_eff_energy1 > 0) lbl_eff_energy1 <- 0 #negative

    #lbl_eff_energy2 <- rnorm(1, mean = -0.064, sd = (-0.004 - (-0.125))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
    # install.packages("rriskDistributions")
    # library(rriskDistributions)
    # q <- c(0.004, 0.064, 0.125)   # lower CI, mean, upper CI
    # p <- c(0.025, 0.5, 0.975) # percentiles: 2.5%, 50%, 97.5%
    #
    # params <- get.beta.par(q = q, p = p)
    # params
    #
    # $par
    # [1]   24.06169 160.75193 #this is alpha, and beta
    lbl_eff_energy2 <- rbeta(1, 6.360689, 88.532728)*(-1)
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

    sp$pop[, energy_diet:= hfss_energy*prop_packaged]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    #Based on different effect sizes 12.9% vs 6.4%, we will apply 6.4% to 75% of the products currently featuring MLT as the effect of NW outperform NL
    sp$pop[year>= sc_year, energy_lbl:= energy_diet*((current_coverage_MLT*lbl_eff_energy2) + (extra_coverage_MLT*lbl_eff_energy1))]
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                 as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_lbl))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_lbl", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
                  )]

    #lbl_eff_salt1 <- rnorm(1, mean = -0.078, sd = (0.168 - (-0.012))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_salt1 <- rbeta(1, 4.584322, 50.638770)*(-1)
    if (lbl_eff_salt1 > 0) lbl_eff_salt1 <- 0 #negative

    #lbl_eff_salt2 <- rnorm(1, mean = -0.013, sd = (0.103 - (-0.077))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_salt2 <- rbeta(1, 0.6301193, 26.3996677)*(-1)
    if (lbl_eff_salt2 > 0) lbl_eff_salt2 <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_salt*prop_packaged] #see above
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_lbl:= salt_diet_g*((current_coverage_MLT*lbl_eff_salt2) + (extra_coverage_MLT*lbl_eff_salt1))]
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_lbl_sn_sensi")

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
#       #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
#       #https://github.com/INSP-RH/bw
#
#       change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
#                                   as.data.table = TRUE)
#
#       sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
#                     calories=as.integer(energy_lbl))]
#
#       sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
#       sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
#       sp$pop[calories < (-200), calories:= -200] #only few individuals
#
#       #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps
#
#       years <- sc_year: max(sp$pop$year)
#       for (yr in years) {
#         syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
#         # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
#         syn_year <- merge(syn_year, change_bmi_Hall,
#                           by = c("age", "sex", "calories", "bmi"),
#                           all.x = TRUE)
#
#         setorder(syn_year, pid, year)
#         syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
#         syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
#         syn_year[, change_bmi_last := NULL]
#
#         # Update the original syn table with new BMI for this year
#         sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
#       }
#       sp$pop[, c("energy_diet", "energy_lbl", "calories", "bmi") := NULL]
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
################# Scenario 1: Taxation 8% and Song et al. (Consumer response only)
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

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- runif(1, min = 0.6, max = 1.0) #per 1% increase in price #based on discussion with Chris 20250805
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    # label Song et al.
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_energy <- rbeta(1, 8.893301, 217.128711)*(-1) #per 100 grams
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

    sp$pop[, energy_diet:= hfss_energy]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~5.8 kcal for energy compensation in 2019 (10% increase in 10 years)
    # The combined effect can be calculated as
    # 1 – [(1+PRi)*(1+PRj)] PR= percent reduction < 0
    sp$pop[year>= sc_year, energy_combined:= ((energy_diet*(1+tax_eff_energy)+ 5.8 * (1 + 0.10)^((year - 19L) / 10L))*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy))) - energy_diet]
    #sp$pop[year>= sc_year, energy_combined:= energy_diet*(((1+tax_eff_energy)*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy)))-1)] # -1 to make it negative # the same as above
    sp$pop[is.na(energy_combined) | is.infinite(energy_combined) | energy_combined > 0, energy_combined := 0] #negative

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
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                                  as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_combined))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_combined", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
    )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_salt <- rbeta(1, 1.454906, 53.419114)*(-1) #per 100 grams
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_salt]
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~0.015 g for salt compensation in 2019 (18% reduction in 10 years)
    sp$pop[year>= sc_year, salt_combined:= ((salt_diet_g*(1+tax_eff_salt) + 0.015 * (1 - 0.18)^((year - 19L) / 10L))*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_salt))) - salt_diet_g]
    #sp$pop[year>= sc_year, salt_combined:= salt_diet_g*(((1+tax_eff_salt)*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_salt)))-1)] # -1 to make it negative
    sp$pop[is.na(salt_combined) | is.infinite(salt_combined) | salt_combined > 0, salt_combined := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change :=  (salt_combined*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_fruit <- rbeta(1, 137.4731, 549.0456)*(-1)

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax rate in %
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_veg <- rbeta(1, 84.65086, 84.65087)*(-1)

    tax_eff_veg <- tax_rate*elasticity_veg #tax rate in %
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative

    # Restore the original seed state
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_lbl_sn")

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
#       #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
#       #https://github.com/INSP-RH/bw
#
#       change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
#                                   as.data.table = TRUE)
#
#       sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
#                     calories=as.integer(energy_combined))]
#
#       sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
#       sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
#       sp$pop[calories < (-200), calories:= -200] #only few individuals
#
#       #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps
#
#       years <- sc_year: max(sp$pop$year)
#       for (yr in years) {
#         syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
#         # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
#         syn_year <- merge(syn_year, change_bmi_Hall,
#                           by = c("age", "sex", "calories", "bmi"),
#                           all.x = TRUE)
#
#         setorder(syn_year, pid, year)
#         syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
#         syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
#         syn_year[, change_bmi_last := NULL]
#
#         # Update the original syn table with new BMI for this year
#         sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
#       }
#       sp$pop[, c("energy_diet", "energy_combined", "calories", "bmi") := NULL]
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
#       #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
#       #https://github.com/INSP-RH/bw
#
#       change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
#                                   as.data.table = TRUE)
#
#       sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
#                     calories=as.integer(energy_combined))]
#
#       sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
#       sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
#       sp$pop[calories < (-200), calories:= -200] #only few individuals
#
#       #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps
#
#       years <- sc_year: max(sp$pop$year)
#       for (yr in years) {
#         syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
#         # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
#         syn_year <- merge(syn_year, change_bmi_Hall,
#                           by = c("age", "sex", "calories", "bmi"),
#                           all.x = TRUE)
#
#         setorder(syn_year, pid, year)
#         syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
#         syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
#         syn_year[, change_bmi_last := NULL]
#
#         # Update the original syn table with new BMI for this year
#         sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
#       }
#       sp$pop[, c("energy_diet", "energy_tax", "energy_lbl", "energy_combined", "calories", "bmi") := NULL]
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
################# Scenario 1: Taxation 8% (Consumer response and reformulation)
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

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- runif(1, min = 0.6, max = 1.0) #per 1% increase in price #based on discussion with Chris 20250805
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    # Reformulation effect
    #ref_eff_energy <- rnorm(1, mean = -0.06, sd = ((-0.06*0.85) - (-0.06*1.15))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    ref_eff_energy <- rbeta(1, 162.5123, 2541.8669)*(-1)
    if (ref_eff_energy > 0) ref_eff_energy <- 0 #negative

    # We assume that 40% of the manufactures reformulate the products to avoid taxation
    prop_ref <- 0.40
    prop_tax <- 1 - prop_ref #the proportion of products with taxation

    sp$pop[, energy_diet:= hfss_energy]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~5.8 kcal for energy compensation in 2019 (10% increase in 10 years)
    sp$pop[year>= sc_year, energy_tax:= energy_diet*((prop_tax*tax_eff_energy) + (prop_ref*ref_eff_energy)) + (prop_tax* (5.8 * (1 + 0.10)^((year - 19L) / 10L)))]
    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                 as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_tax))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_tax", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
                  )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    # Reformulation effect
    #ref_eff_salt <- rnorm(1, mean = -0.15, sd = ((-0.15*0.85) - (-0.15*1.15))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    ref_eff_salt <- rbeta(1, 146.4041, 828.3033)*(-1)
    if (ref_eff_salt > 0) ref_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_salt]
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~0.015 g for salt compensation in 2019 (18% reduction in 10 years)
    sp$pop[year>= sc_year, salt_tax:= salt_diet_g*((prop_tax*tax_eff_salt) + (prop_ref*ref_eff_salt)) + prop_tax*(0.015 * (1 - 0.18)^((year - 19L) / 10L))]
    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_fruit <- rbeta(1, 137.4731, 549.0456)*(-1)

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax rate in %
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_veg <- rbeta(1, 84.65086, 84.65087)*(-1)

    tax_eff_veg <- tax_rate*elasticity_veg #tax rate in %
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative

    # Restore the original seed state
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_ref_sc1_tax")

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


######## NUTRITION WARNING LABELLING ##################################################################################################################################################
################# Scenario 1: Song et al. (Consumer response and reformulation)
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

 ## Reduction in energy
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_energy <- rbeta(1, 8.893301, 217.128711)*(-1) #per 100 grams
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    # Reformulation effect 0.9% (95% CI: [–3.1, 4.9])
    #ref_eff_energy <- rnorm(1, mean = -0.009, sd = ((0.049- (-0.031))/3.92))
    # using Betta distribution for percentage, as suggested by Chris:
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
    # using Betta distribution for percentage, as suggested by Chris:
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
    ref_eff_energy <- rbeta(1, 1.932238, 40.036484)*(-1)
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

    sp$pop[, energy_diet:= hfss_energy*prop_packaged]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    sp$pop[year>= sc_year, energy_lbl:= energy_diet*(((extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy)*prop_lbl) + (ref_eff_energy*prop_ref))]
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                 as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_lbl))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_lbl", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
                  )]

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_salt <- rbeta(1, 1.454906, 53.419114)*(-1) #per 100 grams
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative


    # Reformulation effect 8.9% (95% CI: [–17.3%, –0.6%]
    #ref_eff_salt <- rnorm(1, mean = -0.089, sd = ((-0.006- (-0.173))/3.92))
    # using Betta distribution for percentage, as suggested by Chris:
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
    # using Betta distribution for percentage, as suggested by Chris:
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
    ref_eff_salt <- rbeta(1, 3.940064, 35.985975)*(-1)
    if (ref_eff_salt > 0) ref_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_salt*prop_packaged] #see above
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_lbl:= salt_diet_g*(((extra_coverage_MLT*hfss_food_weight/100*lbl_eff_salt)*prop_lbl) + (ref_eff_salt*prop_ref))]
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_ref_sc1_lbl_sn")

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
#       #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
#       #https://github.com/INSP-RH/bw
#
#       change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
#                  as.data.table = TRUE)
#
#       sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
#                     calories=as.integer(energy_lbl))]
#
#       sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
#       sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
#       sp$pop[calories < (-200), calories:= -200] #only few individuals
#
#       #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps
#
#       years <- sc_year: max(sp$pop$year)
#       for (yr in years) {
#         syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
#         # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
#         syn_year <- merge(syn_year, change_bmi_Hall,
#                           by = c("age", "sex", "calories", "bmi"),
#                           all.x = TRUE)
#
#         setorder(syn_year, pid, year)
#         syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
#         syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
#         syn_year[, change_bmi_last := NULL]
#
#         # Update the original syn table with new BMI for this year
#         sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
#       }
#       sp$pop[, c("energy_diet", "energy_lbl", "calories", "bmi") := NULL]
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
################# Scenario 1: Taxation 8% and Song et al. (Consumer response and reformulation)
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

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- runif(1, min = 0.6, max = 1.0) #per 1% increase in price #based on discussion with Chris 20250805
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    # label Song et al.
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_energy <- rbeta(1, 8.893301, 217.128711)*(-1) #per 100 grams
    if (lbl_eff_energy > 0) lbl_eff_energy <- 0 #negative

    # Reformulation effect (We used from the taxation only as the coverage is higher)
    #ref_eff_energy <- rnorm(1, mean = -0.06, sd = ((-0.06*0.85) - (-0.06*1.15))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    ref_eff_energy <- rbeta(1, 162.5123, 2541.8669)*(-1)
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

    sp$pop[, energy_diet:= hfss_energy]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~5.8 kcal for energy compensation in 2019 (10% increase in 10 years)
    # The combined effect can be calculated as
    # 1 – [(1+PRi)*(1+PRj)] PR= percent reduction < 0
    sp$pop[year>= sc_year, energy_combined:= prop_tax*(((energy_diet*(1+tax_eff_energy) + 5.8 * (1 + 0.10)^((year - 19L) / 10L))*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy)))-energy_diet)
                                             + (energy_diet*prop_ref*ref_eff_energy)] # -1 to make it negative
    sp$pop[is.na(energy_combined) | is.infinite(energy_combined) | energy_combined > 0, energy_combined := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                                  as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_combined))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_combined", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
    )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_salt <- rbeta(1, 1.454906, 53.419114)*(-1) #per 100 grams
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    # Reformulation effect
    #ref_eff_salt <- rnorm(1, mean = -0.15, sd = ((-0.15*0.85) - (-0.15*1.15))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    ref_eff_salt <- rbeta(1, 146.4041, 828.3033)*(-1)
    if (ref_eff_salt > 0) ref_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_salt]
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    #Updated 20260805 #adding compensation for grains and dairy 1.6% for 8% tax ~0.015 g for salt compensation in 2019 (18% reduction in 10 years)
    sp$pop[year>= sc_year, salt_combined:= prop_tax*(((salt_diet_g*(1+tax_eff_salt) + 0.015 * (1 - 0.18)^((year - 19L) / 10L))*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_salt)))-salt_diet_g)
                                           + (salt_diet_g*prop_ref*ref_eff_salt)] # -1 to make it negative

    sp$pop[is.na(salt_combined) | is.infinite(salt_combined) | salt_combined > 0, salt_combined := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change :=  (salt_combined*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_fruit <- rbeta(1, 137.4731, 549.0456)*(-1)

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax rate in %
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_veg <- rbeta(1, 84.65086, 84.65087)*(-1)

    tax_eff_veg <- tax_rate*elasticity_veg #tax rate in %
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative

    # Restore the original seed state
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_ref_sc1_tax_lbl_sn")

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
#       #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
#       #https://github.com/INSP-RH/bw
#
#       change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
#                                   as.data.table = TRUE)
#
#       sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
#                     calories=as.integer(energy_combined))]
#
#       sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
#       sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
#       sp$pop[calories < (-200), calories:= -200] #only few individuals
#
#       #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps
#
#       years <- sc_year: max(sp$pop$year)
#       for (yr in years) {
#         syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
#         # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
#         syn_year <- merge(syn_year, change_bmi_Hall,
#                           by = c("age", "sex", "calories", "bmi"),
#                           all.x = TRUE)
#
#         setorder(syn_year, pid, year)
#         syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
#         syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
#         syn_year[, change_bmi_last := NULL]
#
#         # Update the original syn table with new BMI for this year
#         sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
#       }
#       sp$pop[, c("energy_diet", "energy_combined", "calories", "bmi") := NULL]
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
######## 11% COMPENSATION ON ENERGY ONLY!!! https://doi.org/10.1186/s12966-022-01287-z

comps <- 0.11

######## TAXATION ##################################################################################################################################################
################# Scenario 1: Taxation 8% (Consumer response only)
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

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- runif(1, min = 0.6, max = 1.0) #per 1% increase in price #based on discussion with Chris 20250805
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    sp$pop[, energy_diet:= hfss_energy]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    sp$pop[year>= sc_year, energy_tax:= (energy_diet*tax_eff_energy + 5.8 * (1 + 0.10)^((year - 19L) / 10L))*(1-comps)] #the effect was reduced for compensation
    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                 as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_tax))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_tax", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
           ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
                  )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    sp$pop[, salt_diet_g:= hfss_salt]
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_tax:= (salt_diet_g*tax_eff_salt + 0.015 * (1 - 0.18)^((year - 19L) / 10L))*(1-comps)] #the effect was reduced for compensation
    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_fruit <- rbeta(1, 137.4731, 549.0456)*(-1)

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax rate in %
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_veg <- rbeta(1, 84.65086, 84.65087)*(-1)

    tax_eff_veg <- tax_rate*elasticity_veg #tax rate in %
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative

    # Restore the original seed state
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_comps")
#
# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


######## NUTRITION WARNING LABELLING ##################################################################################################################################################
################# Scenario 1: Song et al. (Consumer response only)
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

 ## Reduction in energy
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_energy <- rbeta(1, 8.893301, 217.128711)*(-1) #per 100 grams
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

    sp$pop[, energy_diet:= hfss_energy*prop_packaged]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    sp$pop[year>= sc_year, energy_lbl:= energy_diet*(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy)*(1-comps)] #the effect was reduced for compensation
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                 as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_lbl))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_lbl", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
                  )]

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_salt <- rbeta(1, 1.454906, 53.419114)*(-1) #per 100 grams
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_salt*prop_packaged] #see above
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_lbl:= salt_diet_g*(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_salt)*(1-comps)] #the effect was reduced for compensation
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_lbl_sn_comps")

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


######## COMBINED EFFECTS OF TAXATION AND NUTRITION WARNING LABELLING ##################################################################################################################################################
################# Scenario 1: Taxation 8% and Song et al. (Consumer response only)
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

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- runif(1, min = 0.6, max = 1.0) #per 1% increase in price #based on discussion with Chris 20250805
    tax_eff_energy <- tax_rate*elasticity_HFSS*(-1)
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    # label Song et al.
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_energy <- rbeta(1, 8.893301, 217.128711)*(-1) #per 100 grams
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

    sp$pop[, energy_diet:= hfss_energy]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L

    # The combined effect can be calculated as
    # 1 – [(1+PRi)*(1+PRj)] PR= percent reduction < 0
    #sp$pop[year>= sc_year, energy_combined:= (energy_diet*(1+tax_eff_energy)*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy))) - energy_diet]
    sp$pop[year>= sc_year, energy_combined:= (((energy_diet*(1+tax_eff_energy) + 5.8 * (1 + 0.10)^((year - 19L) / 10L))*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy)))-energy_diet)*(1-comps)] # -1 to make it negative # the same as above #the effect was reduced for compensation
    sp$pop[is.na(energy_combined) | is.infinite(energy_combined) | energy_combined > 0, energy_combined := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                                  as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_combined))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_combined", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
    )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_salt <- rbeta(1, 1.454906, 53.419114)*(-1) #per 100 grams
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_salt]
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_combined:= (((salt_diet_g*(1+tax_eff_salt) + 0.015 * (1 - 0.18)^((year - 19L) / 10L))*(1+(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_salt)))-salt_diet_g)*(1-comps)] # -1 to make it negative # the same as above #the effect was reduced for compensation
    sp$pop[is.na(salt_combined) | is.infinite(salt_combined) | salt_combined > 0, salt_combined := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change :=  (salt_combined*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_fruit <- rbeta(1, 137.4731, 549.0456)*(-1)

    tax_eff_fruit <- tax_rate*elasticity_fruit #tax rate in %
    if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric


    #elasticity_veg <- rnorm(1, mean = -0.5, sd = ((-0.5*0.85) - (-0.5*1.15))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    elasticity_veg <- rbeta(1, 84.65086, 84.65087)*(-1)

    tax_eff_veg <- tax_rate*elasticity_veg #tax rate in %
    if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative

    # Restore the original seed state
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_tax_lbl_sn_comps")

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

################# Scenario 1: Song et al. (Consumer response only)
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

 ## Reduction in energy
    #lbl_eff_energy <- rnorm(1, mean = -0.038, sd = (0.068 - 0.009)/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_energy <- rbeta(1, 8.893301, 217.128711)*(-1) #per 100 grams
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

    sp$pop[, energy_diet:= hfss_energy*prop_packaged]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    sp$pop[year>= sc_year, energy_lbl:= energy_diet*(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_energy)]
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                 as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_lbl))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_lbl", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
                  )]

    #lbl_eff_salt <- rnorm(1, mean = -0.021, sd = (0.082 - (-0.041))/3.92) #per 100 grams
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_salt <- rbeta(1, 1.454906, 53.419114)*(-1) #per 100 grams
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_salt*prop_packaged] #see above
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_lbl:= salt_diet_g*(extra_coverage_MLT*hfss_food_weight/100*lbl_eff_salt)]
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_lbl_sn_cover")

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)




######## SENSITIVITY ANALYSES WITH DIFFERENT EFFECTS FOR NW ##################################################################################################################################################
#We conducted a sensitivity analysis based on real-world evidence from Chile. https://doi.org/10.1371/journal.pmed.1004463
#We used evidence from the Phase 1 (less strict nutrient thresholds for classifying high-in products compared to Phase 2) in Chile,
#and conservatively assumed that energy and salt in purchases of HFSS products to be reduced by 17% and 15%, respectively,
#applied to 25% HFSS products without FOPNL in the UK to avoid overestimation.
#These effects are greater than the effects from the meta-analysis we used to inform the main modelling analysis.

################# Scenario 1: Evidence from Chile Phase 1 (Consumer response only)
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

 ## Reduction in energy
    # Effect = 17.1% (15.3% - 18.9%) https://doi.org/10.1371/journal.pmed.1004463
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_energy <- rbeta(1, 288.4971, 1397.5370)*(-1)
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

    sp$pop[, energy_diet:= hfss_energy*prop_packaged]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    sp$pop[year>= sc_year, energy_lbl:= energy_diet*(extra_coverage_MLT*lbl_eff_energy)]
    sp$pop[is.na(energy_lbl) | is.infinite(energy_lbl) | energy_lbl > 0, energy_lbl := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                 as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_lbl))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_lbl", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
                  )]

    # Effect = 15.2% (12.9% - 17.6%) https://doi.org/10.1371/journal.pmed.1004463
    # using Betta distribution for percentage, as suggested by Chris:
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
    lbl_eff_salt <- rbeta(1, 136.3641, 759.3622)*(-1)
    if (lbl_eff_salt > 0) lbl_eff_salt <- 0 #negative

    sp$pop[, salt_diet_g:= hfss_salt*prop_packaged] #see above
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_lbl:= salt_diet_g*(extra_coverage_MLT*lbl_eff_salt)]
    sp$pop[is.na(salt_lbl) | is.infinite(salt_lbl) | salt_lbl > 0, salt_lbl := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_lbl :=  (salt_lbl*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_lbl_sensi_eff")

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)


################# Sensitivity Scenario: Taxation 10% (Consumer response only) and Subsidies 15% on FNVs
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

 ## Reduction in energy
    # Based on a study in the Scottish context, a weighted price elasticity = −0.939 (SD = 0.139) for HFSS products
    #tax_rate <- 0.08
    #elasticity_HFSS <- rnorm(1, mean = -0.939, sd = 0.139) #per 1% increase in price
    elasticity_HFSS <- runif(1, min = 0.5, max = 0.9) #per 1% increase in price #based on discussion with Chris 20250805 see scenario 2 in https://doi.org/10.1016/j.jclepro.2023.138564
    tax_eff_energy <- 0.10*elasticity_HFSS*(-1) #sensitivity 10% tax
    if (tax_eff_energy > 0) tax_eff_energy <- 0 #negative

    sp$pop[, energy_diet:= hfss_energy]
    sp$pop[is.na(energy_diet) | is.infinite(energy_diet) | energy_diet < 0, energy_diet := 0] #positive

    #We will estimate the impact from 2026
    sc_year <- 26L
    sp$pop[year>= sc_year, energy_tax:= energy_diet*tax_eff_energy + 7.3 * (1 + 0.10)^((year - 19L) / 10L)] #2% effect of compensation
    sp$pop[is.na(energy_tax) | is.infinite(energy_tax) | energy_tax > 0, energy_tax := 0] #negative

    sp$pop[year>= sc_year, bmi_curr_xps_old := bmi_curr_xps] #it will be used for later calculation

    if (Hall) {
      #Translating reduction in energy (kcal) into reduction in body weight (kg) using Hall model
      #https://github.com/INSP-RH/bw

      change_bmi_Hall <- read_fst("./Data_Hall_bmi_13_94_kcal_1_200.fst",
                 as.data.table = TRUE)

      sp$pop[, `:=`(bmi=as.integer(bmi_curr_xps),
                    calories=as.integer(energy_tax))]

      sp$pop[bmi< 13, bmi:= 13] #minimum available bmi 13
      sp$pop[bmi> 94, bmi:= 94] #maximum available bmi 94
      sp$pop[calories < (-200), calories:= -200] #only few individuals

      #sp$pop$bmi_curr_xps_ <- sp$pop$bmi_curr_xps

      years <- sc_year: max(sp$pop$year)
      for (yr in years) {
        syn_year <- copy(sp$pop[year == yr, c("age", "sex", "calories", "bmi", "pid", "year", "bmi_curr_xps")])
        # syn_year <- absorb_dt(syn_year, change_bmi_Hall)
        syn_year <- merge(syn_year, change_bmi_Hall,
                          by = c("age", "sex", "calories", "bmi"),
                          all.x = TRUE)

        setorder(syn_year, pid, year)
        syn_year[calories == 0, change_bmi_last := 0] #no change in bmi if calories reduction == 0
        syn_year[, bmi_curr_xps := bmi_curr_xps + change_bmi_last] #change_bmi is negative
        syn_year[, change_bmi_last := NULL]

        # Update the original syn table with new BMI for this year
        sp$pop[year == yr, bmi_curr_xps := syn_year$bmi_curr_xps]
      }
      sp$pop[, c("energy_diet", "energy_tax", "calories", "bmi") := NULL]

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
      # Restore the original seed state
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
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed

    tchol_sd  <- (0.305 - 0.154)/(2*(qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -0.229, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed

    sp$pop[bmi_curr_xps_old >= 25 & year >= sc_year, #bmi baseline?
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
           ]

    sp$pop[, c("bmi_curr_xps_old", "bmi_chnge") := NULL]

 ## Reduction in salt
    sp$pop[, `:=`(hypertensive = as.integer(sbp_curr_xps > 140),
                  black = as.integer(ethnicity=="black african")
                  )]

    tax_eff_salt <- tax_eff_energy # the effect was assumed the same (the proportional change)

    sp$pop[, salt_diet_g:= hfss_salt]
    sp$pop[is.na(salt_diet_g) | is.infinite(salt_diet_g) | salt_diet_g < 0, salt_diet_g := 0] #positive

    sp$pop[year>= sc_year, salt_tax:= salt_diet_g*tax_eff_salt + 0.018 * (1 - 0.18)^((year - 19L) / 10L)] #2% effect of compensation
    sp$pop[is.na(salt_tax) | is.infinite(salt_tax) | salt_tax > 0, salt_tax := 0] #negative

    #Translating reduction in salt (gram) into reduction in SBP (mmHG)
    #coefficient & SE for 100 mmol sodium or 5.85 g/d salt = α = -3.735 (0.730); B1 = -0.105 (0.029); B2 = -1.874 (0.884); and B3 = -2.489 (1.188 )
    r_new_sbp <- runif(1, min = 0, max = 1)
    a <- qnorm(r_new_sbp, mean = -3.735, sd = 0.730) # coefficient per g/d
    b1 <- qnorm(r_new_sbp, mean = -0.105, sd = 0.029)
    b2 <- qnorm(r_new_sbp, mean = -1.874, sd = 0.884)
    b3 <- qnorm(r_new_sbp, mean = -2.489, sd = 1.188)

    sp$pop[, sbp_change_tax :=  (salt_tax*(a + ((age-50)*b1) + (b2*hypertensive) +  (b3*black)))/5.85] #positive value
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
    # using Betta distribution for percentage, as suggested by Chris:
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
    tax_eff_fruit <- rbeta(1, 153.5565, 1240.4041)*(1) #+ due to subsidies

    # tax_eff_fruit <- 0.15*elasticity_fruit #tax rate in % #15% subsidies
    # if (tax_eff_fruit > 0) tax_eff_fruit <- 0 #negative

    sp$pop[year >= sc_year,
           fruit_curr_xps := as.integer(round(fruit_curr_xps*(1+tax_eff_fruit)),0)] #tax_eff_fruit is negative  #ask Chris whether I can save it as numeric


    # mean = 0.07, sd = ((0.07*1.15) - (0.07*0.85))/3.92)
    # using Betta distribution for percentage, as suggested by Chris:
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
    tax_eff_veg <- rbeta(1, 158.5147, 2102.5560)*(1) #+ due to subsidies

    # tax_eff_veg <- 0.15*elasticity_veg #tax rate in % #15% subsidies
    # if (tax_eff_veg > 0) tax_eff_veg <- 0 #negative

    sp$pop[year >= sc_year,
           veg_curr_xps := as.integer(round(veg_curr_xps*(1+tax_eff_veg)),0)] #tax_eff_veg is negative

    # Restore the original seed state
    .Random.seed <- saved_seed
  }
)

IMPACTncd$
  run(start_runs:n_runs, multicore = TRUE, scenario_nam = "cons_sc1_10tax_15sub")

# print("SUMMARY EXPORTS!!!")
#
# IMPACTncd$export_summaries(multicore = TRUE)



print("SUMMARY EXPORTS NOW!!!")

IMPACTncd$export_summaries(multicore = TRUE)

print("CREATE TABLES NOW!!!")

IMPACTncd$export_tables(
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
)
