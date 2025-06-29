#--------------------------------------------------------------------------------
#---------------------------------------------------------------LD/LQs LOOP

rm(list=ls())

# Load packages
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(ggplot2)
library(ggpubr)
library(shiny)
library(BBmisc)
library(chemCal)
library(openxlsx)
library(writexl)

#------------------------------------------------------------------------------- Load data
load("../Data/Original_Data_Serum_Method.Rdata")
rm(list = setdiff(ls(), c("LDs_LQs_R2_OCs__ME_batch4", "LDs_LQs_R2_OCs_A8A_batch1",
                          "LDs_LQs_R2_OCs_A8A_batch2", "LDs_LQs_R2_OCs_A8A_batch3",
                          "LDs_LQs_R2_OCs_A8A_batch5", "LDs_LQs_R2_OCs_A8A_batch6")))

##Do it only the first time
DB <- NULL 

#------------------------------------------------------------------------------- Batch 1
##Specify all info per batch
##CC1
batch <- 1    #240902_AST8A
d <- LDs_LQs_R2_OCs_A8A_batch1
tbb.inicial <- 52.45
pcb209.inicial <- 50.98
octachloronaphthalene.inicial <- 10.0869

##Calculate the concentrations of standards in samples
##TBB
TBB.conc <- (tbb.inicial * 25)/100
##PCB209
PCB209.conc <- (pcb209.inicial * 25)/100
##ISTD
Octachloronaphthalene.conc <- (octachloronaphthalene.inicial * 100) / 100

##Other details to be processed
##Compounds to be adjusted by TBB
recoveries.tbb <- c("PeCB", "Tecnazene", "a-HCH", "HCB", "b-HCH", "g-HCH",
                    "Quintozene", "d-HCH", "e-HCH", "PCB28", "VIN", "Hepta-Cl",
                    "PCB52", "Aldrin", "Octachlorostyrene", "Isodrin",
                    "B-Hepta-Cl", "Oxy-Chlordane", "A-Hepta-Cl", "TBB")
##Compounds to be adjusted by PCB209
recoveries.pcb209 <- c("Trans-Chlordane", "opDDE", "PCB101",
                       "a-Endosulfan",  "Cis-Chlordane", "ppDDE", "Dieldrin",
                       "opDDD", "Endrin", "PCB118", "b-Endosulfan", "ppDDD",
                       "opDDT", "PCB153", "ppDDT", "PCB138",
                       "Endosulfan-sulfate", "Methoxychlor", "PCB180", "Mirex", "PCB209")

##Defining the parameters to build the M.parameters matrix for the loop
d.interesting <- d %>% filter(!Compound %in% c("TBB", "PCB209", "Octachloronaphthalene")) %>%
  dplyr::select(Compound) %>% unique()
nC <- length(unique(d.interesting$Compound))
compound.labels <- as.character(unique(d.interesting$Compound))

M.parameters <- array(NA, dim=c(nC, 7),
                      dimnames=list(Compound=compound.labels,
                                    parameter=c("R2", "R2-IS",
                                                "LD", "LQ",
                                                "Batch", "Injection Date",
                                                "Compound"))) 
M.parameters


##Loop starts
for(c in 1:nC){
  C <- as.character(unique(d$Compound))[c]
  d.now <- d %>%
    filter(Compound == C) %>%
    droplevels()
  
  #Table with the Blanks and Recoveries (%)
  model.IS <- lm(Area ~ Calibration, data=d %>%
                   filter(Type %in%  "Calibration") %>%
                   filter(!is.na(Calibration)) %>%
                   filter(Compound %in% "Octachloronaphthalene"))
  model.TBB <- lm(Area ~ Calibration, data=d %>%
                    filter(Type %in%  "Calibration") %>%
                    filter(!is.na(Calibration)) %>%
                    filter(Compound %in% "TBB"))
  model.PCB209 <- lm(Area ~ Calibration, data=d %>%
                       filter(Type %in%  "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% "PCB209"))
  
  data.ld <- d %>% 
    filter(Compound %in% c(C, "Octachloronaphthalene")) %>%
    filter(Type %in% "Calibration") %>%
    filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")) %>%
    pivot_wider(names_from="Compound", values_from=c("Area", "Calibration")) %>%
    rename(Area = paste0("Area_", C),
           Calibration = paste0("Calibration_", C)) %>%
    mutate(y = Area / `Area_Octachloronaphthalene`,
           x = Calibration / `Calibration_Octachloronaphthalene`)
  
  #No corrected model by ISTD
  model <- lm(Area ~ Calibration, data=d %>%
                filter(Type %in% "Calibration") %>%
                filter(!is.na(Calibration)) %>%
                filter(Compound %in% C) %>%
                filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")))
  sum <- summary(model)
  print(sum)
  list(sum)
  M.parameters[C, 1] <- sum$adj.r.squared
  
  quan <- d %>%
    filter(Compound %in% c(C, "Octachloronaphthalene", "TBB", "PCB209")) %>%
    mutate(`Injected Conc.`= case_when(
      Compound %in% C ~ (Area - model$coefficients[[1]]) / model$coefficients[[2]],
      Compound %in% "Octachloronaphthalene" ~ (Area - model.IS$coefficients[[1]]) / model.IS$coefficients[[2]],
      Compound %in% "TBB" ~ (Area - model.TBB$coefficients[[1]]) / model.TBB$coefficients[[2]],
      Compound %in% "PCB209" ~ (Area - model.PCB209$coefficients[[1]]) / model.PCB209$coefficients[[2]],
      TRUE ~ 1000000)) %>%
    mutate(`Vial Conc.` = `Injected Conc.` / 2) 
  
  #Corrected model by ISTD (octachloronaphthalene)
  model.ld <- lm(y~x, data = data.ld)
  sum.ld <- summary(model.ld)
  print(sum.ld)
  list(sum.ld)
  M.parameters[C, 2] <- sum.ld$adj.r.squared
  #Now we have computed the adj R2 of the models  
  
  
  #Now we will compute LD/LQ
  IS.correction <- quan %>%
    dplyr::select(Sample, Type, Compound, Vol, `Vial Conc.`) %>%
    spread(Compound, `Vial Conc.`) %>%
    gather(Compound, `Vial Conc.`, -c(Sample, Type, Vol, Octachloronaphthalene)) %>%
    mutate(`Vial Conc. Corrected` = `Vial Conc.` * (Octachloronaphthalene.conc/Octachloronaphthalene)) %>%
    dplyr::select(-`Vial Conc.`)
  Recoveries <- IS.correction %>%
    dplyr::select(-Octachloronaphthalene) %>%
    spread(Compound, `Vial Conc. Corrected`) %>%
    mutate(Recovery = case_when(
      C %in% recoveries.pcb209 ~ (PCB209/PCB209.conc)*100,
      C %in% recoveries.tbb ~ (TBB/TBB.conc)*100,
      TRUE ~ 1000000)) %>%
    dplyr::select(Sample, Vol, Recovery)
  
  t <- IS.correction %>%
    filter(Compound %in% C) %>%
    dplyr::select(-Octachloronaphthalene) %>%
    filter(Type %in% "Calibration") %>%
    left_join(Recoveries) %>%
    mutate(`Sample Conc.` = `Vial Conc. Corrected` * (100 / Vol)) %>%
    mutate(`Sample Conc. with Recovery` = (`Sample Conc.`/Recovery)*100) %>%
    mutate(LD = lod(model.ld)$x * (125/500)/2) %>%
    mutate(LQ = loq(model.ld)$x * (125/500)/2) %>%
    dplyr::select(Sample, Compound, LD, LQ) %>%
    left_join(quan %>% dplyr::select(Sample, Compound, Batch, `Injection Date`)) %>%
    filter(Sample %in% 1) %>%
    unique()
  t
  
  M.parameters[C, 3] <- t$LD
  M.parameters[C, 4] <- t$LQ
  M.parameters[C, 5] <- t$Batch
  M.parameters[C, 6] <- t$`Injection Date`
  M.parameters[C, 7] <- t$Compound
}
M.parameters

##Merge M.parameters to DB
DB <- bind_rows(DB,as.data.frame(M.parameters)) %>%
  dplyr::select(Compound, R2, `R2-IS`, LD, LQ, Batch, `Injection Date`)
DB

#------------------------------------------------------------------------------- Batch 2
##Specify all info per batch
##CC1
batch <- 2.1    #240916_AST8A
d <- LDs_LQs_R2_OCs_A8A_batch2
tbb.inicial <- 52.45
pcb209.inicial <- 50.98
octachloronaphthalene.inicial <- 10.0667

##Calculate the concentrations of standards in samples
##TBB
TBB.conc <- (tbb.inicial * 25)/100
##PCB209
PCB209.conc <- (pcb209.inicial * 25)/100
##ISTD
Octachloronaphthalene.conc <- (octachloronaphthalene.inicial * 100) / 100

##Other details to be processed
##Compounds to be adjusted by TBB
recoveries.tbb <- c("PeCB", "Tecnazene", "a-HCH", "HCB", "b-HCH", "g-HCH",
                    "Quintozene", "d-HCH", "e-HCH", "PCB28", "VIN", "Hepta-Cl",
                    "PCB52", "Aldrin", "Octachlorostyrene", "Isodrin",
                    "B-Hepta-Cl", "Oxy-Chlordane", "A-Hepta-Cl", "TBB")
##Compounds to be adjusted by PCB209
recoveries.pcb209 <- c("Trans-Chlordane", "opDDE", "PCB101",
                       "a-Endosulfan",  "Cis-Chlordane", "ppDDE", "Dieldrin",
                       "opDDD", "Endrin", "PCB118", "b-Endosulfan", "ppDDD",
                       "opDDT", "PCB153", "ppDDT", "PCB138",
                       "Endosulfan-sulfate", "Methoxychlor", "PCB180", "Mirex", "PCB209")

##Defining the parameters to build the M.parameters matrix for the loop
d.interesting <- d %>% filter(!Compound %in% c("TBB", "PCB209", "Octachloronaphthalene")) %>%
  dplyr::select(Compound) %>% unique()
nC <- length(unique(d.interesting$Compound))
compound.labels <- as.character(unique(d.interesting$Compound))

M.parameters <- array(NA, dim=c(nC, 7),
                      dimnames=list(Compound=compound.labels,
                                    parameter=c("R2", "R2-IS",
                                                "LD", "LQ",
                                                "Batch", "Injection Date",
                                                "Compound"))) 
M.parameters


##Loop starts
for(c in 1:nC){
  C <- as.character(unique(d$Compound))[c]
  d.now <- d %>%
    filter(Compound == C) %>%
    droplevels()
  
  #Table with the Blanks and Recoveries (%)
  model.IS <- lm(Area ~ Calibration, data=d %>%
                   filter(Type %in%  "Calibration") %>%
                   filter(!is.na(Calibration)) %>%
                   filter(Compound %in% "Octachloronaphthalene"))
  model.TBB <- lm(Area ~ Calibration, data=d %>%
                    filter(Type %in%  "Calibration") %>%
                    filter(!is.na(Calibration)) %>%
                    filter(Compound %in% "TBB"))
  model.PCB209 <- lm(Area ~ Calibration, data=d %>%
                       filter(Type %in%  "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% "PCB209"))
  
  data.ld <- d %>% 
    filter(Compound %in% c(C, "Octachloronaphthalene")) %>%
    filter(Type %in% "Calibration") %>%
    filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")) %>%
    pivot_wider(names_from="Compound", values_from=c("Area", "Calibration")) %>%
    rename(Area = paste0("Area_", C),
           Calibration = paste0("Calibration_", C)) %>%
    mutate(y = Area / `Area_Octachloronaphthalene`,
           x = Calibration / `Calibration_Octachloronaphthalene`)
  
  #No corrected model by ISTD
  model <- lm(Area ~ Calibration, data=d %>%
                filter(Type %in% "Calibration") %>%
                filter(!is.na(Calibration)) %>%
                filter(Compound %in% C) %>%
                filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")))
  sum <- summary(model)
  print(sum)
  list(sum)
  M.parameters[C, 1] <- sum$adj.r.squared
  
  quan <- d %>%
    filter(Compound %in% c(C, "Octachloronaphthalene", "TBB", "PCB209")) %>%
    mutate(`Injected Conc.`= case_when(
      Compound %in% C ~ (Area - model$coefficients[[1]]) / model$coefficients[[2]],
      Compound %in% "Octachloronaphthalene" ~ (Area - model.IS$coefficients[[1]]) / model.IS$coefficients[[2]],
      Compound %in% "TBB" ~ (Area - model.TBB$coefficients[[1]]) / model.TBB$coefficients[[2]],
      Compound %in% "PCB209" ~ (Area - model.PCB209$coefficients[[1]]) / model.PCB209$coefficients[[2]],
      TRUE ~ 1000000)) %>%
    mutate(`Vial Conc.` = `Injected Conc.` / 2) 
  
  #Corrected model by ISTD (octachloronaphthalene)
  model.ld <- lm(y~x, data = data.ld)
  sum.ld <- summary(model.ld)
  print(sum.ld)
  list(sum.ld)
  M.parameters[C, 2] <- sum.ld$adj.r.squared
  #Now we have computed the adj R2 of the models  
  
  
  #Now we will compute LD/LQ
  IS.correction <- quan %>%
    dplyr::select(Sample, Type, Compound, Vol, `Vial Conc.`) %>%
    spread(Compound, `Vial Conc.`) %>%
    gather(Compound, `Vial Conc.`, -c(Sample, Type, Vol, Octachloronaphthalene)) %>%
    mutate(`Vial Conc. Corrected` = `Vial Conc.` * (Octachloronaphthalene.conc/Octachloronaphthalene)) %>%
    dplyr::select(-`Vial Conc.`)
  Recoveries <- IS.correction %>%
    dplyr::select(-Octachloronaphthalene) %>%
    spread(Compound, `Vial Conc. Corrected`) %>%
    mutate(Recovery = case_when(
      C %in% recoveries.pcb209 ~ (PCB209/PCB209.conc)*100,
      C %in% recoveries.tbb ~ (TBB/TBB.conc)*100,
      TRUE ~ 1000000)) %>%
    dplyr::select(Sample, Vol, Recovery)
  
  t <- IS.correction %>%
    filter(Compound %in% C) %>%
    dplyr::select(-Octachloronaphthalene) %>%
    filter(Type %in% "Calibration") %>%
    left_join(Recoveries) %>%
    mutate(`Sample Conc.` = `Vial Conc. Corrected` * (100 / Vol)) %>%
    mutate(`Sample Conc. with Recovery` = (`Sample Conc.`/Recovery)*100) %>%
    mutate(LD = lod(model.ld)$x * (125/500)/2) %>%
    mutate(LQ = loq(model.ld)$x * (125/500)/2) %>%
    dplyr::select(Sample, Compound, LD, LQ) %>%
    left_join(quan %>% dplyr::select(Sample, Compound, Batch, `Injection Date`)) %>%
    filter(Sample %in% 1) %>%
    unique()
  t
  
  M.parameters[C, 3] <- t$LD
  M.parameters[C, 4] <- t$LQ
  M.parameters[C, 5] <- t$Batch
  M.parameters[C, 6] <- t$`Injection Date`
  M.parameters[C, 7] <- t$Compound
}
M.parameters

##Merge M.parameters to DB
DB <- bind_rows(DB,as.data.frame(M.parameters)) %>%
  dplyr::select(Compound, R2, `R2-IS`, LD, LQ, Batch, `Injection Date`)
DB

#------------------------------------------------------------------------------- Batch 3
##Specify all info per batch
##CC1
batch <- 3.1    #240930_AST8A
d <- LDs_LQs_R2_OCs_A8A_batch3
tbb.inicial <- 54.6762
pcb209.inicial <- 50.6693
octachloronaphthalene.inicial <- 10.0667

##Calculate the concentrations of standards in samples
##TBB
TBB.conc <- (tbb.inicial * 25)/100
##PCB209
PCB209.conc <- (pcb209.inicial * 25)/100
##ISTD
Octachloronaphthalene.conc <- (octachloronaphthalene.inicial * 100) / 100

##Other details to be processed
##Compounds to be adjusted by TBB
recoveries.tbb <- c("PeCB", "Tecnazene", "a-HCH", "HCB", "b-HCH", "g-HCH",
                    "Quintozene", "d-HCH", "e-HCH", "PCB28", "VIN", "Hepta-Cl",
                    "PCB52", "Aldrin", "Octachlorostyrene", "Isodrin",
                    "B-Hepta-Cl", "Oxy-Chlordane", "A-Hepta-Cl", "TBB")
##Compounds to be adjusted by PCB209
recoveries.pcb209 <- c("Trans-Chlordane", "opDDE", "PCB101",
                       "a-Endosulfan",  "Cis-Chlordane", "ppDDE", "Dieldrin",
                       "opDDD", "Endrin", "PCB118", "b-Endosulfan", "ppDDD",
                       "opDDT", "PCB153", "ppDDT", "PCB138",
                       "Endosulfan-sulfate", "Methoxychlor", "PCB180", "Mirex", "PCB209")

##Defining the parameters to build the M.parameters matrix for the loop
d.interesting <- d %>% filter(!Compound %in% c("TBB", "PCB209", "Octachloronaphthalene")) %>%
  dplyr::select(Compound) %>% unique()
nC <- length(unique(d.interesting$Compound))
compound.labels <- as.character(unique(d.interesting$Compound))

M.parameters <- array(NA, dim=c(nC, 7),
                      dimnames=list(Compound=compound.labels,
                                    parameter=c("R2", "R2-IS",
                                                "LD", "LQ",
                                                "Batch", "Injection Date",
                                                "Compound"))) 
M.parameters


##Loop starts
for(c in 1:nC){
  C <- as.character(unique(d$Compound))[c]
  d.now <- d %>%
    filter(Compound == C) %>%
    droplevels()
  
  #Table with the Blanks and Recoveries (%)
  model.IS <- lm(Area ~ Calibration, data=d %>%
                   filter(Type %in%  "Calibration") %>%
                   filter(!is.na(Calibration)) %>%
                   filter(Compound %in% "Octachloronaphthalene"))
  model.TBB <- lm(Area ~ Calibration, data=d %>%
                    filter(Type %in%  "Calibration") %>%
                    filter(!is.na(Calibration)) %>%
                    filter(Compound %in% "TBB"))
  model.PCB209 <- lm(Area ~ Calibration, data=d %>%
                       filter(Type %in%  "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% "PCB209"))
  
  data.ld <- d %>% 
    filter(Compound %in% c(C, "Octachloronaphthalene")) %>%
    filter(Type %in% "Calibration") %>%
    filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")) %>%
    pivot_wider(names_from="Compound", values_from=c("Area", "Calibration")) %>%
    rename(Area = paste0("Area_", C),
           Calibration = paste0("Calibration_", C)) %>%
    mutate(y = Area / `Area_Octachloronaphthalene`,
           x = Calibration / `Calibration_Octachloronaphthalene`)
  
  #No corrected model by ISTD
  model <- lm(Area ~ Calibration, data=d %>%
                filter(Type %in% "Calibration") %>%
                filter(!is.na(Calibration)) %>%
                filter(Compound %in% C) %>%
                filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")))
  sum <- summary(model)
  print(sum)
  list(sum)
  M.parameters[C, 1] <- sum$adj.r.squared
  
  quan <- d %>%
    filter(Compound %in% c(C, "Octachloronaphthalene", "TBB", "PCB209")) %>%
    mutate(`Injected Conc.`= case_when(
      Compound %in% C ~ (Area - model$coefficients[[1]]) / model$coefficients[[2]],
      Compound %in% "Octachloronaphthalene" ~ (Area - model.IS$coefficients[[1]]) / model.IS$coefficients[[2]],
      Compound %in% "TBB" ~ (Area - model.TBB$coefficients[[1]]) / model.TBB$coefficients[[2]],
      Compound %in% "PCB209" ~ (Area - model.PCB209$coefficients[[1]]) / model.PCB209$coefficients[[2]],
      TRUE ~ 1000000)) %>%
    mutate(`Vial Conc.` = `Injected Conc.` / 2) 
  
  #Corrected model by ISTD (octachloronaphthalene)
  model.ld <- lm(y~x, data = data.ld)
  sum.ld <- summary(model.ld)
  print(sum.ld)
  list(sum.ld)
  M.parameters[C, 2] <- sum.ld$adj.r.squared
  #Now we have computed the adj R2 of the models  
  
  
  #Now we will compute LD/LQ
  IS.correction <- quan %>%
    dplyr::select(Sample, Type, Compound, Vol, `Vial Conc.`) %>%
    spread(Compound, `Vial Conc.`) %>%
    gather(Compound, `Vial Conc.`, -c(Sample, Type, Vol, Octachloronaphthalene)) %>%
    mutate(`Vial Conc. Corrected` = `Vial Conc.` * (Octachloronaphthalene.conc/Octachloronaphthalene)) %>%
    dplyr::select(-`Vial Conc.`)
  Recoveries <- IS.correction %>%
    dplyr::select(-Octachloronaphthalene) %>%
    spread(Compound, `Vial Conc. Corrected`) %>%
    mutate(Recovery = case_when(
      C %in% recoveries.pcb209 ~ (PCB209/PCB209.conc)*100,
      C %in% recoveries.tbb ~ (TBB/TBB.conc)*100,
      TRUE ~ 1000000)) %>%
    dplyr::select(Sample, Vol, Recovery)
  
  t <- IS.correction %>%
    filter(Compound %in% C) %>%
    dplyr::select(-Octachloronaphthalene) %>%
    filter(Type %in% "Calibration") %>%
    left_join(Recoveries) %>%
    mutate(`Sample Conc.` = `Vial Conc. Corrected` * (100 / Vol)) %>%
    mutate(`Sample Conc. with Recovery` = (`Sample Conc.`/Recovery)*100) %>%
    mutate(LD = lod(model.ld)$x * (125/500)/2) %>%
    mutate(LQ = loq(model.ld)$x * (125/500)/2) %>%
    dplyr::select(Sample, Compound, LD, LQ) %>%
    left_join(quan %>% dplyr::select(Sample, Compound, Batch, `Injection Date`)) %>%
    filter(Sample %in% 1) %>%
    unique()
  t
  
  M.parameters[C, 3] <- t$LD
  M.parameters[C, 4] <- t$LQ
  M.parameters[C, 5] <- t$Batch
  M.parameters[C, 6] <- t$`Injection Date`
  M.parameters[C, 7] <- t$Compound
}
M.parameters

##Merge M.parameters to DB
DB <- bind_rows(DB,as.data.frame(M.parameters)) %>%
  dplyr::select(Compound, R2, `R2-IS`, LD, LQ, Batch, `Injection Date`)
DB

#------------------------------------------------------------------------------- Batch 4
##Specify all info per batch
##CC2
batch <- 4    #250429_ME_OCs
d <- LDs_LQs_R2_OCs__ME_batch4
tbb.inicial <- 54.8409
pcb209.inicial <- 50.8220
octachloronaphthalene.inicial <- 9.9887


##Calculate the concentrations of standards in samples
##TBB
TBB.conc <- (tbb.inicial * 25)/100
##PCB209
PCB209.conc <- (pcb209.inicial * 25)/100
##ISTD
Octachloronaphthalene.conc <- (octachloronaphthalene.inicial * 100) / 100

##Other details to be processed
##Compounds to be adjusted by TBB
recoveries.tbb <- c("PeCB", "Tecnazene", "a-HCH", "HCB", "b-HCH", "g-HCH",
                    "Quintozene", "d-HCH", "e-HCH", "PCB28", "VIN", "Hepta-Cl",
                    "PCB52", "Aldrin", "Octachlorostyrene", "Isodrin",
                    "B-Hepta-Cl", "Oxy-Chlordane", "A-Hepta-Cl", "TBB")
##Compounds to be adjusted by PCB209
recoveries.pcb209 <- c("Trans-Chlordane", "opDDE", "PCB101",
                       "a-Endosulfan",  "Cis-Chlordane", "ppDDE", "Dieldrin",
                       "opDDD", "Endrin", "PCB118", "b-Endosulfan", "ppDDD",
                       "opDDT", "PCB153", "ppDDT", "PCB138",
                       "Endosulfan-sulfate", "Methoxychlor", "PCB180", "Mirex", "PCB209")

##Defining the parameters to build the M.parameters matrix for the loop
d.interesting <- d %>% filter(!Compound %in% c("TBB", "PCB209", "Octachloronaphthalene")) %>%
  dplyr::select(Compound) %>% unique()
nC <- length(unique(d.interesting$Compound))
compound.labels <- as.character(unique(d.interesting$Compound))

M.parameters <- array(NA, dim=c(nC, 7),
                      dimnames=list(Compound=compound.labels,
                                    parameter=c("R2", "R2-IS",
                                                "LD", "LQ",
                                                "Batch", "Injection Date",
                                                "Compound"))) 
M.parameters


##Loop starts
for(c in 1:nC){
  C <- as.character(unique(d$Compound))[c]
  d.now <- d %>%
    filter(Compound == C) %>%
    droplevels()

  #Table with the Blanks and Recoveries (%)
  model.IS <- lm(Area ~ Calibration, data=d %>%
                   filter(Type %in%  "Calibration") %>%
                   filter(!is.na(Calibration)) %>%
                   filter(Compound %in% "Octachloronaphthalene"))
  model.TBB <- lm(Area ~ Calibration, data=d %>%
                    filter(Type %in%  "Calibration") %>%
                    filter(!is.na(Calibration)) %>%
                    filter(Compound %in% "TBB"))
  model.PCB209 <- lm(Area ~ Calibration, data=d %>%
                       filter(Type %in%  "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% "PCB209"))
  
  data.ld <- d %>% 
    filter(Compound %in% c(C, "Octachloronaphthalene")) %>%
    filter(Type %in% "Calibration") %>%
    filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")) %>%
    pivot_wider(names_from="Compound", values_from=c("Area", "Calibration")) %>%
    rename(Area = paste0("Area_", C),
           Calibration = paste0("Calibration_", C)) %>%
    mutate(y = Area / `Area_Octachloronaphthalene`,
           x = Calibration / `Calibration_Octachloronaphthalene`)
  
  #No corrected model by ISTD
  model <- lm(Area ~ Calibration, data=d %>%
                filter(Type %in% "Calibration") %>%
                filter(!is.na(Calibration)) %>%
                filter(Compound %in% C) %>%
                filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")))
  sum <- summary(model)
  print(sum)
  list(sum)
  M.parameters[C, 1] <- sum$adj.r.squared
  
  quan <- d %>%
    filter(Compound %in% c(C, "Octachloronaphthalene", "TBB", "PCB209")) %>%
    mutate(`Injected Conc.`= case_when(
      Compound %in% C ~ (Area - model$coefficients[[1]]) / model$coefficients[[2]],
      Compound %in% "Octachloronaphthalene" ~ (Area - model.IS$coefficients[[1]]) / model.IS$coefficients[[2]],
      Compound %in% "TBB" ~ (Area - model.TBB$coefficients[[1]]) / model.TBB$coefficients[[2]],
      Compound %in% "PCB209" ~ (Area - model.PCB209$coefficients[[1]]) / model.PCB209$coefficients[[2]],
      TRUE ~ 1000000)) %>%
    mutate(`Vial Conc.` = `Injected Conc.` / 2) 
  
  #Corrected model by ISTD (octachloronaphthalene)
  model.ld <- lm(y~x, data = data.ld)
  sum.ld <- summary(model.ld)
  print(sum.ld)
  list(sum.ld)
  M.parameters[C, 2] <- sum.ld$adj.r.squared
  #Now we have computed the adj R2 of the models  
  
  
  #Now we will compute LD/LQ
  IS.correction <- quan %>%
    dplyr::select(Sample, Type, Compound, Vol, `Vial Conc.`) %>%
    spread(Compound, `Vial Conc.`) %>%
    gather(Compound, `Vial Conc.`, -c(Sample, Type, Vol, Octachloronaphthalene)) %>%
    mutate(`Vial Conc. Corrected` = `Vial Conc.` * (Octachloronaphthalene.conc/Octachloronaphthalene)) %>%
    dplyr::select(-`Vial Conc.`)
  Recoveries <- IS.correction %>%
    dplyr::select(-Octachloronaphthalene) %>%
    spread(Compound, `Vial Conc. Corrected`) %>%
    mutate(Recovery = case_when(
      C %in% recoveries.pcb209 ~ (PCB209/PCB209.conc)*100,
      C %in% recoveries.tbb ~ (TBB/TBB.conc)*100,
      TRUE ~ 1000000)) %>%
    dplyr::select(Sample, Vol, Recovery)
  
  t <- IS.correction %>%
    filter(Compound %in% C) %>%
    dplyr::select(-Octachloronaphthalene) %>%
    filter(Type %in% "Calibration") %>%
    left_join(Recoveries) %>%
    mutate(`Sample Conc.` = `Vial Conc. Corrected` * (100 / Vol)) %>%
    mutate(`Sample Conc. with Recovery` = (`Sample Conc.`/Recovery)*100) %>%
    mutate(LD = lod(model.ld)$x * (125/500)/2) %>%
    mutate(LQ = loq(model.ld)$x * (125/500)/2) %>%
    dplyr::select(Sample, Compound, LD, LQ) %>%
    left_join(quan %>% dplyr::select(Sample, Compound, Batch, `Injection Date`)) %>%
    filter(Sample %in% 1) %>%
    unique()
  t
  
  M.parameters[C, 3] <- t$LD
  M.parameters[C, 4] <- t$LQ
  M.parameters[C, 5] <- t$Batch
  M.parameters[C, 6] <- t$`Injection Date`
  M.parameters[C, 7] <- t$Compound
}
M.parameters

##Merge M.parameters to DB
DB <- bind_rows(DB,as.data.frame(M.parameters)) %>%
  dplyr::select(Compound, R2, `R2-IS`, LD, LQ, Batch, `Injection Date`)
DB


#------------------------------------------------------------------------------- Batch 5
##Specify all info per batch
##CC1
batch <- 5    #250210_AST8A
d <- LDs_LQs_R2_OCs_A8A_batch5
tbb.inicial <- 55.8358
pcb209.inicial <- 51.7439
octachloronaphthalene.inicial <- 9.9683

##Calculate the concentrations of standards in samples
##TBB
TBB.conc <- (tbb.inicial * 25)/100
##PCB209
PCB209.conc <- (pcb209.inicial * 25)/100
##ISTD
Octachloronaphthalene.conc <- (octachloronaphthalene.inicial * 100) / 100

##Other details to be processed
##Compounds to be adjusted by TBB
recoveries.tbb <- c("PeCB", "Tecnazene", "a-HCH", "HCB", "b-HCH", "g-HCH",
                    "Quintozene", "d-HCH", "e-HCH", "PCB28", "VIN", "Hepta-Cl",
                    "PCB52", "Aldrin", "Octachlorostyrene", "Isodrin",
                    "B-Hepta-Cl", "Oxy-Chlordane", "A-Hepta-Cl", "TBB")
##Compounds to be adjusted by PCB209
recoveries.pcb209 <- c("Trans-Chlordane", "opDDE", "PCB101",
                       "a-Endosulfan",  "Cis-Chlordane", "ppDDE", "Dieldrin",
                       "opDDD", "Endrin", "PCB118", "b-Endosulfan", "ppDDD",
                       "opDDT", "PCB153", "ppDDT", "PCB138",
                       "Endosulfan-sulfate", "Methoxychlor", "PCB180", "Mirex", "PCB209")

##Defining the parameters to build the M.parameters matrix for the loop
d.interesting <- d %>% filter(!Compound %in% c("TBB", "PCB209", "Octachloronaphthalene")) %>%
  dplyr::select(Compound) %>% unique()
nC <- length(unique(d.interesting$Compound))
compound.labels <- as.character(unique(d.interesting$Compound))

M.parameters <- array(NA, dim=c(nC, 7),
                      dimnames=list(Compound=compound.labels,
                                    parameter=c("R2", "R2-IS",
                                                "LD", "LQ",
                                                "Batch", "Injection Date",
                                                "Compound"))) 
M.parameters


##Loop starts
for(c in 1:nC){
  C <- as.character(unique(d$Compound))[c]
  d.now <- d %>%
    filter(Compound == C) %>%
    droplevels()
  
  #Table with the Blanks and Recoveries (%)
  model.IS <- lm(Area ~ Calibration, data=d %>%
                   filter(Type %in%  "Calibration") %>%
                   filter(!is.na(Calibration)) %>%
                   filter(Compound %in% "Octachloronaphthalene"))
  model.TBB <- lm(Area ~ Calibration, data=d %>%
                    filter(Type %in%  "Calibration") %>%
                    filter(!is.na(Calibration)) %>%
                    filter(Compound %in% "TBB"))
  model.PCB209 <- lm(Area ~ Calibration, data=d %>%
                       filter(Type %in%  "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% "PCB209"))
  
  data.ld <- d %>% 
    filter(Compound %in% c(C, "Octachloronaphthalene")) %>%
    filter(Type %in% "Calibration") %>%
    filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")) %>%
    pivot_wider(names_from="Compound", values_from=c("Area", "Calibration")) %>%
    rename(Area = paste0("Area_", C),
           Calibration = paste0("Calibration_", C)) %>%
    mutate(y = Area / `Area_Octachloronaphthalene`,
           x = Calibration / `Calibration_Octachloronaphthalene`)
  
  #No corrected model by ISTD
  model <- lm(Area ~ Calibration, data=d %>%
                filter(Type %in% "Calibration") %>%
                filter(!is.na(Calibration)) %>%
                filter(Compound %in% C) %>%
                filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")))
  sum <- summary(model)
  print(sum)
  list(sum)
  M.parameters[C, 1] <- sum$adj.r.squared
  
  quan <- d %>%
    filter(Compound %in% c(C, "Octachloronaphthalene", "TBB", "PCB209")) %>%
    mutate(`Injected Conc.`= case_when(
      Compound %in% C ~ (Area - model$coefficients[[1]]) / model$coefficients[[2]],
      Compound %in% "Octachloronaphthalene" ~ (Area - model.IS$coefficients[[1]]) / model.IS$coefficients[[2]],
      Compound %in% "TBB" ~ (Area - model.TBB$coefficients[[1]]) / model.TBB$coefficients[[2]],
      Compound %in% "PCB209" ~ (Area - model.PCB209$coefficients[[1]]) / model.PCB209$coefficients[[2]],
      TRUE ~ 1000000)) %>%
    mutate(`Vial Conc.` = `Injected Conc.` / 2) 
  
  #Corrected model by ISTD (octachloronaphthalene)
  model.ld <- lm(y~x, data = data.ld)
  sum.ld <- summary(model.ld)
  print(sum.ld)
  list(sum.ld)
  M.parameters[C, 2] <- sum.ld$adj.r.squared
  #Now we have computed the adj R2 of the models  
  
  
  #Now we will compute LD/LQ
  IS.correction <- quan %>%
    dplyr::select(Sample, Type, Compound, Vol, `Vial Conc.`) %>%
    spread(Compound, `Vial Conc.`) %>%
    gather(Compound, `Vial Conc.`, -c(Sample, Type, Vol, Octachloronaphthalene)) %>%
    mutate(`Vial Conc. Corrected` = `Vial Conc.` * (Octachloronaphthalene.conc/Octachloronaphthalene)) %>%
    dplyr::select(-`Vial Conc.`)
  Recoveries <- IS.correction %>%
    dplyr::select(-Octachloronaphthalene) %>%
    spread(Compound, `Vial Conc. Corrected`) %>%
    mutate(Recovery = case_when(
      C %in% recoveries.pcb209 ~ (PCB209/PCB209.conc)*100,
      C %in% recoveries.tbb ~ (TBB/TBB.conc)*100,
      TRUE ~ 1000000)) %>%
    dplyr::select(Sample, Vol, Recovery)
  
  t <- IS.correction %>%
    filter(Compound %in% C) %>%
    dplyr::select(-Octachloronaphthalene) %>%
    filter(Type %in% "Calibration") %>%
    left_join(Recoveries) %>%
    mutate(`Sample Conc.` = `Vial Conc. Corrected` * (100 / Vol)) %>%
    mutate(`Sample Conc. with Recovery` = (`Sample Conc.`/Recovery)*100) %>%
    mutate(LD = lod(model.ld)$x * (125/500)/2) %>%
    mutate(LQ = loq(model.ld)$x * (125/500)/2) %>%
    dplyr::select(Sample, Compound, LD, LQ) %>%
    left_join(quan %>% dplyr::select(Sample, Compound, Batch, `Injection Date`)) %>%
    filter(Sample %in% 1) %>%
    unique()
  t
  
  M.parameters[C, 3] <- t$LD
  M.parameters[C, 4] <- t$LQ
  M.parameters[C, 5] <- t$Batch
  M.parameters[C, 6] <- t$`Injection Date`
  M.parameters[C, 7] <- t$Compound
}
M.parameters

##Merge M.parameters to DB
DB <- bind_rows(DB,as.data.frame(M.parameters)) %>%
  dplyr::select(Compound, R2, `R2-IS`, LD, LQ, Batch, `Injection Date`)
DB


#------------------------------------------------------------------------------- Batch 6
##Specify all info per batch
##CC1
batch <- 6    #250217_AST8A
d <- LDs_LQs_R2_OCs_A8A_batch6
tbb.inicial <- 55.8358
pcb209.inicial <- 51.7439
octachloronaphthalene.inicial <- 9.9683

##Calculate the concentrations of standards in samples
##TBB
TBB.conc <- (tbb.inicial * 25)/100
##PCB209
PCB209.conc <- (pcb209.inicial * 25)/100
##ISTD
Octachloronaphthalene.conc <- (octachloronaphthalene.inicial * 100) / 100

##Other details to be processed
##Compounds to be adjusted by TBB
recoveries.tbb <- c("PeCB", "Tecnazene", "a-HCH", "HCB", "b-HCH", "g-HCH",
                    "Quintozene", "d-HCH", "e-HCH", "PCB28", "VIN", "Hepta-Cl",
                    "PCB52", "Aldrin", "Octachlorostyrene", "Isodrin",
                    "B-Hepta-Cl", "Oxy-Chlordane", "A-Hepta-Cl", "TBB")
##Compounds to be adjusted by PCB209
recoveries.pcb209 <- c("Trans-Chlordane", "opDDE", "PCB101",
                       "a-Endosulfan",  "Cis-Chlordane", "ppDDE", "Dieldrin",
                       "opDDD", "Endrin", "PCB118", "b-Endosulfan", "ppDDD",
                       "opDDT", "PCB153", "ppDDT", "PCB138",
                       "Endosulfan-sulfate", "Methoxychlor", "PCB180", "Mirex", "PCB209")

##Defining the parameters to build the M.parameters matrix for the loop
d.interesting <- d %>% filter(!Compound %in% c("TBB", "PCB209", "Octachloronaphthalene")) %>%
  dplyr::select(Compound) %>% unique()
nC <- length(unique(d.interesting$Compound))
compound.labels <- as.character(unique(d.interesting$Compound))

M.parameters <- array(NA, dim=c(nC, 7),
                      dimnames=list(Compound=compound.labels,
                                    parameter=c("R2", "R2-IS",
                                                "LD", "LQ",
                                                "Batch", "Injection Date",
                                                "Compound"))) 
M.parameters


##Loop starts
for(c in 1:nC){
  C <- as.character(unique(d$Compound))[c]
  d.now <- d %>%
    filter(Compound == C) %>%
    droplevels()
  
  #Table with the Blanks and Recoveries (%)
  model.IS <- lm(Area ~ Calibration, data=d %>%
                   filter(Type %in%  "Calibration") %>%
                   filter(!is.na(Calibration)) %>%
                   filter(Compound %in% "Octachloronaphthalene"))
  model.TBB <- lm(Area ~ Calibration, data=d %>%
                    filter(Type %in%  "Calibration") %>%
                    filter(!is.na(Calibration)) %>%
                    filter(Compound %in% "TBB"))
  model.PCB209 <- lm(Area ~ Calibration, data=d %>%
                       filter(Type %in%  "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% "PCB209"))
  
  data.ld <- d %>% 
    filter(Compound %in% c(C, "Octachloronaphthalene")) %>%
    filter(Type %in% "Calibration") %>%
    filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")) %>%
    pivot_wider(names_from="Compound", values_from=c("Area", "Calibration")) %>%
    rename(Area = paste0("Area_", C),
           Calibration = paste0("Calibration_", C)) %>%
    mutate(y = Area / `Area_Octachloronaphthalene`,
           x = Calibration / `Calibration_Octachloronaphthalene`)
  
  #No corrected model by ISTD
  model <- lm(Area ~ Calibration, data=d %>%
                filter(Type %in% "Calibration") %>%
                filter(!is.na(Calibration)) %>%
                filter(Compound %in% C) %>%
                filter(Sample %in% c("0.06", "0.12", "0.25", "0.5", "1", "3", "5", "7", "12")))
  sum <- summary(model)
  print(sum)
  list(sum)
  M.parameters[C, 1] <- sum$adj.r.squared
  
  quan <- d %>%
    filter(Compound %in% c(C, "Octachloronaphthalene", "TBB", "PCB209")) %>%
    mutate(`Injected Conc.`= case_when(
      Compound %in% C ~ (Area - model$coefficients[[1]]) / model$coefficients[[2]],
      Compound %in% "Octachloronaphthalene" ~ (Area - model.IS$coefficients[[1]]) / model.IS$coefficients[[2]],
      Compound %in% "TBB" ~ (Area - model.TBB$coefficients[[1]]) / model.TBB$coefficients[[2]],
      Compound %in% "PCB209" ~ (Area - model.PCB209$coefficients[[1]]) / model.PCB209$coefficients[[2]],
      TRUE ~ 1000000)) %>%
    mutate(`Vial Conc.` = `Injected Conc.` / 2) 
  
  #Corrected model by ISTD (octachloronaphthalene)
  model.ld <- lm(y~x, data = data.ld)
  sum.ld <- summary(model.ld)
  print(sum.ld)
  list(sum.ld)
  M.parameters[C, 2] <- sum.ld$adj.r.squared
  #Now we have computed the adj R2 of the models  
  
  
  #Now we will compute LD/LQ
  IS.correction <- quan %>%
    dplyr::select(Sample, Type, Compound, Vol, `Vial Conc.`) %>%
    spread(Compound, `Vial Conc.`) %>%
    gather(Compound, `Vial Conc.`, -c(Sample, Type, Vol, Octachloronaphthalene)) %>%
    mutate(`Vial Conc. Corrected` = `Vial Conc.` * (Octachloronaphthalene.conc/Octachloronaphthalene)) %>%
    dplyr::select(-`Vial Conc.`)
  Recoveries <- IS.correction %>%
    dplyr::select(-Octachloronaphthalene) %>%
    spread(Compound, `Vial Conc. Corrected`) %>%
    mutate(Recovery = case_when(
      C %in% recoveries.pcb209 ~ (PCB209/PCB209.conc)*100,
      C %in% recoveries.tbb ~ (TBB/TBB.conc)*100,
      TRUE ~ 1000000)) %>%
    dplyr::select(Sample, Vol, Recovery)
  
  t <- IS.correction %>%
    filter(Compound %in% C) %>%
    dplyr::select(-Octachloronaphthalene) %>%
    filter(Type %in% "Calibration") %>%
    left_join(Recoveries) %>%
    mutate(`Sample Conc.` = `Vial Conc. Corrected` * (100 / Vol)) %>%
    mutate(`Sample Conc. with Recovery` = (`Sample Conc.`/Recovery)*100) %>%
    mutate(LD = lod(model.ld)$x * (125/500)/2) %>%
    mutate(LQ = loq(model.ld)$x * (125/500)/2) %>%
    dplyr::select(Sample, Compound, LD, LQ) %>%
    left_join(quan %>% dplyr::select(Sample, Compound, Batch, `Injection Date`)) %>%
    filter(Sample %in% 1) %>%
    unique()
  t
  
  M.parameters[C, 3] <- t$LD
  M.parameters[C, 4] <- t$LQ
  M.parameters[C, 5] <- t$Batch
  M.parameters[C, 6] <- t$`Injection Date`
  M.parameters[C, 7] <- t$Compound
}
M.parameters

##Merge M.parameters to DB
DB <- bind_rows(DB,as.data.frame(M.parameters)) %>%
  dplyr::select(Compound, R2, `R2-IS`, LD, LQ, Batch, `Injection Date`)
DB

##------------------------------------------------------------------------------ FINAL DB
##Making some changes to variables in DB (data base) and deleting the first row
DB <- DB %>%
    mutate(Compound = as.factor(Compound),
      R2 = as.numeric(R2),
      `R2-IS` = as.numeric(`R2-IS`),
      LD = as.numeric(LD),
      LQ = as.numeric(LQ),
      Batch = as.factor(Batch),
      `Injection Date` = as.factor(`Injection Date`))
rownames(DB) <- NULL
DB
#write.xlsx(DB, file="DB.xlsx")


##Cheking there are no LDs > LQs
DB %>% 
  tibble() %>%
  filter(LD > LQ) 

#DB <- read.xlsx("DB.xlsx")
##Calculation of means of R2, R2-IS, LD and LQ
DB_mean <- DB %>%
#filter(!Batch == "6") %>%
  group_by(Compound) %>%
  summarize(R2 = mean(R2),
         `R2-IS` = mean(`R2-IS`),
         LD = mean(LD),
         LQ = mean(LQ)) %>%
  unique()

##Setting the right compound order on DB_mean data base
Order <- c("PeCB", "HCB", "Octachlorostyrene",
           "a-HCH", "b-HCH", "g-HCH", "d-HCH", "e-HCH",
           "Hepta-Cl", "A-Hepta-Cl", "B-Hepta-Cl",
           "Aldrin", "Endrin", "Dieldrin", "Isodrin",
           "Oxy-Chlordane", "Trans-Chlordane", "Cis-Chlordane",
           "a-Endosulfan", "b-Endosulfan", "Endosulfan-sulfate", "Mirex",
           "PCB28", "PCB52", "PCB101", "PCB118", "PCB153", "PCB138", "PCB180",
           "opDDE", "ppDDE", "opDDD", "ppDDD", "opDDT", "ppDDT", "Methoxychlor",
           "Tecnazene", "Quintozene", "VIN")
Table1_1 <- DB_mean%>%
  arrange(factor(Compound, levels = Order)) %>%
  dplyr::select(Compound, LD, LQ, `R2-IS`) %>%
  rename(Linearity = `R2-IS`) %>%
  mutate(LD = round(LD, 4),
         LQ = round(LQ, 4),
         Linearity = round(Linearity, 4))

Table1_1
write.xlsx(Table1_1, file="Table1_1.xlsx")

