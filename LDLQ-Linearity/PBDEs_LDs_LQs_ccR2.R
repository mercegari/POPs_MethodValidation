rm(list=ls())

# Load packages
library(shiny)
library(BBmisc)
library(chemCal)
library(writexl)
library(stringr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(readxl)
library(dplyr)
library(stats)
library(tidyr)
library(ggrepel)
library(kableExtra)
library(formattable)
library(forcats)
library(bookdown)
library(patchwork)
library(scales)
library(psych)
library(see)
library(openxlsx)

#Load data
load("../Data/Original_Data_Serum_Method.Rdata")
rm(list = setdiff(ls(), c("d.BDE118_T1", "d.BDE209c_T1", "d.BDE118_T2", "d.BDE209c_T2")))

# Only do it once at the beginning
DB <- NULL 

###----LDs/LQs
##CC1
#
batch <- 1   #250407_LD_LQ_PBDEs
d.BDE118 <- d.BDE118_T1
d.BDE209c <- d.BDE209c_T1
pcb209.inicial <- 50.822
bde118.inicial <- 10.03678088
bde209c.inicial <- 9.897961209

#----
## Calculate the concentrations of standards in samples
##PCB209
PCB209.conc <- (pcb209.inicial * 25)/50
##BDE118
BDE118.conc <- (bde118.inicial * 50)/50
##BDE209c
BDE209c.conc <- (bde209c.inicial * 50)/50 
  
#Defining the surrogate for injection correction (PCB209)
recoveries.pcb209 <- c("BDE28",   "BDE47",   "BDE100",  "BDE99",   "BDE154",  "BDE153", 
                         "BDE183",  "BDE209",  "BDE118",  "BDE209c", "PCB209")

#Defining a matrix for the compounds thath are corrected by BDE118
d.interesting1 <- d.BDE118 %>%
  filter(!Compound %in% c("PCB209", "BDE118", "BDE209c")) %>%
  dplyr::select(Compound) %>% unique()
nC <- length(unique(d.interesting1$Compound))
compound.labels1 <- as.character(unique(d.interesting1$Compound))

M.parameters1 <- array(NA, dim=c(nC, 7),
                       dimnames=list(Compound=compound.labels1,
                                     parameter=c("R2", "R2-IS",
                                                 "LD", "LQ",
                                                 "Batch", "Injection Date",
                                                 "Compound"))) 
M.parameters1

#Defining a matrix for the compounds thath are corrected by BDE209c
d.interesting2 <- d.BDE209c %>%
  filter(!Compound %in% c("PCB209", "BDE118", "BDE209c")) %>%
  dplyr::select(Compound) %>% unique()
nC2 <- length(unique(d.interesting2$Compound))
compound.labels2 <- as.character(unique(d.interesting2$Compound))

M.parameters2 <- array(NA, dim=c(nC2, 7),
                       dimnames=list(Compound=compound.labels2,
                                     parameter=c("R2", "R2-IS",
                                                 "LD", "LQ",
                                                 "Batch", "Injection Date",
                                                 "Compound"))) 
M.parameters2



## Loop for compounds corrected by BDE118 in batch nº1 starts

for(c in 1:nC){
  C <- as.character(unique(d.BDE118$Compound))[c]
  d.now <- d.BDE118 %>%
    filter(Compound == C) %>%
    droplevels()
  
  # Non-corrected models
  model.PCB209 <- lm(Area ~ Calibration, data=d.BDE118 %>%
                       filter(Type %in%  "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% "PCB209"))
  model.BDE118 <- lm(Area ~ Calibration, data=d.BDE118 %>%
                       filter(Type %in%  "Calibration") %>%
                       filter(Compound %in% "BDE118"))
  
  data.ld.BDE118 <- d.BDE118 %>%
    filter(Compound %in% c(C, "BDE118")) %>%
    filter(Type %in% "Calibration") %>%
    filter(Sample %in% c("0.02", "0.04", "0.16", "0.4", "1.6", "3.33", "8.33", "16", "33")) %>%
    pivot_wider(names_from="Compound", values_from=c("Area", "Calibration")) %>%
    rename(Area = paste0("Area_", C),
           Calibration = paste0("Calibration_", C)) %>%
    mutate(y = Area / `Area_BDE118`,
           x = Calibration / `Calibration_BDE118`)
  
  
  #No corrected model by ISTD
  model.BDE118 <- lm(Area ~ Calibration, data=d.BDE118 %>%
                       filter(Type %in% "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% C) %>%
                       filter(Sample %in% c("0.02", "0.04", "0.16", "0.4", "1.6", "3.33", "8.33", "16", "33")))
  sum.BDE118 <- summary(model.BDE118)
  print(sum.BDE118)
  list(sum.BDE118)
  M.parameters1[C, 1] <- sum.BDE118$adj.r.squared
  
  
  quan.BDE118 <- d.BDE118 %>%
    filter(Compound %in% c(C, "BDE118", "PCB209")) %>%
    mutate(`Injected Conc.`= case_when(
      Compound %in% C ~ (Area - model.BDE118$coefficients[[1]]) / model.BDE118$coefficients[[2]],
      Compound %in% "BDE118" ~ (Area - model.BDE118$coefficients[[1]]) / model.BDE118$coefficients[[2]],
      Compound %in% "PCB209" ~ (Area - model.PCB209$coefficients[[1]]) / model.PCB209$coefficients[[2]],
      TRUE ~ 1000000)) %>%
    mutate(`Vial Conc.` = `Injected Conc.` / 2) 
  
  
  #Corrected model by BDE118 or BDE209c
  model.ld.BDE118 <- lm(y~x, data = data.ld.BDE118)
  sum.ld.BDE118 <- summary(model.ld.BDE118)
  print(sum.ld.BDE118)
  list(sum.ld.BDE118)
  M.parameters1[C, 2] <- sum.ld.BDE118$adj.r.squared
  # Now we have computed the adj R2 of the models  
  
  
  # Now we will compute LD/LQ
  IS.correction.BDE118 <- quan.BDE118 %>%
    dplyr::select(Sample, Type, Compound, `Vial Conc.`) %>%
    spread(Compound, `Vial Conc.`) %>%
    gather(Compound, `Vial Conc.`, -c(Sample, Type, BDE118)) %>%
    mutate(`Vial Conc. Corrected` = `Vial Conc.` * (BDE118.conc/BDE118)) %>%
    dplyr::select(-`Vial Conc.`)
  Recoveries <- IS.correction.BDE118 %>%
    dplyr::select(-BDE118) %>%
    spread(Compound, `Vial Conc. Corrected`) %>%
    mutate(Recovery = case_when(
      C %in% recoveries.pcb209 ~ (PCB209/PCB209.conc)*100,
      TRUE ~ 1000000)) %>%
    dplyr::select(Sample, Recovery)
  
  t.BDE118 <- IS.correction.BDE118 %>%
    filter(Compound %in% C) %>%
    dplyr::select(-BDE118) %>%
    filter(Type %in% "Calibration") %>%
    left_join(Recoveries) %>%
    mutate(LD = lod(model.ld.BDE118)$x * (50/500)/2) %>%
    mutate(LQ = loq(model.ld.BDE118)$x * (50/500)/2) %>%
    dplyr::select(Sample, Compound, LD, LQ) %>%
    left_join(quan.BDE118 %>% dplyr::select(Sample, Compound, Batch, `Injection Date`)) %>%
    filter(Sample %in% 33) %>%
    unique()
  t.BDE118
  
  M.parameters1[C, 3] <- t.BDE118$LD
  M.parameters1[C, 4] <- t.BDE118$LQ
  M.parameters1[C, 5] <- t.BDE118$Batch
  M.parameters1[C, 6] <- t.BDE118$`Injection Date`
  M.parameters1[C, 7] <- t.BDE118$Compound
}
M.parameters1

## Loop for compounds corrected by BDE209c in batch nº1 starts
for(c in 1:nC2){
  C <- as.character(unique(d.BDE209c$Compound))[c]
  d.now <- d.BDE209c %>%
    filter(Compound == C) %>%
    droplevels()
  
  # Non-corrected models
  model.PCB209 <- lm(Area ~ Calibration, data=d.BDE209c%>%
                       filter(Type %in%  "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% "PCB209"))
  model.BDE209c <- lm(Area ~ Calibration, data=d.BDE209c %>%
                        filter(Type %in%  "Calibration") %>%
                        filter(!is.na(Calibration)) %>%
                        filter(Compound %in% "BDE209c"))
  
  
  data.ld.BDE209c <- d.BDE209c %>%
    filter(Compound %in% c(C, "BDE209c")) %>%
    filter(Type %in% "Calibration") %>%
    filter(Sample %in% c("0.02", "0.04", "0.16", "0.4", "1.6", "3.33", "8.33", "16", "33")) %>%
    pivot_wider(names_from="Compound", values_from=c("Area", "Calibration")) %>%
    rename(Area = paste0("Area_", C),
           Calibration = paste0("Calibration_", C)) %>%
    mutate(y = Area / `Area_BDE209c`,
           x = Calibration / `Calibration_BDE209c`)
  
  #No corrected model by ISTD
  model.BDE209c <- lm(Area ~ Calibration, data=d.BDE209c %>%
                        filter(Type %in% "Calibration") %>%
                        filter(!is.na(Calibration)) %>%
                        filter(Compound %in% C) %>%
                        filter(Sample %in% c("0.02", "0.04", "0.16", "0.4", "1.6", "3.33", "8.33", "16", "33")))
  sum.BDE209c <- summary(model.BDE209c)
  print(sum.BDE209c)
  list(sum.BDE209c)
  M.parameters2[C, 1] <- sum.BDE209c$adj.r.squared
  
  quan.BDE209c <- d.BDE209c %>%
    filter(Compound %in% c(C, "BDE209c", "PCB209")) %>%
    mutate(`Injected Conc.`= case_when(
      Compound %in% C ~ (Area - model.BDE209c$coefficients[[1]]) / model.BDE209c$coefficients[[2]],
      Compound %in% "BDE209c" ~ (Area - model.BDE209c$coefficients[[1]]) / model.BDE209c$coefficients[[2]],
      Compound %in% "PCB209" ~ (Area - model.PCB209$coefficients[[1]]) / model.PCB209$coefficients[[2]],
      TRUE ~ 1000000)) %>%
    mutate(`Vial Conc.` = `Injected Conc.` / 2)
  
  #Corrected model by BDE118 or BDE209c
  model.ld.BDE209c <- lm(y~x, data = data.ld.BDE209c)
  sum.ld.BDE209c <- summary(model.ld.BDE209c)
  print(sum.ld.BDE209c)
  list(sum.ld.BDE209c)
  M.parameters2[C, 2] <- sum.ld.BDE209c$adj.r.squared
  # Now we have computed the adj R2 of the models  
  
  
  # Now we will compute LD/LQ
  IS.correction.BDE209c <- quan.BDE209c %>%
    dplyr::select(Sample, Type, Compound, `Vial Conc.`) %>%
    spread(Compound, `Vial Conc.`) %>%
    gather(Compound, `Vial Conc.`, -c(Sample, Type, BDE209c)) %>%
    mutate(`Vial Conc. Corrected` = `Vial Conc.` * (BDE209c.conc/BDE209c)) %>%
    dplyr::select(-`Vial Conc.`)
  Recoveries <- IS.correction.BDE209c %>%
    dplyr::select(-BDE209c) %>%
    spread(Compound, `Vial Conc. Corrected`) %>%
    mutate(Recovery = case_when(
      C %in% recoveries.pcb209 ~ (PCB209/PCB209.conc)*100,
      TRUE ~ 1000000)) %>%
    dplyr::select(Sample, Recovery)
  
  t.BDE209c <- IS.correction.BDE209c %>%
    filter(Compound %in% C) %>%
    dplyr::select(-BDE209c) %>%
    filter(Type %in% "Calibration") %>%
    left_join(Recoveries) %>%
    mutate(LD = lod(model.ld.BDE209c)$x * (50/500)/2) %>%
    mutate(LQ = loq(model.ld.BDE209c)$x * (50/500)/2) %>%
    dplyr::select(Sample, Compound, LD, LQ) %>%
    left_join(quan.BDE209c %>% dplyr::select(Sample, Compound, Batch, `Injection Date`)) %>%
    filter(Sample %in% 33) %>%
    unique()
  t.BDE209c
  
  M.parameters2[C, 3] <- t.BDE209c$LD
  M.parameters2[C, 4] <- t.BDE209c$LQ
  M.parameters2[C, 5] <- t.BDE209c$Batch
  M.parameters2[C, 6] <- t.BDE209c$`Injection Date`
  M.parameters2[C, 7] <- t.BDE209c$Compound
}
M.parameters2


# Merging M.parameters1 and M.parameters2 to M.parameters data base
M.parameters <- bind_rows(as.data.frame(M.parameters1),as.data.frame(M.parameters2)) %>%
  dplyr::select(Compound, R2, `R2-IS`, LD, LQ, Batch, `Injection Date`)
head(M.parameters)

DB <- bind_rows(DB, M.parameters) %>%
  dplyr::select(Compound, R2, `R2-IS`, LD, LQ, Batch, `Injection Date`)


######---------------------------------------------------------------------------------------------------------------

###----LDs/LQs
##CC1
#
#
batch <- 2   #250512_LD_LQ_PBDEs
d.BDE118 <- d.BDE118_T2
d.BDE209c <- d.BDE209c_T2
pcb209.inicial <- 50.822
bde118.inicial <- 10.03678088
bde209c.inicial <- 9.897961209

#----
## Calculate the concentrations of standards in samples

##PCB209
PCB209.conc <- (pcb209.inicial * 25)/50
##BDE118
BDE118.conc <- (bde118.inicial * 50)/50
##BDE209c
BDE209c.conc <- (bde209c.inicial * 50)/50
  
#Defining the surrogate (PCB209)
recoveries.pcb209 <- c("BDE28",   "BDE47",   "BDE100",  "BDE99",   "BDE154",  "BDE153", 
                         "BDE183",  "BDE209",  "BDE118",  "BDE209c", "PCB209")

#Defining a matrix for the compounds thath are corrected by BDE118
d.interesting1 <- d.BDE118 %>%
  filter(!Compound %in% c("PCB209", "BDE118", "BDE209c")) %>%
  dplyr::select(Compound) %>% unique()
nC <- length(unique(d.interesting1$Compound))
compound.labels1 <- as.character(unique(d.interesting1$Compound))

M.parameters1 <- array(NA, dim=c(nC, 7),
                       dimnames=list(Compound=compound.labels1,
                                     parameter=c("R2", "R2-IS",
                                                 "LD", "LQ",
                                                 "Batch", "Injection Date",
                                                 "Compound"))) 
M.parameters1

#Defining a matrix for the compounds thath are corrected by BDE209c
d.interesting2 <- d.BDE209c %>%
  filter(!Compound %in% c("PCB209", "BDE118", "BDE209c")) %>%
  dplyr::select(Compound) %>% unique()
nC2 <- length(unique(d.interesting2$Compound))
compound.labels2 <- as.character(unique(d.interesting2$Compound))

M.parameters2 <- array(NA, dim=c(nC2, 7),
                       dimnames=list(Compound=compound.labels2,
                                     parameter=c("R2", "R2-IS",
                                                 "LD", "LQ",
                                                 "Batch", "Injection Date",
                                                 "Compound"))) 
M.parameters2



## Loop for compounds corrected by BDE118 in batch nº1 starts
# c <- 1
for(c in 1:nC){
  C <- as.character(unique(d.BDE118$Compound))[c]
  d.now <- d.BDE118 %>%
    filter(Compound == C) %>%
    droplevels()
  
  # Non-corrected models
  model.PCB209 <- lm(Area ~ Calibration, data=d.BDE118 %>%
                       filter(Type %in%  "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% "PCB209"))
  model.BDE118 <- lm(Area ~ Calibration, data=d.BDE118 %>%
                       filter(Type %in%  "Calibration") %>%
                       filter(Compound %in% "BDE118"))
  
  data.ld.BDE118 <- d.BDE118 %>%
    filter(Compound %in% c(C, "BDE118")) %>%
    filter(Type %in% "Calibration") %>%
    filter(Sample %in% c("0.02", "0.04", "0.16", "0.4", "1.6", "3.33", "8.33", "16", "33")) %>%
    pivot_wider(names_from="Compound", values_from=c("Area", "Calibration")) %>%
    rename(Area = paste0("Area_", C),
           Calibration = paste0("Calibration_", C)) %>%
    mutate(y = Area / `Area_BDE118`,
           x = Calibration / `Calibration_BDE118`)
  
  
  #No corrected model by ISTD
  model.BDE118 <- lm(Area ~ Calibration, data=d.BDE118 %>%
                       filter(Type %in% "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% C) %>%
                       filter(Sample %in% c("0.02", "0.04", "0.16", "0.4", "1.6", "3.33", "8.33", "16", "33")))
  sum.BDE118 <- summary(model.BDE118)
  print(sum.BDE118)
  list(sum.BDE118)
  M.parameters1[C, 1] <- sum.BDE118$adj.r.squared
  
  
  quan.BDE118 <- d.BDE118 %>%
    filter(Compound %in% c(C, "BDE118", "PCB209")) %>%
    mutate(`Injected Conc.`= case_when(
      Compound %in% C ~ (Area - model.BDE118$coefficients[[1]]) / model.BDE118$coefficients[[2]],
      Compound %in% "BDE118" ~ (Area - model.BDE118$coefficients[[1]]) / model.BDE118$coefficients[[2]],
      Compound %in% "PCB209" ~ (Area - model.PCB209$coefficients[[1]]) / model.PCB209$coefficients[[2]],
      TRUE ~ 1000000)) %>%
    mutate(`Vial Conc.` = `Injected Conc.` / 2) 
  
  
  #Corrected model by BDE118 or BDE209c
  model.ld.BDE118 <- lm(y~x, data = data.ld.BDE118)
  sum.ld.BDE118 <- summary(model.ld.BDE118)
  print(sum.ld.BDE118)
  list(sum.ld.BDE118)
  M.parameters1[C, 2] <- sum.ld.BDE118$adj.r.squared
  # Now we have computed the adj R2 of the models  
  
  
  # Now we will compute LD/LQ
  IS.correction.BDE118 <- quan.BDE118 %>%
    dplyr::select(Sample, Type, Compound, `Vial Conc.`) %>%
    spread(Compound, `Vial Conc.`) %>%
    gather(Compound, `Vial Conc.`, -c(Sample, Type, BDE118)) %>%
    mutate(`Vial Conc. Corrected` = `Vial Conc.` * (BDE118.conc/BDE118)) %>%
    dplyr::select(-`Vial Conc.`)
  Recoveries <- IS.correction.BDE118 %>%
    dplyr::select(-BDE118) %>%
    spread(Compound, `Vial Conc. Corrected`) %>%
    mutate(Recovery = case_when(
      C %in% recoveries.pcb209 ~ (PCB209/PCB209.conc)*100,
      TRUE ~ 1000000)) %>%
    dplyr::select(Sample, Recovery)
  
  t.BDE118 <- IS.correction.BDE118 %>%
    filter(Compound %in% C) %>%
    dplyr::select(-BDE118) %>%
    filter(Type %in% "Calibration") %>%
    left_join(Recoveries) %>%
    mutate(LD = lod(model.ld.BDE118)$x * (50/500)/2) %>%
    mutate(LQ = loq(model.ld.BDE118)$x * (50/500)/2) %>%
    dplyr::select(Sample, Compound, LD, LQ) %>%
    left_join(quan.BDE118 %>% dplyr::select(Sample, Compound, Batch, `Injection Date`)) %>%
    filter(Sample %in% 33) %>%
    unique()
  t.BDE118
  
  M.parameters1[C, 3] <- t.BDE118$LD
  M.parameters1[C, 4] <- t.BDE118$LQ
  M.parameters1[C, 5] <- t.BDE118$Batch
  M.parameters1[C, 6] <- t.BDE118$`Injection Date`
  M.parameters1[C, 7] <- t.BDE118$Compound
}
M.parameters1

## Loop for compounds corrected by BDE209c in batch nº2 starts
for(c in 1:nC2){
  C <- as.character(unique(d.BDE209c$Compound))[c]
  d.now <- d.BDE209c %>%
    filter(Compound == C) %>%
    droplevels()
  
  # Non-corrected models
  model.PCB209 <- lm(Area ~ Calibration, data=d.BDE209c%>%
                       filter(Type %in%  "Calibration") %>%
                       filter(!is.na(Calibration)) %>%
                       filter(Compound %in% "PCB209"))
  model.BDE209c <- lm(Area ~ Calibration, data=d.BDE209c %>%
                        filter(Type %in%  "Calibration") %>%
                        filter(!is.na(Calibration)) %>%
                        filter(Compound %in% "BDE209c"))
  
  
  data.ld.BDE209c <- d.BDE209c %>%
    filter(Compound %in% c(C, "BDE209c")) %>%
    filter(Type %in% "Calibration") %>%
    filter(Sample %in% c("0.02", "0.04", "0.16", "0.4", "1.6", "3.33", "8.33", "16", "33")) %>%
    pivot_wider(names_from="Compound", values_from=c("Area", "Calibration")) %>%
    rename(Area = paste0("Area_", C),
           Calibration = paste0("Calibration_", C)) %>%
    mutate(y = Area / `Area_BDE209c`,
           x = Calibration / `Calibration_BDE209c`)
  
  #No corrected model by ISTD
  model.BDE209c <- lm(Area ~ Calibration, data=d.BDE209c %>%
                        filter(Type %in% "Calibration") %>%
                        filter(!is.na(Calibration)) %>%
                        filter(Compound %in% C) %>%
                        filter(Sample %in% c("0.02", "0.04", "0.16", "0.4", "1.6", "3.33", "8.33", "16", "33")))
  sum.BDE209c <- summary(model.BDE209c)
  print(sum.BDE209c)
  list(sum.BDE209c)
  M.parameters2[C, 1] <- sum.BDE209c$adj.r.squared
  
  quan.BDE209c <- d.BDE209c %>%
    filter(Compound %in% c(C, "BDE209c", "PCB209")) %>%
    mutate(`Injected Conc.`= case_when(
      Compound %in% C ~ (Area - model.BDE209c$coefficients[[1]]) / model.BDE209c$coefficients[[2]],
      Compound %in% "BDE209c" ~ (Area - model.BDE209c$coefficients[[1]]) / model.BDE209c$coefficients[[2]],
      Compound %in% "PCB209" ~ (Area - model.PCB209$coefficients[[1]]) / model.PCB209$coefficients[[2]],
      TRUE ~ 1000000)) %>%
    mutate(`Vial Conc.` = `Injected Conc.` / 2)
  
  #Corrected model by BDE118 or BDE209c
  model.ld.BDE209c <- lm(y~x, data = data.ld.BDE209c)
  sum.ld.BDE209c <- summary(model.ld.BDE209c)
  print(sum.ld.BDE209c)
  list(sum.ld.BDE209c)
  M.parameters2[C, 2] <- sum.ld.BDE209c$adj.r.squared
  # Now we have computed the adj R2 of the models  
  
  
  # Now we will compute LD/LQ
  IS.correction.BDE209c <- quan.BDE209c %>%
    dplyr::select(Sample, Type, Compound, `Vial Conc.`) %>%
    spread(Compound, `Vial Conc.`) %>%
    gather(Compound, `Vial Conc.`, -c(Sample, Type, BDE209c)) %>%
    mutate(`Vial Conc. Corrected` = `Vial Conc.` * (BDE209c.conc/BDE209c)) %>%
    dplyr::select(-`Vial Conc.`)
  Recoveries <- IS.correction.BDE209c %>%
    dplyr::select(-BDE209c) %>%
    spread(Compound, `Vial Conc. Corrected`) %>%
    mutate(Recovery = case_when(
      C %in% recoveries.pcb209 ~ (PCB209/PCB209.conc)*100,
      TRUE ~ 1000000)) %>%
    dplyr::select(Sample, Recovery)
  
  t.BDE209c <- IS.correction.BDE209c %>%
    filter(Compound %in% C) %>%
    dplyr::select(-BDE209c) %>%
    filter(Type %in% "Calibration") %>%
    left_join(Recoveries) %>%
    mutate(LD = lod(model.ld.BDE209c)$x * (50/500)/2) %>%
    mutate(LQ = loq(model.ld.BDE209c)$x * (50/500)/2) %>%
    dplyr::select(Sample, Compound, LD, LQ) %>%
    left_join(quan.BDE209c %>% dplyr::select(Sample, Compound, Batch, `Injection Date`)) %>%
    filter(Sample %in% 33) %>%
    unique()
  t.BDE209c
  
  M.parameters2[C, 3] <- t.BDE209c$LD
  M.parameters2[C, 4] <- t.BDE209c$LQ
  M.parameters2[C, 5] <- t.BDE209c$Batch
  M.parameters2[C, 6] <- t.BDE209c$`Injection Date`
  M.parameters2[C, 7] <- t.BDE209c$Compound
}
M.parameters2



# Merging M.parameters1 and M.parameters2 to M.parameters data base
M.parameters <- bind_rows(as.data.frame(M.parameters1),as.data.frame(M.parameters2)) %>%
  dplyr::select(Compound, R2, `R2-IS`, LD, LQ, Batch, `Injection Date`)
head(M.parameters)

DB <- bind_rows(DB, M.parameters) %>%
  dplyr::select(Compound, R2, `R2-IS`, LD, LQ, Batch, `Injection Date`)

head(DB)


DB_tot <- DB %>%
  mutate(Compound = as.factor(Compound),
         R2 = as.numeric(R2),
         `R2-IS` = as.numeric(`R2-IS`),
         LD = as.numeric(LD),
         LQ = as.numeric(LQ))
rownames(DB_tot) <- NULL
write.xlsx(DB_tot, file = "Table_1_2.xlsx")






