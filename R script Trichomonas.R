### R script written by Manuela Merling de Chapa
setwd("C:/Users/manum/Desktop/Manuskripte/Trichos/submission/Statistik")

## how condition was calculated ############################################
## Condition 
library(readxl)

### goshawk ###############################################################
Condition <- read_excel("./prevalence_model.xlsx")

Condition2 <- Condition[, c("Ring", "Weight", "Wing", "Sex", "Species")]
Condition2 <- droplevels(na.omit(Condition2))
Condition2 <- as.data.frame(Condition2)

goshawk_urban <- subset(Condition2, Species == "Goshawk")
gosawk_urban <- droplevels(na.omit(goshawk_urban))

female <- subset(goshawk_urban, Sex == "female")
male <- subset(goshawk_urban, Sex == "male")

plot(Weight~Wing, data = female)
y <- lm(Weight~Wing, data = female)
abline(y)
title(main = "Condition female goshawks", cex.main = 0.9)

cond<-rstandard(y)

female$Condition2<-cond

female

plot(Weight~Wing, data = male)
x <- lm(Weight~Wing, data = male)
abline(x)
title(main = "Condition male goshawks", cex.main = 0.9)

cond2<-rstandard(x)

male$Condition2 <- cond2

sex <- rbind(female,male)
sex <- sex[, c("Ring", "Condition2")]

#install.packages("data.table")
library(data.table)

setDT(Condition)

setDT(sex)

Condition<-sex[Condition, on=c('Ring')]

## Condition
### sparrowhawk #############################################################

Condition <- read_excel("./prevalence_model.xlsx")

Condition2 <- Condition[, c("Ring", "Weight", "Wing", "Sex", "Species")]
Condition2 <- droplevels(na.omit(Condition2))
Condition2 <- as.data.frame(Condition2)

sparrowhawk <- subset(Condition2, Species == "Sparrowhawk")
sparrowhawk <- droplevels(na.omit(sparrowhawk))

female <- subset(sparrowhawk, Sex == "female")
male <- subset(sparrowhawk, Sex == "male")

plot(Weight~Wing, data = female)
y <- lm(Weight~Wing, data = female)
abline(y)
title(main = "Condition female sparrowhawks", cex.main = 0.9)

cond<-rstandard(y)

female$Condition2<-cond

female

plot(Weight~Wing, data = male)
x <- lm(Weight~Wing, data = male)
abline(x)
title(main = "Condition male sparrowhawks", cex.main = 0.9)

cond2<-rstandard(x)

male$Condition2 <- cond2

sex <- rbind(female,male)
sex <- sex[, c("Ring", "Condition2")]

#install.packages("data.table")
library(data.table)

setDT(Condition)

setDT(sex)

Condition<-sex[Condition, on=c('Ring')]


## GLMMs #################################################################

## prevalence ###########################################################
library(readxl)
library(spaMM)
library(DHARMa)
library(lme4) 
library(car)
library(purrr)
library(binom)
library(ggplot2)
library(maps)
library(mapdata)
library(dplyr)
library(doSNOW)
spaMM.options(nb_cores = 3, separation_max = 1)
boot.repl <- 1000
nb_cores <- 3

strain_master <- read_excel("./prevalence_model.xlsx")

strain2 <- strain_master[, c("Year", "Location","No_nestlings", "Sex",
                             "Territory", "Habitat", "Species", "Species_Habitat",
                             "Age", "Condition2", "Prevalence", "Age_cat")]

strain2$Year <- as.factor(strain2$Year)
strain2$Location <- as.factor(strain2$Location)
strain2$Territory <- as.factor(strain2$Territory)
strain2$Sex <- as.factor(strain2$Sex)
strain2$Habitat <- as.factor(strain2$Habitat)
strain2$Species <- as.factor(strain2$Species)
strain2$Age_cat <- as.factor(strain2$Age_cat)
strain2$Species_Habitat <- as.factor(strain2$Species_Habitat)
strain2 <- as.data.frame(strain2)
strain2 <- droplevels(na.omit(strain2)) 
str(strain2)
head(strain2)
range(strain2$Condition2)

glmm_prevalence <- fitme(Prevalence ~ Species_Habitat + Condition2 + Age_cat + Sex + Year
                       + (1|Location/Territory),
                       family = binomial(link = "logit"), 
                       data = strain2, method = "PQL/L")

### odds
#exp(glmm_prevalence$fixef["SpeciesSparrowhawk"])
exp(glmm_prevalence$fixef)

#### test effects
glmm_prevalence_noSpecies <- fitme(Prevalence ~  Year + Sex + Age_cat + Condition2 + 
                                   (1|Location/Territory),
                                 family = binomial(link = "logit"), data = strain2, 
                                 method = "PQL/L")

glmm_prevalence_noYear <- fitme(Prevalence ~  Species_Habitat + Sex + Age_cat + Condition2 + 
                                (1|Location/Territory),
                              family = binomial(link = "logit"), data = strain2, 
                              method = "PQL/L")

glmm_prevalence_noSex <- fitme(Prevalence ~  Species_Habitat + Year + Age_cat + Condition2 + 
                               (1|Location/Territory),
                             family = binomial(link = "logit"), data = strain2, 
                             method = "PQL/L")

glmm_prevalence_noAge <- fitme(Prevalence ~  Species_Habitat + Sex + Year + Condition2 + 
                               (1|Location/Territory),
                             family = binomial(link = "logit"), data = strain2, 
                             method = "PQL/L")

glmm_prevalence_noCondition <- fitme(Prevalence ~  Species_Habitat + Year + Sex + Age_cat + 
                                     (1|Location/Territory),
                                   family = binomial(link = "logit"), data = strain2, 
                                   method = "PQL/L")

anova(glmm_prevalence, glmm_prevalence_noSpecies, boot.repl = boot.repl, nb_cores =nb_cores)
anova(glmm_prevalence, glmm_prevalence_noYear, boot.repl = boot.repl, nb_cores =nb_cores)
anova(glmm_prevalence, glmm_prevalence_noSex, boot.repl = boot.repl, nb_cores =nb_cores)
anova(glmm_prevalence, glmm_prevalence_noAge, boot.repl = boot.repl, nb_cores =nb_cores)
anova(glmm_prevalence, glmm_prevalence_noCondition, boot.repl = boot.repl, nb_cores =nb_cores)


#### glmer Species effect
GLMM_infect2_bis <- glmer(Prevalence~Species_Habitat+Age_cat+Sex+Condition2+Year+(1|Location/Territory), family=binomial(link="logit"), 
                          nAGQ = 0, data=strain2, control = glmerControl(optimizer = c("bobyqa")))

logLik(glmm_prevalence)
logLik(GLMM_infect2_bis)

contr <- c("Goshawk - Sparrowhawk" = c(1/2, 1/2, -1))
library(multcomp)
summary(glht(GLMM_infect2_bis))
summary(glht(GLMM_infect2_bis, linfct = mcp(Species_Habitat = contr)))
CI <- confint(glht(GLMM_infect2_bis, linfct = mcp(Species_Habitat = contr))) ##  on the logit scale
exp(CI$confint) ## On the odd ratio scale, exp calculate odd ratio (chance of beeing infected urban compared rural), plus CI


#### plots prevalence
### year betwen species
plot_data <- strain2 %>%
  group_by(Year, Species) %>%
  summarise(n = n(), 
            n_inf = sum(Prevalence == 1, na.rm = T), 
            n_not_inf = sum(Prevalence == 0, na.rm = T), 
            prop_inf = n_inf / n) %>% 
  ungroup() %>% 
  mutate(CI_lwr = map2_dbl(n_inf, n, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 5])), ## method exact
         CI_upr = map2_dbl(n_inf, n, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 6]))) ## method exact

plot_data$Species <- relevel(as.factor(plot_data$Species), ref = "Goshawk")

ggplot(plot_data, aes(y = prop_inf, x = Year, ymin = CI_lwr, ymax = CI_upr, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.8) +
  geom_errorbar(width = 0, position =  position_dodge(width =  0.9)) +
  scale_fill_grey() +
  scale_y_continuous("Proportion of infected nestlings") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 

ggsave(filename = "./year_plot.pdf", width = 8, height = 8)

### age cat between species
plot_data <- strain2 %>%
  group_by(Age_cat, Species) %>%
  summarise(n = n(), 
            n_inf = sum(Prevalence == 1, na.rm = T), 
            n_not_inf = sum(Prevalence == 0, na.rm = T), 
            prop_inf = n_inf / n) %>% 
  ungroup() %>% 
  mutate(CI_lwr = map2_dbl(n_inf, n, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 5])), ## method exact
         CI_upr = map2_dbl(n_inf, n, ~ as.numeric(binom.confint( x = .x, n = .y)[5, 6]))) ## method exact

plot_data$Species <- relevel(as.factor(plot_data$Species), ref = "Goshawk")
plot_data$Age_cat <- relevel(as.factor(plot_data$Age_cat), ref = "young")

ggplot(plot_data, aes(y = prop_inf, x = Age_cat, ymin = CI_lwr, ymax = CI_upr, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.8) +
  geom_errorbar(width = 0, position =  position_dodge(width =  0.9)) +
  scale_fill_grey() +
  scale_y_continuous("Proportion of infected nestlings") +
  scale_x_discrete("Age category") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 

ggsave(filename = "./age_species_plot.pdf", width = 8, height = 8)


## clinical signs #########################################################
library(readxl)
library(spaMM)
library(DHARMa)
library(lme4) 
library(car)
library(purrr)
library(binom)
library(ggplot2)
library(maps)
library(mapdata)
library(dplyr)
library(doSNOW)

### fisher test 
### goshawks
clinical_goshawk <- read_excel("./Symptoms_goshawk.xlsx")
clinical_goshawk$Strain <- as.factor(clinical_goshawk$Strain)
str(clinical_goshawk)

clinical_test2 <- clinical_goshawk[, c("yes", "no")]

fisher.test(clinical_test2, simulate.p.value = TRUE)

### sparrowhawk
clinical_ACNI <- read_excel("./symptoms_ACNI.xlsx")
clinical_ACNI$Strain <- as.factor(clinical_ACNI$Strains)
str(clinical_ACNI)

clinical_test <- clinical_ACNI[, c("yes", "no")]

fisher.test(clinical_test, simulate.p.value = TRUE)

### double 
x <- matrix(c(8,18,17,179),2,2)                          
dimnames(x) <-  list(c("Clinical signs", "no"), c("double yes", "no"))
x
fisher.test(x)

### GLMM  ##############
strain_master <- read_excel("./clinical_model.xlsx")

strain3 <- strain_master[, c("Year", "Location","No_nestlings", "Sex", "Condition2",
                             "Territory", "Habitat", "Species", "Species_Habitat", "Age_cat",
                             "Clinical_binary", "Strain_cat4","Ring", "double_infection")]

strain3$Year <- as.factor(strain3$Year)
strain3$Location <- as.factor(strain3$Location)
strain3$Territory <- as.factor(strain3$Territory)
strain3$Sex <- as.factor(strain3$Sex)
strain3$Habitat <- as.factor(strain3$Habitat)
strain3$Species <- as.factor(strain3$Species)
strain3$Age_cat <- as.factor(strain3$Age_cat)
strain3$Strain_cat4 <- as.factor(strain3$Strain_cat4)
strain3$Species_Habitat <- as.factor(strain3$Species_Habitat)
strain3 <- as.data.frame(strain3)
strain3 <- droplevels(na.omit(strain3))
str(strain3)
head(strain3)
table(strain3$Strain_cat4)
head(strain3$Strain_cat4)
table(strain3$Species)

glmm_clinical <- fitme(Clinical_binary ~ Species + Strain_cat4
                       + Year + Sex + (1|Territory),
                       family = binomial(link = "logit"), 
                       data = strain3, method = "PQL/L")

### test effects
glmm_clinical_noSpecies <- fitme(Clinical_binary ~ Strain_cat4 + Year
                                 + Sex + (1|Territory),
                                 family = binomial(link = "logit"), data = strain3, 
                                 method = "PQL/L")

glmm_clinical_noYear <- fitme(Clinical_binary ~ Species + Strain_cat4 + Sex                               + (1|Territory),
                              family = binomial(link = "logit"), data = strain3, 
                              method = "PQL/L")

glmm_clinical_nostrain <- fitme(Clinical_binary ~ Species + Year + Sex
                                + (1|Territory),
                                family = binomial(link = "logit"), data = strain3, 
                                method = "PQL/L")

glmm_clinical_nosex <- fitme(Clinical_binary ~ Species + Year + Strain_cat4
                             + (1|Territory),
                             family = binomial(link = "logit"), data = strain3, 
                             method = "PQL/L")

anova(glmm_clinical, glmm_clinical_noSpecies, boot.repl = boot.repl, nb_cores =nb_cores)
anova(glmm_clinical, glmm_clinical_noYear, boot.repl = boot.repl, nb_cores =nb_cores)
anova(glmm_clinical, glmm_clinical_nostrain, boot.repl = boot.repl, nb_cores =nb_cores) 
anova(glmm_clinical, glmm_clinical_nosex, boot.repl = boot.repl, nb_cores =nb_cores)

#### figure
pred_goshawk_A1 <- predict(glmm_clinical,
                           newdata = expand.grid(Species = "Goshawk",
                                                 Year = "2015", Sex ="female",
                                                 Strain_cat4 = "A1"),
                           re.form = NA, intervals = "predVar") 
pred_1 <- data.frame(pred = pred_goshawk_A1, attr(pred_goshawk_A1, "intervals"), 
                     Species = "Goshawk", Strain_cat4 = "A1")

pred_sparrowhawk_A1 <- predict(glmm_clinical,
                               newdata = expand.grid(Species = "Sparrowhawk",
                                                     Year = "2015", Sex ="female",
                                                     Strain_cat4 = "A1"),
                               re.form = NA, intervals = "predVar")
pred_2 <- data.frame(pred = pred_sparrowhawk_A1, attr(pred_sparrowhawk_A1, "intervals"), 
                     Species = "Sparrowhawk", Strain_cat4 = "A1")


pred_goshawk_C4 <- predict(glmm_clinical,
                           newdata = expand.grid(Species = "Goshawk",
                                                 Sex ="female",
                                                 Year = "2015",
                                                 Strain_cat4 = "C4"),
                           re.form = NA, intervals = "predVar")
pred_3 <- data.frame(pred = pred_goshawk_C4, attr(pred_goshawk_C4, "intervals"), 
                     Species = "Goshawk", Strain_cat4 = "C4")

pred_sparrowhawk_C4 <- predict(glmm_clinical,
                               newdata = expand.grid(Species = "Sparrowhawk",
                                                     Year = "2015", Sex ="female",
                                                     Strain_cat4 = "C4"),
                               re.form = NA, intervals = "predVar")
pred_4 <- data.frame(pred = pred_sparrowhawk_C4, attr(pred_sparrowhawk_C4, "intervals"), 
                     Species = "Sparrowhawk", Strain_cat4 = "C4")


pred_goshawk_anders <- predict(glmm_clinical,
                               newdata = expand.grid(Species = "Goshawk",
                                                     Year = "2015", Sex ="female",
                                                     Strain_cat4 = "different"),
                               re.form = NA, intervals = "predVar") 
pred_5 <- data.frame(pred = pred_goshawk_anders, attr(pred_goshawk_anders, "intervals"), 
                     Species = "Goshawk", Strain_cat4 = "different")

pred_sparrowhawk_anders <- predict(glmm_clinical,
                                   newdata = expand.grid(Species = "Sparrowhawk",
                                                         Year = "2015", Sex ="female",
                                                         Strain_cat4 = "different"),
                                   re.form = NA, intervals = "predVar")
pred_6 <- data.frame(pred = pred_sparrowhawk_anders, attr(pred_sparrowhawk_anders, "intervals"), 
                     Species = "Sparrowhawk", Strain_cat4 = "different")



pred_goshawk_double <- predict(glmm_clinical,
                               newdata = expand.grid(Species = "Goshawk",
                                                     Sex ="female",
                                                     Year = "2015",
                                                     Strain_cat4 = "double"),
                               re.form = NA, intervals = "predVar") 
pred_7 <- data.frame(pred = pred_goshawk_double, attr(pred_goshawk_double, "intervals"), 
                     Species = "Goshawk", Strain_cat4 = "double")

pred_sparrowhawk_double <- predict(glmm_clinical,
                                   newdata = expand.grid(Species = "Sparrowhawk",
                                                         Year = "2015", Sex ="female",
                                                         Strain_cat4 = "double"),
                                   re.form = NA, intervals = "predVar")
pred_8 <- data.frame(pred = pred_sparrowhawk_double, attr(pred_sparrowhawk_double, "intervals"), 
                     Species = "Sparrowhawk", Strain_cat4 = "double")


pred_Species_strain <- rbind(pred_1, pred_2, pred_3, pred_4, pred_5, pred_6, pred_7, pred_8)
pred_Species_strain$Species <- factor(pred_Species_strain$Species, levels = c("Goshawk", "Sparrowhawk"))
pred_Species_strain$Strain_cat4 <- factor(pred_Species_strain$Strain_cat4, levels = c("A1", "C4", "different", "double"))

ggplot(pred_Species_strain, aes(x = Strain_cat4, y = pred, fill = Species, ymin = predVar_0.025, ymax = predVar_0.975)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.8) +
  geom_errorbar(width = 0, position = position_dodge(preserve = "total", width = 0.9)) +
  scale_fill_grey() +
  scale_y_continuous("Proportion of clinical signs") +
  scale_x_discrete("Strain category") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank()) 

ggsave(filename = "./strain_clinical.pdf", width = 8, height = 8)
