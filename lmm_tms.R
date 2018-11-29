# lmm_tms.R
# https://osf.io/r73xg/

# Modified for OSF: 10/30/2018 By: Tyler Santander, Allison Shapiro, Evan Layher

## --- LICENSE INFORMATION --- ##
## Modified BSD-2 License - for Non-Commercial Use Only

## Copyright (c) 2017-18, The Regents of the University of California
## All rights reserved.

## Redistribution and use in source and binary forms, with or without modification, are 
## permitted for non-commercial use only provided that the following conditions are met:

## 1. Redistributions of source code must retain the above copyright notice, this list 
##    of conditions and the following disclaimer.

## 2. Redistributions in binary form must reproduce the above copyright notice, this list 
##    of conditions and the following disclaimer in the documentation and/or other 
##    materials provided with the distribution.

## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
## EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
## OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
## SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
## INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
## TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
## OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY 
## WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## For permission to use for commercial purposes, please contact UCSB’s Office of 
## Technology & Industry Alliances at 805-893-5180 or info@tia.ucsb.edu.
## --------------------------- ##

# File paths
statsDir <- ''
codeDir <- ''
outDir <- ''

# Stat files
dataFile1 <- paste(statsDir, 'tms1_sdt.csv', sep = '', collapse = NULL)
dataFile2 <- paste(statsDir, 'tms2_sdt.csv', sep = '', collapse = NULL)

# Code files
source(paste(codeDir, 'summarizeWithin.R', sep = '', collapse = NULL))

# Output files
fig5 <- paste(outDir, 'fig5.tiff', sep = '', collapse = NULL)
fig6 <- paste(outDir, 'fig6.tiff', sep = '', collapse = NULL)
fig7 <- paste(outDir, 'fig7.tiff', sep = '', collapse = NULL)
fig8 <- paste(outDir, 'fig8.tiff', sep = '', collapse = NULL)
tab2 <- paste(outDir, 'tab1.tiff', sep = '', collapse = NULL)
tab3 <- paste(outDir, 'tab2.tiff', sep = '', collapse = NULL)

# Load packages
require(broom)
require(car)
require(dplyr)
require(effects)
require(forcats)
require(ggplot2)
require(gridExtra)
require(lme4)
require(merTools)
require(MuMIn)
require(predictmeans)
require(simr)
require(sjmisc)
require(sjPlot)
require(tidyr)
require(stargazer)

# Import data files
data1 <- read.csv(dataFile1)
data2 <- read.csv(dataFile2)

# Rename factors
data1 <- data1 %>%
  mutate(sub = factor(sub)) %>%
  select(sub = sub,
         cond = conVlib,
         time = preVstm,
         loc = stim,
         c = normC,
         d = d)

data2 <- data2 %>%
  mutate(sub = factor(sub)) %>%
  select(sub = sub,
         cond = conVlib,
         time = preVstm,
         loc = stim,
         c = normC,
         d = d)

# Obtain normalized c and d' means
cMean1 <- summarySEwithin(data1, measurevar = 'c', withinvars = c("cond", "time"), idvar = 'sub')
cMean2 <- summarySEwithin(data2, measurevar = 'c', withinvars = c("cond", "time"), idvar = 'sub')
dMean1 <- summarySEwithin(data1, measurevar = 'd', withinvars = c("time"), idvar = 'sub')
dMean2 <- summarySEwithin(data2, measurevar = 'd', withinvars = c("time"), idvar = 'sub')

# Choose references
data1$cond <- relevel(data1$cond, ref = 'Con') # Reference (Con/Lib)
data1$time <- relevel(data1$time, ref = 'Stm') # Reference (Pre/Stm)
data2$cond <- relevel(data2$cond, ref = 'Con') # Reference (Con/Lib)
data2$time <- relevel(data2$time, ref = 'Stm') # Reference (Pre/Stm)

# Change contrasts from default of 'treatment coding' to 'deviation coding'
contrasts(data1$time) <- contr.sum(2)
contrasts(data1$loc) <- contr.sum(3)
contrasts(data2$time) <- contr.sum(2)
contrasts(data2$loc) <- contr.sum(3)

colnames(contrasts(data1$loc)) <- c('Rifg', 'Rmfg')
colnames(contrasts(data2$loc)) <- c('Rifg', 'Rdlpfc')

# 3 way interaction linear mixed effects (sub as random effect)
lmm1 <- lmer(c ~ cond * time * loc + (1|sub), data = data1) # Exp. 1 normalized c
lmm2 <- lmer(c ~ cond * time * loc + (1|sub), data = data2) # Exp. 2 normalized c

lmm1d <- lmer(d ~ cond * time * loc + (1|sub), data = data1) # Exp. 1 d'
lmm2d <- lmer(d ~ cond * time * loc + (1|sub), data = data2) # Exp. 2 d'

# LSD correction
lsdEst1 <- predictmeans(model = lmm1,
                        modelterm = 'cond:time:loc',
                        pairwise = TRUE,
                        adj = 'tukey',
                        mplot = FALSE, 
                        pplot = FALSE, 
                        bkplot = FALSE)
lsdFit1 <- as.data.frame(lsdEst1[[1]])
lsd1    <- as.numeric(lsdEst1[[4]][3])
lsdFit1 <- lsdFit1 %>%
  mutate(lwr = Freq - lsd1,
         upr = Freq + lsd1)

lsdEst2 <- predictmeans(model = lmm2,
                       modelterm = 'cond:time:loc',
                       pairwise = TRUE,
                       adj = 'tukey',
                       mplot = FALSE, 
                       pplot = FALSE, 
                       bkplot = FALSE)
lsdFit2 <- as.data.frame(lsdEst2[[1]])
lsd2    <- as.numeric(lsdEst2[[4]][3])
lsdFit2 <- lsdFit2 %>%
  mutate(lwr = Freq - lsd2,
         upr = Freq + lsd2)

lsdEst1d <- predictmeans(model = lmm1d,
                        modelterm = 'cond:time:loc',
                        pairwise = TRUE,
                        adj = 'tukey',
                        mplot = FALSE, 
                        pplot = FALSE, 
                        bkplot = FALSE)
lsdFit1d <- as.data.frame(lsdEst1d[[1]])
lsd1d    <- as.numeric(lsdEst1d[[4]][3])
lsdFit1d <- lsdFit1d %>%
  mutate(lwr = Freq - lsd1d,
         upr = Freq + lsd1d)

lsdEst2d <- predictmeans(model = lmm2d,
                        modelterm = 'cond:time:loc',
                        pairwise = TRUE,
                        adj = 'tukey',
                        mplot = FALSE, 
                        pplot = FALSE, 
                        bkplot = FALSE)
lsdFit2d <- as.data.frame(lsdEst2d[[1]])
lsd2d    <- as.numeric(lsdEst2d[[4]][3])
lsdFit2d <- lsdFit2d %>%
  mutate(lwr = Freq - lsd2d,
         upr = Freq + lsd2d)

# FIGURE 5: Normalized c, experiments 1 & 2
tiff(fig5, height = 6, width = 12, units = 'in', res = 300)
groupdat <- merge(data1, lsdFit1, all = TRUE)

# Rename factor levels
groupdat <- groupdat %>% mutate(time = fct_recode(time, "pre" = "Pre", "post" = "Stm"),
                            loc = fct_recode(loc, "rIFG" = "Rifg", "rMFG" = "Rmfg", "SHAM" = "Sham"),
                            cond = fct_recode(cond, "Liberal" = "Lib", "Conservative" = "Con"))

# Reorder location factor
groupdat$loc <- factor(groupdat$loc, levels = c("SHAM", "rIFG", "rMFG"))

# Plot data
ex1plot <- ggplot(groupdat, aes(x = forcats::fct_rev(time), y = c, group = sub, colour = sub)) +
  geom_line() + scale_colour_grey(start = .2, end = .8, guide = FALSE) +
  geom_line(aes(y = Freq), colour = "firebrick", size = 1, alpha = .2) +
  scale_y_continuous(limits = c(-2.1, 2.1)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = .2, colour = "firebrick", size = 1, alpha = .2) +
  facet_grid(cond~loc) +
  theme_minimal() +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(face = "bold", size = 14),
        legend.title = element_text(size = 10),
        strip.text.x = element_text(size = 20, colour = "black", face = "bold"),
        strip.text.y = element_text(size = 0, colour = "black", face = "bold")) +
  ylab('Normalized c') + xlab('Test Time') +
  ggtitle("Experiment 1") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = '0.5', margin = margin(0,0,10,0)))

exp2dat <- merge(data2, lsdFit2, all = TRUE)

# Rename factor levels
exp2dat <- exp2dat %>% mutate(time = fct_recode(time, "pre" = "Pre", "post" = "Stm"),
                            loc = fct_recode(loc, "rIFG" = "Rifg", "rDLPFC" = "Rdlpfc", "SHAM" = "Sham"),
                            cond = fct_recode(cond, "Liberal" = "Lib", "Conservative" = "Con"))

# Reorder location factor
exp2dat$loc <- factor(exp2dat$loc, levels = c("SHAM", "rIFG", "rDLPFC"))

# Plot data
ex2plot <- ggplot(exp2dat, aes(x = forcats::fct_rev(time), y = c, group = sub, colour = sub)) +
  geom_line() + scale_colour_grey(start = .2, end = .8, guide = FALSE) +
  geom_line(aes(y = Freq), colour = "firebrick", size = 1, alpha = .2) +
  scale_y_continuous(limits = c(-2.1,2.1)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), colour = "firebrick", width = .2, size =1 , alpha = .2) +
  facet_grid(cond~loc) +
  theme_minimal() +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(face = "bold", size = 14),
        legend.title = element_text(size = 10),
        strip.text.x = element_text(size = 20, colour = "black", face = "bold"),
        strip.text.y = element_text(size = 20, colour = "black", face = "bold")) +
  ylab('') + 
  xlab('Test Time') +
  ggtitle("Experiment 2") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = '0.5', margin = margin(0,0,10,0)))
  
grid.arrange(ex1plot, ex2plot, nrow = 1)
dev.off()

# FIGURE 6: d', experiments 1 & 2
tiff(fig6, height = 6, width = 12, units = 'in', res=300)
groupdat <- merge(data1, lsdFit1d, all = TRUE)

# Rename factor levels
groupdat <- groupdat %>% mutate(time = fct_recode(time, "pre" = "Pre", "post" = "Stm"),
                            loc =fct_recode(loc, "rIFG" = "Rifg", "rMFG" = "Rmfg", "SHAM" = "Sham"),
                            cond=fct_recode(cond, "Liberal" = "Lib", "Conservative" = "Con"))

# Reorder location factor
groupdat$loc <- factor(groupdat$loc, levels = c("SHAM", "rIFG", "rMFG"))

# Plot data
ex1plot <- ggplot(groupdat, aes(x = forcats::fct_rev(time), y = d, group = sub, colour = sub)) +
  geom_line() + scale_colour_grey(start = .2, end = .8, guide = FALSE) +
  geom_line(aes(y = Freq), colour = "firebrick", size = 1, alpha = .2) +
  scale_y_continuous(limits = c(-0.6,1.6)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = .2, colour = "firebrick", size = 1, alpha = .2) +
  facet_grid(cond~loc) +
  theme_minimal() +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(face = "bold", size = 14),
        legend.title = element_text(size = 10),
        strip.text.x = element_text(size = 20, colour = "black", face = "bold"),
        strip.text.y = element_text(size = 0, colour = "black", face = "bold")) +
  ylab("d'") + xlab('Test Time') +
  ggtitle("Experiment 1") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = '0.5', margin = margin(0,0,10,0)))

exp2dat<-merge(data2, lsdFit2d, all = TRUE)

# Rename factor levels
exp2dat <- exp2dat %>% mutate(time = fct_recode(time, "pre" = "Pre", "post" = "Stm"),
                          loc = fct_recode(loc, "rIFG" = "Rifg", "rDLPFC" = "Rdlpfc", "SHAM" = "Sham"),
                          cond = fct_recode(cond, "Liberal" = "Lib", "Conservative" = "Con"))

# Reorder location factor
exp2dat$loc <- factor(exp2dat$loc, levels = c("SHAM", "rIFG", "rDLPFC"))

# Plot data
ex2plot <- ggplot(exp2dat, aes(x = forcats::fct_rev(time), y = d, group = sub, colour = sub)) +
  geom_line() + scale_colour_grey(start = .2, end = .8, guide = FALSE) +
  geom_line(aes(y = Freq), colour = "firebrick", size = 1, alpha = .2) +
  scale_y_continuous(limits = c(-1,1.6)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), colour = "firebrick", width = .2, size = 1, alpha = .2) +
  facet_grid(cond~loc) +
  theme_minimal() +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(face = "bold", size = 14),
        legend.title = element_text(size = 10),
        strip.text.x = element_text(size = 20, colour = "black", face = "bold"),
        strip.text.y = element_text(size = 20, colour = "black", face = "bold")) +
  ylab('') + 
  xlab('Test Time') +
  ggtitle("Experiment 2") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = '0.5', margin = margin(0,0,10,0)))

grid.arrange(ex1plot, ex2plot, nrow=1)
dev.off()

# FIGURE 7:

# Obtain posterior mean & 95% CIs
fe1C <- FEsim(lmm1, 1000)
fe1C <- fe1C %>%
  mutate(lwr = mean - 1.96 * sd,
         upr = mean + 1.96 * sd)

# Remove intercept, clean up labels
feNoInt1C      <- fe1C[-1,]
feNoInt1C$term <- factor(feNoInt1C$term, 
                              levels = c("condLib",
                                         "time1",
                                         "locRifg",
                                         "locRmfg",
                                         "condLib:time1",
                                         "condLib:locRifg",
                                         "condLib:locRmfg",
                                         "time1:locRifg",
                                         "time1:locRmfg",
                                         "condLib:time1:locRifg",
                                         "condLib:time1:locRmfg"),
                              labels = c("Lib > Con",
                                         "Post > Pre",
                                         "rIFG > Sham",
                                         "rMFG > Sham",
                                         "(Lib > Con) * (Post > Pre)",
                                         "(Lib > Con) * (rIFG > Sham)",
                                         "(Lib > Con) * (rMFG > Sham)",
                                         "(Post > Pre) * (rIFG > Sham)",
                                         "(Post > Pre) * (rMFG > Sham)",
                                         "(Lib > Con) * (Post > Pre) * (rIFG > Sham)",
                                         "(Lib > Con) * (Post > Pre) * (rMFG > Sham)"))

# Plot
ggplot(feNoInt1C, aes(x = term, y = mean, ymin = lwr, ymax = upr)) +
  geom_pointrange() + 
  geom_hline(yintercept = 0, color = 'red') +
  xlab("Term") + ylab("Posterior mean") + theme_Publication() +
  coord_flip()

# Rinse and repeat for d'
fe1D <- FEsim(lmm1d, 1000)
fe1D <- fe1D %>%
  mutate(lwr = mean - 1.96 * sd,
         upr = mean + 1.96 * sd)

feNoInt1D      <- fe1D[-1,]
feNoInt1D$term <- factor(feNoInt1D$term, 
                         levels = c("condLib",
                                    "time1",
                                    "locRifg",
                                    "locRmfg",
                                    "condLib:time1",
                                    "condLib:locRifg",
                                    "condLib:locRmfg",
                                    "time1:locRifg",
                                    "time1:locRmfg",
                                    "condLib:time1:locRifg",
                                    "condLib:time1:locRmfg"),
                         labels = c("Lib > Con",
                                    "Post > Pre",
                                    "rIFG > Sham",
                                    "rMFG > Sham",
                                    "(Lib > Con) * (Post > Pre)",
                                    "(Lib > Con) * (rIFG > Sham)",
                                    "(Lib > Con) * (rMFG > Sham)",
                                    "(Post > Pre) * (rIFG > Sham)",
                                    "(Post > Pre) * (rMFG > Sham)",
                                    "(Lib > Con) * (Post > Pre) * (rIFG > Sham)",
                                    "(Lib > Con) * (Post > Pre) * (rMFG > Sham)"))

ggplot(feNoInt1D, aes(x = term, y = mean, ymin = lwr, ymax = upr)) +
  geom_pointrange() + 
  geom_hline(yintercept = 0, color = 'red') +
  xlab("Term") + ylab("Posterior mean") + theme_Publication() +
  coord_flip()

# FIGURE 8:

# Obtain posterior mean & 95% CIs
fe2C <- FEsim(lmm2, 1000)
fe2C <- fe2C %>%
  mutate(lwr = mean - 1.96 * sd,
         upr = mean + 1.96 * sd)

# Remove intercept, clean up labels
feNoInt2C      <- fe2C[-1,]
feNoInt2C$term <- factor(feNoInt2C$term, 
                         levels = c("condLib",
                                    "time1",
                                    "locRifg",
                                    "locRdlpfc",
                                    "condLib:time1",
                                    "condLib:locRifg",
                                    "condLib:locRdlpfc",
                                    "time1:locRifg",
                                    "time1:locRdlpfc",
                                    "condLib:time1:locRifg",
                                    "condLib:time1:locRdlpfc"),
                         labels = c("Lib > Con",
                                    "Post > Pre",
                                    "rIFG > Sham",
                                    "rDLPFC > Sham",
                                    "(Lib > Con) * (Post > Pre)",
                                    "(Lib > Con) * (rIFG > Sham)",
                                    "(Lib > Con) * (rDLPFC > Sham)",
                                    "(Post > Pre) * (rIFG > Sham)",
                                    "(Post > Pre) * (rDLPFC > Sham)",
                                    "(Lib > Con) * (Post > Pre) * (rIFG > Sham)",
                                    "(Lib > Con) * (Post > Pre) * (rDLPFC > Sham)"))

# Plot
ggplot(feNoInt2C, aes(x = term, y = mean, ymin = lwr, ymax = upr)) +
  geom_pointrange() + 
  geom_hline(yintercept = 0, color = 'red') +
  xlab("Term") + ylab("Posterior mean") + theme_Publication() +
  coord_flip()

# Rinse and repeat for d'
fe2D <- FEsim(lmm2d, 1000)
fe2D <- fe2D %>%
  mutate(lwr = mean - 1.96 * sd,
         upr = mean + 1.96 * sd)

feNoInt2D      <- fe2D[-1,]
feNoInt2D$term <- factor(feNoInt2D$term, 
                         levels = c("condLib",
                                    "time1",
                                    "locRifg",
                                    "locRdlpfc",
                                    "condLib:time1",
                                    "condLib:locRifg",
                                    "condLib:locRdlpfc",
                                    "time1:locRifg",
                                    "time1:locRdlpfc",
                                    "condLib:time1:locRifg",
                                    "condLib:time1:locRdlpfc"),
                         labels = c("Lib > Con",
                                    "Post > Pre",
                                    "rIFG > Sham",
                                    "rDLPFC > Sham",
                                    "(Lib > Con) * (Post > Pre)",
                                    "(Lib > Con) * (rIFG > Sham)",
                                    "(Lib > Con) * (rDLPFC > Sham)",
                                    "(Post > Pre) * (rIFG > Sham)",
                                    "(Post > Pre) * (rDLPFC > Sham)",
                                    "(Lib > Con) * (Post > Pre) * (rIFG > Sham)",
                                    "(Lib > Con) * (Post > Pre) * (rDLPFC > Sham)"))

ggplot(feNoInt2D, aes(x = term, y = mean, ymin = lwr, ymax = upr)) +
  geom_pointrange() + 
  geom_hline(yintercept = 0, color = 'red') +
  xlab("Term") + ylab("Posterior mean") + theme_Publication() +
  coord_flip()

# TABLE 2:

# Spit out model-level statistics for Exp 1 c
stargazer(lmm1, type = "html",
          covariate.labels = c("Lib > Con",
                               "Post > Pre",
                               "rIFG > Sham",
                               "rMFG > Sham",
                               "(Lib > Con) * (Post > Pre)",
                               "(Lib > Con) * (rIFG > Sham)",
                               "(Lib > Con) * (rMFG > Sham)",
                               "(Post > Pre) * (rIFG > Sham)",
                               "(Post > Pre) * (rMFG > Sham)",
                               "(Lib > Con) * (Post > Pre) * (rIFG > Sham)",
                               "(Lib > Con) * (Post > Pre) * (rMFG > Sham)"),
          out = "exp1-criterion.htm")

# Now Exp 1 d'
stargazer(lmm1d, type = "html",
          covariate.labels = c("Lib > Con",
                               "Post > Pre",
                               "rIFG > Sham",
                               "rMFG > Sham",
                               "(Lib > Con) * (Post > Pre)",
                               "(Lib > Con) * (rIFG > Sham)",
                               "(Lib > Con) * (rMFG > Sham)",
                               "(Post > Pre) * (rIFG > Sham)",
                               "(Post > Pre) * (rMFG > Sham)",
                               "(Lib > Con) * (Post > Pre) * (rIFG > Sham)",
                               "(Lib > Con) * (Post > Pre) * (rMFG > Sham)"),
          out = "exp1-dprime.htm")

# TABLE 3:

# Same as above - spit out model tables for c and d' in Exp 2
stargazer(lmm2, type = "html",
          covariate.labels = c("Lib > Con",
                               "Post > Pre",
                               "rIFG > Sham",
                               "rDLPFC > Sham",
                               "(Lib > Con) * (Post > Pre)",
                               "(Lib > Con) * (rIFG > Sham)",
                               "(Lib > Con) * (rDLPFC > Sham)",
                               "(Post > Pre) * (rIFG > Sham)",
                               "(Post > Pre) * (rDLPFC > Sham)",
                               "(Lib > Con) * (Post > Pre) * (rIFG > Sham)",
                               "(Lib > Con) * (Post > Pre) * (rDLPFC > Sham)"),
          out = "exp2-criterion.htm")

stargazer(lmm2d, type = "html",
          covariate.labels = c("Lib > Con",
                               "Post > Pre",
                               "rIFG > Sham",
                               "rDLPFC > Sham",
                               "(Lib > Con) * (Post > Pre)",
                               "(Lib > Con) * (rIFG > Sham)",
                               "(Lib > Con) * (rDLPFC > Sham)",
                               "(Post > Pre) * (rIFG > Sham)",
                               "(Post > Pre) * (rDLPFC > Sham)",
                               "(Lib > Con) * (Post > Pre) * (rIFG > Sham)",
                               "(Lib > Con) * (Post > Pre) * (rDLPFC > Sham)"),
          out = "exp2-dprime.htm")
