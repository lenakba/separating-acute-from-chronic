library(tidyverse)
library(dlnm)
library(slider) # for EWMA

folder_base = "my_file_location//"
l_mids_qsl = readRDS(paste0(folder_base, "my_imputed_file.rds"))

# We will make the figures based on 1 of the 5 models run on imputed data
# this is just to make things easier to code and estimate
# however, we will make the final model based on all 5 datasets (see futher down), estimated by Ruben's Rules
# The coefficients will thus be reported with the correct standard errors and CIs using this method.
d_daily_1 = l_mids_qsl[[1]]
d_daily_1 = d_daily_1 %>% group_by(player_id_num) %>% 
            mutate(day_nr = 1:n()) %>% ungroup()

# add variables for acute and overuse injuries
# we used the more proper terms sudden and gradual in the main article
d_daily_1 = d_daily_1 %>% mutate(acute = ifelse(injury == 1 & onset == "Sudden" & !is.na(onset), 1, 0),
                                 overuse = ifelse(injury == 1 & onset == "Gradual" & !is.na(onset), 1, 0))

# to report sample sizes
n_players = nrow(d_daily_1 %>% distinct(player_id_num))
n_injuries = d_daily_1 %>% count(injury)

#----------------------------- logistic regression with interaction

# define min and max lag for DLNM
# remember that day 0 (current day) 
# will be modelled separately from the 4 weeks prior
# we will lag the exposure values, so that "0" is the day before current day
lag_min = 0
lag_max = 26

# DLNM is to be run on past training load only
# we will lag the exposure data by 1 day
d_daily_1 = d_daily_1 %>% mutate(exposure_daily_lag = lag(exposure_daily))

# finding the Q matrix for DLNM
# Don't worry if there are any missing data (NA) in the past
# DLNM will model based on the data available
q_mat = tsModel::Lag(d_daily_1$exposure_daily_lag, lag_min:lag_max)

# calculate crossproducts for DLNM
# using subjectively chosen locations for knots
# based on the range of the data (hist(d_daily_1$exposure_daily))
cb_past_exposure = crossbasis(q_mat, lag=c(lag_min, lag_max), 
                              argvar = list(fun="ns", knots = c(15, 25, 80)),
                              arglag = list(fun="ns", knots = 4))

# remove days without exposure - when the players are not at risk
# note that the crossproducts have to be calculated before this step,
# and the rows in which players are not at risk removed manually from the matrix
zero_pos = which(d_daily_1$exposure_daily==0)
d_daily_exposed = d_daily_1 %>% slice(-zero_pos)
cb_past_exposure_exposed = cb_past_exposure[-zero_pos,]

# sample size after changes?
n_players_exposed = nrow(d_daily_exposed %>% distinct(player_id_num))
n_injuries_exposed = d_daily_exposed %>% count(injury)

# model with interaction (main result)
# knots are placed at same locations as for the DLNM on chronic load
# we denote the object exposureonly because it only looks at observations in which
# the players were exposed
fit_interaction_exposureonly = glm(injury ~ rms::rcs(exposure_daily, c(15, 25, 80)) + 
                        cb_past_exposure_exposed + 
                          exposure_daily*cb_past_exposure_exposed, 
                        family = "binomial", 
                        data = d_daily_exposed)

# Code for running and saving the same model with random effects using the lme4 package
# NB! This may take quite  a bit of time, depending on computing power
library(lme4)
fit_interaction_glmm = glmer(injury ~ rms::rcs(exposure_daily, c(15, 25, 80)) +
                        cb_past_exposure_exposed +
                        exposure_daily*cb_past_exposure_exposed + (1|player_id_num), family = "binomial", data = d_daily_exposed)

# save glmm object so the previous step can be skipped later
saveRDS(fit_interaction_glmm, file = paste0(folder_base, "fit_interaction_glmm.rds"))
fit_interaction_glmm = readRDS(paste0(folder_base, "fit_interaction_glmm.rds"))

# wald is only for faster computation
# use bootstrap for final model
parameters::parameters(fit_interaction_exposureonly, ci_method = "wald", exponentiate = TRUE)
parameters::parameters(fit_interaction_glmm, ci_method = "wald", exponentiate = TRUE)

#---------------------------------------- predict values

# function to fill into a matrix with the same size as the number of exposure_daily predvalues
fill_matrix = function(cb_values, n_values = 13){
  m_predvalues = matrix(nrow = n_values, ncol = ncol(cb_past_exposure))
  m_predvalues[1, 1:(ncol(cb_past_exposure))] = cb_values
  colnames(m_predvalues) = names(cb_values)
  m_predvalues %>% as_tibble() %>% fill(!!names(cb_values), .direction = "down") %>% as.matrix()
}

# helper function that uses function above
# to find examples in the observed data
find_predvalues = function(q_mat, cb_mat, exposure_preds, n_values = 13){
  
  # find examples of low, medium and high levels of chronic load
  d_qmat = q_mat %>% as_tibble()
  d_mean_exposure = rowSums(d_qmat) %>% enframe(name = NULL, value = "sum_27day")
  
  pos_zero = (which(d_mean_exposure$sum_27day == min(d_mean_exposure$sum_27day, na.rm = TRUE)))[1]
  pos_min = (which(d_mean_exposure$sum_27day == 180))[1] # this is 30 minutes, 6 days a week
  pos_median = (which(d_mean_exposure$sum_27day == median(d_mean_exposure$sum_27day, na.rm = TRUE)))[1]
  # the 75% quantile was chosen for "high" exposure
  quantiles = quantile(d_mean_exposure$sum_27day, na.rm = TRUE)
  quantile_75 = quantiles[4]
  pos_high = (which(d_mean_exposure$sum_27day == quantile_75))[1]
  
  zero_chronic = cb_mat[pos_zero, ]
  min_chronic = cb_mat[pos_min, ]
  median_chronic = cb_mat[pos_median, ]
  high_chronic = cb_mat[pos_high, ]
  
  m_zero = fill_matrix(zero_chronic)
  m_min = fill_matrix(min_chronic)
  m_median = fill_matrix(median_chronic)
  m_high = fill_matrix(high_chronic)
  
  # make the example datasets
  d_predvalues_zero = tibble(
    exposure_daily = exposure_preds,
    cb_past_exposure = m_zero
  )
  
  d_predvalues_min = tibble(
    exposure_daily = exposure_preds,
    cb_past_exposure = m_min
  )
  
  d_predvalues_median = tibble(
    exposure_daily = exposure_preds,
    cb_past_exposure = m_median
  )
  
  d_predvalues_high = tibble(
    exposure_daily = exposure_preds,
    cb_past_exposure = m_high
  )
  
  d_predvalues = bind_rows(d_predvalues_zero, d_predvalues_min, d_predvalues_median, d_predvalues_high)
  
  labels = c(rep("Zero", n_values), 
             rep("Low", n_values), 
             rep("Medium", n_values), 
             rep("High", n_values))
  
  d_predvalues = d_predvalues %>% mutate(labels = labels)
  d_predvalues
}

# function to predict values based on example data and a chosen model
pred_interaction = function(fit, d_predvalues, exposure_preds, n_values = 13){
  
d_preds_zero = predict(fit, newdata = d_predvalues %>% 
                         filter(labels == "Zero"), type = "response")
d_preds_min = predict(fit, newdata = d_predvalues %>% 
                        filter(labels == "Low"), type = "response")
d_preds_median = predict(fit, newdata = d_predvalues %>% 
                           filter(labels == "Medium"), type = "response")
d_preds_high = predict(fit, newdata = d_predvalues %>% 
                         filter(labels == "High"), type = "response")

preds = c(d_preds_zero,
          d_preds_min, 
          d_preds_median, 
          d_preds_high)

# convert to dataframe
labels = c(rep("Zero", n_values), 
           rep("Low", n_values),
           rep("Medium", n_values), 
           rep("High", n_values))

n_rep = 4
d_preds = rep(exposure_preds, n_rep) %>% 
  enframe(value = "exposure_daily", name = NULL) %>% 
  mutate(preds = preds, labels = labels)
d_preds
}

# calculate predictions
q_mat_exposed = q_mat[-zero_pos,]
exposure_preds_exposed = seq(10, 130, 10)
d_predvalues_exposed = find_predvalues(q_mat_exposed, cb_past_exposure_exposed, exposure_preds_exposed)
names(d_predvalues_exposed) = c("exposure_daily", "cb_past_exposure_exposed", "labels")
d_preds_exposed = pred_interaction(fit_basic_exposureonly, d_predvalues_exposed, exposure_preds_exposed)
d_preds_interaction_exposed = pred_interaction(fit_interaction_exposureonly, d_predvalues_exposed, exposure_preds_exposed)

# make figures
library(directlabels)
text_size = 16
# note that nih_distinct is an object from a 
# local R package with OSTRC themed colors
# coloring commands are included for transparency
# the Trebuchet MS font family has to be installed with a font package
ostrc_theme =  theme(panel.border = element_blank(), 
                     panel.background = element_blank(),
                     panel.grid = element_blank(),
                     # axis.line = element_line(color = nih_distinct[4]), 
                     strip.background = element_blank(),
                     # axis.ticks = element_line(color = nih_distinct[4]),
                     strip.text.x = element_text(size = text_size, 
                                                 family="Trebuchet MS", 
                                                 colour="black", 
                                                 face = "bold", hjust = -0.01))

plot_interaction_exposed = 
  ggplot(d_preds_interaction_exposed, aes(x = exposure_daily, y = preds, group = labels, color = labels))  + 
  geom_line(size = 0.8) +
  # scale_color_manual(values = nih_distinct[1:4]) +
  theme_line(text_size) +
  ostrc_theme +
  scale_x_continuous(breaks = scales::breaks_width(10, 0)) +
  scale_y_continuous(labels = axis_percent, limits = c(NA, 0.053), breaks = scales::breaks_width(0.01, 0)) +
  ylab("Probability of injury") +
  xlab("Minutes in activity")

devEMF::emf("plot_interaction_exposed.emf", width = 6, height = 6)
direct.label.ggplot(plot_interaction_exposed, method = list("top.bumpup", cex = 1.2))
dev.off()

setEPS()
postscript("plot_interaction_exposed.eps", width = 6, height = 6)
direct.label.ggplot(plot_interaction_exposed, method = list("top.bumpup", cex = 1.2))
dev.off()

#--------------------------------- Running analysis on all imputed datasets using Ruben's rules

# this is for final model estimates reported in the main table in the paper

l_mids_tbl_qsl = l_mids_qsl %>%  mice::complete("all") %>% 
  map(. %>% as_tibble() %>% mutate(exposure_daily_lag = lag(exposure_daily),
                                   exposure_daily_lag_imputed = ifelse(is.na(exposure_daily_lag), 
                                                                       mean(exposure_daily, na.rm = TRUE), exposure_daily_lag)))

# finding the Q matrix
l_q_mat = l_mids_tbl_qsl %>% map(~tsModel::Lag(.$exposure_daily_lag, lag_min:lag_max))

# calculate crosspredictions for DLNM
l_cb_past_exposure = l_q_mat %>% map(~crossbasis(., lag=c(lag_min, lag_max), 
                                                 argvar = list(fun="ns", knots = c(15, 25, 80)),
                                                 arglag = list(fun="ns", knots = 4)))

# removing exposure days after calculation of DLNM crosspredictions
l_zero_pos = l_mids_tbl_qsl %>% map(~which(.$exposure_daily==0))
l_daily_exposed = map2(l_mids_tbl_qsl,
                       l_zero_pos,
                       ~.x[-.y,])
l_cb_past_exposure_exposed = map2(l_cb_past_exposure,
                                  l_zero_pos,
                                  ~.x[-.y,])

l_fit_interaction = map2(l_daily_exposed,
                         l_cb_past_exposure_exposed,
                         ~glm(injury ~ rms::rcs(exposure_daily, c(15, 25, 80)) +.y + exposure_daily*.y, 
                              family = "binomial", 
                              data = .x))

d_pooled = summary(l_fit_interaction %>% mice::pool(), conf.int = TRUE, exponentiate = TRUE) %>% as_tibble() %>% mutate_if(is.numeric, ~round(., 3))

# writing model coefficients to server
#write_excel_csv(d_pooled, "modelcoefs_qsl.csv", delim = ";", na = "")

#--------------------------- to report the example profiles used in supplementary table

# find examples of low, medium and high levels of chronic load
d_qmat = q_mat_exposed %>% as_tibble()
d_mean_exposure = rowSums(d_qmat) %>% enframe(name = NULL, value = "sum_27day")

pos_zero = (which(d_mean_exposure$sum_27day == min(d_mean_exposure$sum_27day, na.rm = TRUE)))[1]
pos_min = (which(d_mean_exposure$sum_27day == 180))[1] # this is 30 minutes, 6 days a week
pos_median = (which(d_mean_exposure$sum_27day == median(d_mean_exposure$sum_27day, na.rm = TRUE)))[1]
# the 75% quantile was chosen for "high" exposure
quantiles = quantile(d_mean_exposure$sum_27day, na.rm = TRUE)
quantile_75 = quantiles[4]
pos_high = (which(d_mean_exposure$sum_27day == quantile_75))[1]

# examples
d_qmat[pos_zero,]
d_qmat[pos_min,] %>% View()
d_qmat[pos_median,] %>% View()
d_qmat[pos_high,] %>% View()

#------------------------------------ GLMM figure (ended up in supplementary)
d_predvalues_exposed_glmm = d_predvalues_exposed %>% mutate(player_id_num = "MRN-01043934") # example player
d_preds_interaction_exposed_glmm = pred_interaction(fit_interaction_glmm, d_predvalues_exposed_glmm, exposure_preds_exposed)

plot_interaction_exposed_glmm = 
  ggplot(d_preds_interaction_exposed_glmm, aes(x = exposure_daily, y = preds, group = labels, color = labels))  + 
  geom_line(size = 0.8) +
  # scale_color_manual(values = nih_distinct[1:4]) +
  theme_line(text_size) +
  ostrc_theme +
  scale_x_continuous(breaks = scales::breaks_width(10, 0)) +
  scale_y_continuous(labels = axis_percent) +
  ylab("Probability of injury") +
  xlab("Minutes in activity")

devEMF::emf("plot_interaction_exposed_glmm.emf", width = 6, height = 6)
direct.label.ggplot(plot_interaction_exposed_glmm, method = list("top.bumpup", cex = 1.2))
dev.off()

setEPS()
postscript("plot_interaction_exposed_glmm.eps", width = 6, height = 6)
direct.label.ggplot(plot_interaction_exposed_glmm, method = list("top.bumpup", cex = 1.2))
dev.off()

#---------------------------------- same analysis, same figure, but with EWMA instead
# calc number of subjects (players)
nsub = nrow(d_daily_1 %>% distinct(player_id_num))

# function for calculating exponentially waited moving averages
# using similar syntax as the RA-function
# an exponential smoothing ratio of 2/(n+1)
# same as in williams et al. 2016
ewma = function(x, n_days){
  TTR::EMA(x, n = n_days, wilder = FALSE)
}

# function calculates ewma on a sliding window of 21 days
slide_ewma = function(x){
  l = slide(x, ~ewma(., lag_max+1), .before = lag_max, step = 1, .complete = TRUE) %>% map(last)
  l = compact(l)
  l = unlist(l)
  l
}

# function to nest the exposure history data by each individual, 
# and run a user-specified function on each of their datasets in the list
function_on_list = function(d, FUN = NULL, day_start){
  nested_list = d %>% group_by(player_id_num) %>% nest()
  nested_list$data = nested_list$data %>% map(., ~FUN(.$exposure_daily_lag_imputed))
  l_unnest = unnest(nested_list, cols = c(data)) %>% group_by(player_id_num) %>% 
    mutate(day = 1:n()) %>% ungroup()
  l_unnest
}

# EWMA can't handle a single missing value, and when we lag the TL observations
# the first will be missing. We will use mean imputation for these first value so that EWMA may be calculated
# on the first 28 days of data.
d_daily_1 = d_daily_1 %>% mutate(exposure_daily_lag_chronic = lag(exposure_daily, 7),
                                 exposure_daily_lag_imputed_chronic = 
                                   ifelse(is.na(exposure_daily_lag_chronic), 
                                          mean(exposure_daily, na.rm = TRUE), exposure_daily_lag_chronic))

d_ewma = function_on_list(d_daily_1, slide_ewma, lag_max+1) %>% rename(exposure_daily_ewma = data)

# remove the first 28 rows for to be able to join the EWMA data
d_daily_mods = d_daily_1 %>% group_by(player_id_num) %>% filter(day_nr >= lag_max+1)
d_daily_mods = d_daily_mods %>% left_join(d_ewma, by = c("player_id_num", "day_nr" = "day"))

fit_ewma = glm(injury ~ rms::rcs(exposure_daily, df = c(25, 50, 100)) + 
                 exposure_daily_ewma + 
                 exposure_daily*exposure_daily_ewma, 
               family = "binomial", data = d_daily_mods %>% filter(exposure_daily != 0))

# wald is only for faster computation
# use bootstrap for final model
parameters::parameters(fit_ewma, ci_method="wald", exponentiate = TRUE)

# EWMA figure has a much easier approach with ggeffects
library(ggeffects)
# use quantiles as example levels
quantile(d_daily_mods$exposure_daily_ewma, na.rm = TRUE)

d_preds_ewma = ggpredict(fit_ewma, terms = 
                           c("exposure_daily [10:130]", 
                             "exposure_daily_ewma [0, 25, 50, 100]"))

preds_interaction_ewma = plot(d_preds_ewma, ci = FALSE) +
  # scale_color_manual(values = nih_distinct[1:4]) +
  theme_line(text_size) +
  ostrc_theme +
  scale_x_continuous(breaks = scales::breaks_width(10, 0))  +
  ylab("Probability of injury") +
  xlab("Minutes in activity") +
  theme(legend.position = "bottom") +
  ggtitle(NULL)

devEMF::emf("plot_interaction_ewma.emf", width = 6, height = 6)
preds_interaction_ewma
dev.off()

setEPS()
postscript("plot_interaction_ewma_qsl.eps", width = 6, height = 6)
preds_interaction_ewma
dev.off()

#----------------------------------Figures

# vector of tl values used in visualizations of predictions
lag_seq = lag_min:lag_max 
predvalues = seq(min(d_daily$exposure_daily), max(d_daily$exposure_daily), 10)

# predict hazards 
cp_preds = crosspred(cb_past_exposure, fit_interaction, at = predvalues, cen = 0, cumul = TRUE)

cumul_freq = ggplot(d_asympt_preds_freq_cumul, aes(x = predvalue, y = coef, group = 1)) +
  geom_ribbon(aes(min = ci_low, max = ci_high), alpha = 0.3, fill = nih_distinct[1]) +
  geom_hline(yintercept = 1, alpha = 0.3, size = 1) +
  geom_line(size = 0.75, color = nih_distinct[4]) +
  theme_base(text_size) +
  ostrc_theme +
  xlab("Daily number of jumps") +
  ylab("Cumulative HR on Day 0") 

#--------------------------------------Crossproduct predictions of chronic load

# model without removing days of no exposure
fit_interaction_alldata = glm(injury ~ cb_past_exposure, 
                              family = "binomial", 
                              data = d_daily_1)

parameters::parameters(fit_interaction_alldata, ci_method = "wald", exponentiate = TRUE)

exposure_preds_cb = seq(0, 130, 10)
pred_cb = crosspred(cb_past_exposure, fit_interaction_alldata, at = exposure_preds_cb, cen = 0, cumul = TRUE)

lag_example = pred_cb$matRRfit[7,] %>% 
  enframe(name = "Lag", value = "OR") %>% 
  mutate(Lag = as.numeric(str_replace_all(Lag, "lag", ""))+1,
         Lag = -Lag,
         low = pred_cb$matRRlow[7,],
         high = pred_cb$matRRhigh[7,])

alpha_evel = 0.3
plot_lag = ggplot(lag_example, aes(x = Lag, y = OR, ymin = low, ymax = high)) +
  geom_hline(yintercept = 1.0, col = nih_distinct[4], alpha = alpha_evel, size = 1) +
  geom_ribbon(fill = nih_distinct[1], alpha = alpha_evel) + 
  geom_line(size = 0.8, color = nih_distinct[4]) +
  theme_line(text_size) +
  ostrc_theme +
  scale_x_continuous(breaks = scales::breaks_width(5, -1)) +
  ylab("Odds Ratio") +
  xlab("Number of days since activity")

# extract confidence intervals for ORs reported in the text
# remember that 0 is the day prior to the current day
# 26 is 27 days prior to current day
# 0 is in column number 1 in the matrix, just to confuse us further.
pred_cb$matRRfit[7,1]
pred_cb$matRRlow[7,1]
pred_cb$matRRhigh[7,1]

pred_cb$matRRfit[7, 6]
pred_cb$matRRlow[7, 6]
pred_cb$matRRhigh[7, 6]

pred_cb$matRRfit[7,19]
pred_cb$matRRlow[7,19]
pred_cb$matRRhigh[7,19]

pred_cb$matRRfit[7,20]
pred_cb$matRRlow[7,20]
pred_cb$matRRhigh[7,20]

pred_cb$matRRfit[7,21]
pred_cb$matRRlow[7,21]
pred_cb$matRRhigh[7,21]

pred_cb$matRRfit[7,22]
pred_cb$matRRlow[7,22]
pred_cb$matRRhigh[7,22]

pred_cb$matRRfit[7,27]
pred_cb$matRRlow[7,27]
pred_cb$matRRhigh[7,27]

# make figures per level of activity duration in minutes

level_example = pred_cb$matRRfit[,1] %>% 
  enframe(name = "Minutes", value = "OR") %>% 
  mutate(Minutes = as.numeric(Minutes),
         low = pred_cb$matRRlow[,1],
         high = pred_cb$matRRhigh[,1])

level_example2 = pred_cb$matRRfit[,11] %>% 
  enframe(name = "Minutes", value = "OR") %>% 
  mutate(Minutes = as.numeric(Minutes),
         low = pred_cb$matRRlow[,11],
         high = pred_cb$matRRhigh[,11])

level_example3 = pred_cb$matRRfit[,27] %>% 
  enframe(name = "Minutes", value = "OR") %>% 
  mutate(Minutes = as.numeric(Minutes),
         low = pred_cb$matRRlow[,27],
         high = pred_cb$matRRhigh[,27])

plot_dlnm = function(data, day){
  
  ggplot(data, aes(x = Minutes, y = OR, group = 1, ymin = low, ymax = high))  +
    geom_hline(yintercept = 1.0, col = nih_distinct[4], alpha = alpha_evel, size = 1) +
    geom_ribbon(fill = nih_distinct[1], alpha = alpha_evel) + 
    geom_line(size = 0.8, color = nih_distinct[4]) +
    theme_line(text_size) +
    ostrc_theme +
    scale_x_continuous(breaks = scales::breaks_width(20, 0)) +
    ylab("Odds Ratio") +
    xlab(paste0("Minutes in activity (Day ",day,")"))
}

plot_level1 = plot_dlnm(level_example, "-1")
plot_level2 = plot_dlnm(level_example2, "-10")
plot_level3 = plot_dlnm(level_example3, "-27")

# eps-files don't support semi-transparency
cairo_pdf("plot_dlnm_qsl.pdf", width = 9, height = 8, fallback_resolution = 600)
ggpubr::ggarrange(plot_lag, plot_level1, plot_level2, plot_level3)
dev.off()

#----------------------------------------- Analyzing acute and overload injuries separately

# models with interaction (main result)
fit_interaction_sudden = glm(acute ~ rms::rcs(exposure_daily, c(15, 25, 80)) + 
                                     cb_past_exposure_exposed + 
                                     exposure_daily*cb_past_exposure_exposed, 
                                   family = "binomial", 
                                   data = d_daily_exposed)

fit_interaction_gradual = glm(overuse ~ rms::rcs(exposure_daily, c(15, 25, 80)) + 
                                     cb_past_exposure_exposed + 
                                     exposure_daily*cb_past_exposure_exposed, 
                                   family = "binomial", 
                                   data = d_daily_exposed)


parameters::parameters(fit_interaction_sudden, ci_method = "wald", exponentiate = TRUE)
parameters::parameters(fit_interaction_gradual, ci_method = "wald", exponentiate = TRUE)

# calculate predictions
d_preds_sudden = pred_interaction(fit_interaction_sudden, d_predvalues_exposed, exposure_preds_exposed)
d_preds_gradual = pred_interaction(fit_interaction_gradual, d_predvalues_exposed, exposure_preds_exposed)

plot_sudden = 
  ggplot(d_preds_sudden, aes(x = exposure_daily, y = preds, group = labels, color = labels))  + 
  geom_line(size = 0.8) +
  scale_color_manual(values = nih_distinct[1:4]) +
  theme_line(text_size) +
  ostrc_theme +
  scale_x_continuous(breaks = scales::breaks_width(20, 0)) +
  scale_y_continuous(labels = axis_percent, limits = c(NA, 0.04), breaks = scales::breaks_width(0.01, 0)) +
  ylab("Probability of injury") +
  xlab("Minutes in activity")

plot_gradual = 
  ggplot(d_preds_gradual, aes(x = exposure_daily, y = preds, group = labels, color = labels))  + 
  geom_line(size = 0.8) +
  scale_color_manual(values = nih_distinct[1:4]) +
  theme_line(text_size) +
  ostrc_theme +
  scale_x_continuous(breaks = scales::breaks_width(20, 0)) +
  scale_y_continuous(labels = axis_percent, limits = c(NA, 0.04), breaks = scales::breaks_width(0.01, 0)) +
  ylab("Probability of injury") +
  xlab("Minutes in activity")

plot_sudden = direct.label.ggplot(plot_sudden, method = list("top.bumpup", cex = 1.2))
plot_gradual = direct.label.ggplot(plot_gradual, method = list("top.bumpup", cex = 1.2))

setEPS()
postscript("plot_acute_vs_overuse.eps", width = 10, height = 6)
ggpubr::ggarrange(plot_sudden, plot_gradual)
dev.off()
