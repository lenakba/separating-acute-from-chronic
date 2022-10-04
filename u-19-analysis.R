library(tidyverse)
library(dlnm)
library(slider) # for EWMA

# read data
folder_base = paste0("O:\\Prosjekter\\Dalen-Lorentsen - Belastningsstyring - prosjekt 1 - metode\\Data\\")
folder_nor_soccer = paste0("soccer-norwegian-u19\\") # location of U19 soccer data
l_u19 = readRDS(paste0(folder_base, folder_nor_soccer,"u19_imputed.RDS"))

# define min and max lag for DLNM
# remember that day 0 (current day) 
# will be modelled separately from the 4 weeks prior
# we will lag the exposure values, so that "0" is the day before current day
lag_min = 0
lag_max = 26

# We will make the figures based on 1 of the 5 models run on imputed data
# this is just to make things easier to code and estimate
# however, we will make the final model based on all 5 datasets (see futher down), estimated by Ruben's Rules
# The coefficients will thus be reported with the correct standard errors and CIs using this method.
d_u19_1 = l_u19[[1]]
d_u19 = d_u19_1 %>% dplyr::select(p_id, date_training, load, minutes, intensity, age, sex, injury, hproblem_new)

n_players = nrow(d_u19 %>% distinct(p_id))
n_injuries = d_u19 %>% count(injury)

#--------------------------------DLNM

# DLNM is to be run on past training load only
# we will lag the exposure data by 1 day
d_u19 = d_u19 %>% mutate(load_lag = lag(load),
                         load_lag_imputed = ifelse(is.na(load_lag), mean(load, na.rm = TRUE), load))

# finding the Q matrix for DLNM
# Don't worry if there are any missing data (NA) in the past
# DLNM will model based on the data available
q_mat = tsModel::Lag(d_u19$load, lag_min:lag_max)

# calculate crossproducts for DLNM
# using subjectively chosen locations for knots
# based on the range of the data (hist(d_u19$load))
cb_past_load = crossbasis(q_mat, lag=c(lag_min, lag_max), 
                              argvar = list(fun="ns", knots = c(250, 500, 700)),
                              arglag = list(fun="ns", knots = 3))

# remove days without exposure - when the players are not at risk
# note that the crossproducts have to be calculated before this step,
# and the rows in which players are not at risk removed manually from the matrix
zero_pos = which(d_u19$minutes==0)
d_u19_exposed = d_u19 %>% slice(-zero_pos)
cb_past_load_exposed = cb_past_load[-zero_pos,]
d_u19_exposed = d_u19_exposed %>% as_tibble()

# sample size after changes?
n_players = nrow(d_u19_exposed %>% distinct(p_id))
n_injuries = d_u19_exposed %>% count(injury)
n_tl_obs = nrow(d_u19_exposed)

# model with interaction (main result)
# knots are placed at same locations as for the DLNM on chronic load
# we denote the object exposureonly because it only looks at observations in which
# the players were exposed
fit_interaction_exposureonly = glm(injury ~ rms::rcs(load, c(250, 500, 700)) + 
                                     cb_past_load_exposed + 
                                     load*cb_past_load_exposed, 
                                   family = "binomial", 
                                   data = d_u19_exposed)

# Code for running and saving the same model with random effects using the lme4 package
# NB! This may take quite  a bit of time, depending on computing power
library(lme4)
fit_interaction_glmm = glmer(injury ~ rms::rcs(load, c(250, 500, 700)) +
                               cb_past_load_exposed +
                               load*cb_past_load_exposed + (1|p_id),
                               family = "binomial", data = d_u19_exposed)

# We saved the glmm fit object as a RDS file so that we could read it later
# without having to run the fit again.
saveRDS(fit_interaction_glmm, file = paste0(folder_base, "fit_interaction_glmm_u19.rds"))
fit_interaction_glmm = readRDS(paste0(folder_base, "fit_interaction_glmm_u19.rds"))

# wald is only for faster computation
# use bootstrap for final model
parameters::parameters(fit_interaction_exposureonly, ci_method = "wald", exponentiate = TRUE)
parameters::parameters(fit_interaction_glmm, ci_method = "wald", exponentiate = TRUE)

#------------------------------------- Predictions 

# fill into a matrix with the same size as the number of exposure_daily predvalues
fill_matrix = function(cb_values, n_values = 13){
  m_predvalues = matrix(nrow = n_values, ncol = ncol(cb_past_load))
  m_predvalues[1, 1:(ncol(cb_past_load))] = cb_values
  colnames(m_predvalues) = names(cb_values)
  m_predvalues %>% as_tibble() %>% fill(!!names(cb_values), .direction = "down") %>% as.matrix()
}

# helper function that uses functions above
# to find examples and predict values
find_predvalues = function(q_mat, cb_mat, exposure_preds, n_values = 13){
  
  # find examples of low, medium and high levels of chronic load
  d_qmat = q_mat %>% as_tibble()
  d_mean_exposure = rowSums(d_qmat) %>% enframe(name = NULL, value = "sum_27day")
  
  pos_min = (which(d_mean_exposure$sum_27day == min(d_mean_exposure$sum_27day, na.rm = TRUE)))[1]
  pos_median = (which(d_mean_exposure$sum_27day == 7163))[1] # the median
  # the 75% quantile was chosen for "high" exposure
  quantiles = quantile(d_mean_exposure$sum_27day, na.rm = TRUE)
  quantile_75 = quantiles[4]
  pos_high = (which(d_mean_exposure$sum_27day == quantile_75))[1]
  
  min_chronic = cb_mat[pos_min, ]
  median_chronic = cb_mat[pos_median, ]
  high_chronic = cb_mat[pos_high, ]
  
  m_min = fill_matrix(min_chronic, n_values = n_values)
  m_median = fill_matrix(median_chronic, n_values = n_values)
  m_high = fill_matrix(high_chronic, n_values = n_values)
  
  # make the example datasets
  
  d_predvalues_min = tibble(
    load = exposure_preds,
    cb_past_load = m_min
  )
  
  d_predvalues_median = tibble(
    load = exposure_preds,
    cb_past_load = m_median
  )
  
  d_predvalues_high = tibble(
    load = exposure_preds,
    cb_past_load = m_high
  )
  
  d_predvalues = bind_rows(d_predvalues_min, d_predvalues_median, d_predvalues_high)
  
  labels = c(rep("Low", n_values), 
             rep("Medium", n_values), 
             rep("High", n_values))
  
  d_predvalues = d_predvalues %>% mutate(labels = labels)
  d_predvalues
}

pred_interaction = function(fit, d_predvalues, exposure_preds, n_values = 13){
  
  d_preds_min = predict(fit, newdata = d_predvalues %>% 
                          filter(labels == "Low"), type = "response")
  d_preds_median = predict(fit, newdata = d_predvalues %>% 
                             filter(labels == "Medium"), type = "response")
  d_preds_high = predict(fit, newdata = d_predvalues %>% 
                           filter(labels == "High"), type = "response")
  
  preds = c(d_preds_min, 
            d_preds_median, 
            d_preds_high)
  
  # convert to dataframe
  labels = c(rep("Low", n_values),
             rep("Medium", n_values), 
             rep("High", n_values))
  
  n_rep = 3
  d_preds = rep(exposure_preds, n_rep) %>% 
    enframe(value = "load", name = NULL) %>% 
    mutate(preds = preds, labels = labels)
  d_preds
}

# calculate predictions
q_mat_exposed = q_mat[-zero_pos,]
exposure_preds_exposed = seq(10, 1000, 10)
d_predvalues_exposed = find_predvalues(q_mat_exposed, cb_past_load_exposed, exposure_preds_exposed, n_values = length(exposure_preds_exposed))
names(d_predvalues_exposed) = c("load", "cb_past_load_exposed", "labels")
d_preds_exposed = pred_interaction(fit_basic_exposureonly, d_predvalues_exposed, exposure_preds_exposed, n_values = length(exposure_preds_exposed))
d_preds_interaction_exposed = pred_interaction(fit_interaction_exposureonly, d_predvalues_exposed, exposure_preds_exposed, n_values = length(exposure_preds_exposed))

#--------------------------------------------Figures
library(directlabels)
text_size = 16
# note that nih_distinct is an object from a 
# local R package with OSTRC themed colors
# running the code will result in aesthetically different figures from the article
# coloring commands are included for transparency
# the Trebuchet MS font family has to be installed with a font package
ostrc_theme =  theme(panel.border = element_blank(), 
                     panel.background = element_blank(),
                     panel.grid = element_blank(),
                     # axis.line = element_line(color = nih_distinct[4]),
                     # axis.ticks = element_line(color = nih_distinct[4])
                     strip.background = element_blank(),
                     strip.text.x = element_text(size = text_size, 
                                                 family="Trebuchet MS", 
                                                 colour="black", 
                                                 face = "bold", 
                                                 hjust = -0.01)
                     )

# theme_line is also a local package function
# with a standard OSTRC theme for line graphs
plot_basic_exposed = 
  ggplot(d_preds_exposed, aes(x = load, y = preds, group = labels, color = labels))  + 
  geom_line(size = 0.6) +
  # scale_color_manual(values = nih_distinct[1:4]) +
  # theme_line(text_size) +
  ostrc_theme 

devEMF::emf("plot_basic_exposed_u19.emf", width = 6, height = 6)
direct.label.ggplot(plot_basic_exposed, "top.bumpup")
dev.off()

plot_interaction_exposed = 
  ggplot(d_preds_interaction_exposed, aes(x = load, y = preds, group = labels, color = labels))  + 
  geom_line(size = 0.8) +
  scale_color_manual(values = nih_distinct[1:4]) +
  theme_line(text_size) +
  ostrc_theme +
  scale_y_continuous(labels = axis_percent, limits = c(NA, 0.053), breaks = scales::breaks_width(0.01, 0)) +
  ylab("Probability of injury") +
  xlab("session rating of percieved exertion (arb. u)")

devEMF::emf("plot_interaction_exposed_u19.emf", width = 6, height = 6)
direct.label.ggplot(plot_interaction_exposed, "top.bumpup")
dev.off()

setEPS()
postscript("plot_interaction_exposed_u19.eps", width = 6, height = 6)
direct.label.ggplot(plot_interaction_exposed, method = list("top.bumpup", cex = 1.2))
dev.off()

#--------------------------------- Running analysis on all imputed datasets using Ruben's rules

# this is for final model estimates reported in the main table in the paper
l_tbl_u19 = l_u19  %>% 
  map(. %>% as_tibble() %>% mutate(load_lag = lag(load),
                                   load_lag_imputed = ifelse(is.na(load_lag), mean(load, na.rm = TRUE), load)))

# finding the Q matrix
l_q_mat = l_tbl_u19 %>% map(~tsModel::Lag(.$load, lag_min:lag_max))

# calculate crosspredictions for DLNM
l_cb_past_exposure = l_q_mat %>% map(~crossbasis(q_mat, lag=c(lag_min, lag_max), 
                                                 argvar = list(fun="ns", knots = c(250, 500, 700)),
                                                 arglag = list(fun="ns", knots = 3)))

# removing exposure days after calculation of DLNM crosspredictions
l_zero_pos = l_tbl_u19 %>% map(~which(.$minutes==0))
l_daily_exposed = map2(l_tbl_u19,
                       l_zero_pos,
                       ~.x[-.y,])
l_cb_past_exposure_exposed = map2(l_cb_past_exposure,
                                  l_zero_pos,
                                  ~.x[-.y,])

l_fit_interaction = map2(l_daily_exposed,
                         l_cb_past_exposure_exposed,
                         ~glm(injury ~ rms::rcs(load, c(250, 500, 700)) + .y +  load*.y, 
                              family = "binomial", 
                              data = .x))

d_pooled = summary(l_fit_interaction %>% mice::pool(), conf.int = TRUE, exponentiate = TRUE) %>% 
  as_tibble() %>% mutate_if(is.numeric, ~round(., 3))

# writing model coefficients to server
write_excel_csv(d_pooled, "modelcoefs_u19.csv", delim = ";", na = "")

#--------------------------- to report the example profiles used
# find examples of low, medium and high levels of chronic load
d_qmat = q_mat_exposed %>% as_tibble()
d_mean_exposure = rowSums(d_qmat) %>% enframe(name = NULL, value = "sum_27day")

pos_min = (which(d_mean_exposure$sum_27day == min(d_mean_exposure$sum_27day, na.rm = TRUE)))[1]
pos_median = (which(d_mean_exposure$sum_27day == 7163))[1] # the median
# the 75% quantile was chosen for "high" exposure
quantiles = quantile(d_mean_exposure$sum_27day, na.rm = TRUE)
quantile_75 = quantiles[4]
pos_high = (which(d_mean_exposure$sum_27day == quantile_75))[1]

# examples
d_qmat[pos_min,] %>% View()
d_qmat[pos_median,] %>% View()
d_qmat[pos_high,] %>% View()

#-------------------------------------- GLMM figure (ended up in supplementary)
d_predvalues_exposed_glmm = d_predvalues_exposed %>% mutate(p_id = "101")
d_preds_interaction_exposed_glmm = pred_interaction(fit_interaction_glmm, d_predvalues_exposed_glmm, exposure_preds_exposed, n_values = length(exposure_preds_exposed))

plot_interaction_exposed_glmm = 
  ggplot(d_preds_interaction_exposed_glmm, aes(x = load, y = preds, group = labels, color = labels))  + 
  geom_line(size = 0.8) +
  # scale_color_manual(values = nih_distinct[1:4]) +
  # theme_line(text_size) +
  ostrc_theme +
  scale_y_continuous(labels = axis_percent, limits = c(NA, 0.053), breaks = scales::breaks_width(0.01, 0)) +
  ylab("Probability of injury") +
  xlab("session rating of percieved exertion (arb. u)")

# to save figure
setEPS()
postscript("plot_interaction_exposed_u19_glmm.eps", width = 6, height = 6)
direct.label.ggplot(plot_interaction_exposed_glmm, method = list("top.bumpup", cex = 1.2))
dev.off()

#---------------------------------- same analysis, same figure, but with EWMA instead
# calc number of subjects (players)
nsub = nrow(d_u19 %>% distinct(p_id))

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
  nested_list = d %>% group_by(p_id) %>% nest()
  nested_list$data = nested_list$data %>% map(., ~FUN(.$load_lag_imputed))
  l_unnest = unnest(nested_list, cols = c(data)) %>% group_by(p_id) %>% 
    mutate(day = 1:n()) %>% ungroup()
  l_unnest
}

# EWMA can't handle a single missing value, and when we lag the TL observations
# the first will be missing. We will use mean imputation for these first value so that EWMA may be calculated
# on the first 28 days of data.
d_u19 = d_u19 %>% mutate(load_lag_chronic = lag(load, 7),
                         load_lag_imputed_chronic = 
                                   ifelse(is.na(load_lag_chronic), 
                                          mean(load, na.rm = TRUE), load_lag_chronic))

d_ewma = function_on_list(d_u19, slide_ewma, lag_max+1) %>% rename(load_ewma = data)

# remove the first 28 rows for to be able to join the EWMA data
d_u19_mods = d_u19 %>% group_by(p_id) %>% mutate(day_nr = 1:n()) %>% filter(day_nr >= lag_max+1)
d_u19_mods = d_u19_mods %>% left_join(d_ewma, by = c("p_id", "day_nr" = "day"))

fit_ewma = glm(injury ~ rms::rcs(load, df = c(250, 500, 700)) + 
                 load_ewma + 
                 load*load_ewma, 
               family = "binomial", data = d_u19_mods %>% filter(minutes != 0))

# wald is only for faster computation
# use bootstrap for final model
parameters::parameters(fit_ewma, ci_method="wald", exponentiate = TRUE)

# EWMA figure has a much easier approach with ggeffects
library(ggeffects)
# look at quantiles for example levels
quantile(d_u19_mods$load_ewma, na.rm = TRUE)

d_preds_ewma = ggpredict(fit_ewma, terms = 
                           c("load [10:800]", 
                             "load_ewma [0, 100, 200, 300]"))

preds_interaction_ewma = plot(d_preds_ewma, ci = FALSE) +
  # scale_color_manual(values = nih_distinct[1:4]) +
  # theme_line(text_size) +
  ostrc_theme  +
  ylab("Probability of injury") +
  xlab("session rating of percieved exertion (arb. u)") +
  ggtitle(NULL)

setEPS()
postscript("plot_interaction_ewma_u19.eps", width = 6, height = 6)
direct.label.ggplot(preds_interaction_ewma, method = list("top.bumpup", cex = 1.2))
dev.off()
