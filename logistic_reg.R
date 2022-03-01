#! /urs/bin/env Rscript
# Title     : logistic
# Objective : logistic model by easystats
# Created by: congliu
# Created on: 2021/8/25

library(tidymodels)

# data
data(bivariate, package = 'modeldata')


# model & feature enginer
lg_model <- logistic_reg(
  mode = 'classification',
  engine = 'glm'
)

lg_recipe <- recipe()

lg_wf <- workflow() %>%
  add_model(lg_model) %>%
  add_formula(Class ~.)

lg_fit <- lg_wf %>% fit(data = bivariate_train)

# check model
names(lg_fit)
lg_fit %>% performance::check_model()
lg_fit$fit$linear.predictors
# Box–Tidwell
car::boxTidwell(lg_fit)
multi_metric <- metric_set()
tidy(lg_fit)
# odds ratio
exp(coef(lg_fit$fit$fit$fit))
# C-statistic
train_res <- stats::predict(lg_fit$fit$fit$fit, type = 'response') %>%
  data.frame() %>%
  bind_cols(bivariate_train) %>% rename('predvalue'='.') %>%
  mutate(Class = if_else(Class=='One', 0, 1))
Hmisc::somers2(train_res$predvalue, train_res$Class)

# test data
test_res <- predict(lg_fit, new_data = bivariate_test) %>%
  bind_cols(bivariate_test) %>%
  bind_cols(predict(lg_fit, new_data = bivariate_test, type = 'prob'))
# test_res <- augment(lg_fit, new_data = bivariate_test)
conf_mat(test_res, truth = Class, estimate = .pred_class) %>%
  autoplot(type = 'heatmap')
f_meas(test_res, truth = Class, estimate = .pred_class)
autoplot(roc_curve(test_res, Class, .pred_One))
roc_auc(test_res, Class, .pred_One)
pROC::ci.auc(response = test_res$Class, predictor = test_res$.pred_One)
pROC::ci.auc(pROC::auc(test_res$Class, test_res$.pred_One))
