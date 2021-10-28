# SVM by tidymodels

library(tidymodels)


# fake data ---------------------------------------------------------------

set.seed(1)
sim_data <- tibble(
  x1 = rnorm(40),
  x2 = rnorm(40),
  y  = factor(rep(c(-1, 1), 20))
) %>%
  mutate(x1 = ifelse(y == 1, x1 + 1.5, x1),
         x2 = ifelse(y == 1, x2 + 1.5, x2))

set.seed(2)
sim_data_test <- tibble(
  x1 = rnorm(20),
  x2 = rnorm(20),
  y  = factor(rep(c(-1, 1), 10))
) %>%
  mutate(x1 = ifelse(y == 1, x1 + 1.5, x1),
         x2 = ifelse(y == 1, x2 + 1.5, x2))

# initate data ------------------------------------------------------------

data(ames)
report::report_table(ames$Sale_Price)
DataExplorer::plot_missing(ames)
DataExplorer::create_report(ames)

set.seed(42)
df_split <- initial_split(ames, prop = 0.8)
df_train <- training(df_split)
df_test <- testing(df_split)

glue::glue('Training data shape: {dim(df_train)}')


# model -------------------------------------------------------------------

# parsnip_addin() could help you generate the model
svm_C <- svm_poly(cost = tune(), degree = tune()) %>% 
  set_mode('regression') %>% 
  set_engine('kernlab')

# recipes for feature engineering
svm_recipes <- 
  recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type,
         data = df_train) %>%
  step_log(Gr_Liv_Area, base = 10) %>% 
  step_dummy(all_nominal_predictors())
tidy(svm_recipes)

# 交叉验证，网格搜索超参数
svm_wf <- workflow() %>% 
  add_model(svm_C %>% set_args(cost = tune())) %>% 
  add_recipe(svm_recipes)

set.seed(42)
sim_data_fold <- vfold_cv(
  df_train,
  v = 10,
  strata = Sale_Price
)

param_grid <- grid_regular(cost(), levels = 10)

tune_res <- tune_grid(
  svm_wf,
  resamples = sim_data_fold,
  grid = param_grid
)
autoplot(tune_res)

best_param <- select_best(tune_res,
                          metric = 'accuracy'
                          )

svm_C_final <- finalize_workflow(svm_wf, best_param)

svm_C_fit <- svm_C_final %>% fit(sim_data) # this is the final model


# 如何查看模型自身的参数和performance
svm_C_fit %>% 
  purrr::pluck('fit')
svm_C_fit %>% 
  predict(new_data = sim_data_test)
svm_C_fit %>% 
  extract_fit_engine() %>% 
  plot()

# 在test data performance
augment(svm_C_fit, 
        new_data = sim_data_test
        ) %>% 
  conf_mat(truth = y, estimate = .pred_class)

augment(svm_C_fit, 
        new_data = sim_data_test) %>%
  roc_curve(truth = y, estimate = .pred_1) %>% 
  autoplot()

svm_C_fit %>% predict(sim_data_test)





