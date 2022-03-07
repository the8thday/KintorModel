# bayesian analysis

# a bayesian analysis is always means

library(easystats)
library(rstanarm)


# Bayesian linear regression
model <- stan_glm(Sepal.Length ~ Petal.Length, data = iris,
                  chains = 4,
                  iter = 2000,
                  warmup = 1000
                  )

parameters::parameters(model = model)
model_performance(model)

# posteriors distribution is all we want
posteriors <- describe_posterior(model)
print_md(posteriors, digits = 2) # 和频率学派的结果相比, 结果都来自后验分布

posterior2 <- insight::get_parameters(model) # 两个参数(intercept&Petal.Length)的后验分布
map_estimate(posterior2$Petal.Length) # 即是众数


# correlation & t-test ----------------------------------------------------

library(BayesFactor)

result <- correlationBF(iris$Sepal.Width, iris$Sepal.Length)

describe_posterior(result)
bayesfactor(result)

plot(bayesfactor(result)) +
  scale_fill_pizza()


# ttest
t_res <- ttestBF(formula = Sepal.Width ~ Species,
                 data = data)






