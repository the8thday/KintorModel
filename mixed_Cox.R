# Mixed effects cox regression models are used to model survival data when there
# are repeated measures on an individual, individuals nested within some other
# hierarchy, or some other reason to have both fixed and random effects



library(survival)
library(coxme)

m <- coxme(Surv(time, status) ~ age + sex + (1|ph.ecog), 
           lung)

