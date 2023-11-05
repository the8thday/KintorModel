# 罗辑回归模型，以及模型评价

library(tidyverse)
library(tidymodels)
library(foreign)
library(rms)


# tidy data ---------------------------------------------------------------

mydata <- read.spss('/Users/congliu/OneDrive/kintor/Clinic_Predict_Model_R_Code_Data/Lowweight.sav') %>%
  as.data.frame()

# 对于二分类数据，一般赋值较小者作为参照
mydata$low <- ifelse(mydata$low =="low weight",1,0)
# 对race 哑变量设置，哑变量的选择应该注意referrnce的选择
mydata$race1 <- ifelse(mydata$race =="white",1,0)
mydata$race2 <- ifelse(mydata$race =="black",1,0)
mydata$race3 <- ifelse(mydata$race =="other",1,0)

# 逻辑回归满足的条件(LINE)
performance::check_collinearity(fit2)

# use rms package ---------------------------------------------------------

attach(mydata)
dd <- datadist(mydata)
options(datadist='dd')

fit <- rms::lrm(
  low ~ age+ftv+ht+lwt+ptl+smoke+ui+race1+race2,
  data = mydata,
  x = T,
  y = T
)
summary(fit)

fit2 <- glm(
  low ~ age+ftv+ht+lwt+ptl+smoke+ui+race1+race2,
  data = mydata,
  family = binomial,
  x = TRUE
)
summary(fit2)

# 特征选择 --------------------------------------------------------------------
# 特征选择在模型的构建 尤其临床模型构建中更是重要


# 评价模型 --------------------------------------------------------------------

# 对于有train or test 数据集的模型而言，可以对test数据集进行评价
prob <- fit %>% predict(test.data,
                        type = 'response'
                        )
# confusion matrix & ROC curve
library(ROCR)

mydata$predvalue <- predict(fit,
                            type = 'response'
                            )

pred <- prediction(mydata$predvalue, mydata$low)
perf <- performance(pred,
                    measure = "tpr",
                    x.measure = "fpr")
plot(perf, colorize=TRUE,
     main = 'ROC Curve'
     )
abline(0,1, col = 3, lty = 2)
plot(perf,
     lty=3,
     col="grey78",
     add=TRUE)
auc <- performance(pred,"auc")
auc@y.values

# for C-statistic
Hmisc::somers2(mydata$predvalue, mydata$low)
# confuse matrix
yardstick::conf_mat()
p <- predict(glm_model, test, type = "response")
hd_or_nohd <- ifelse(p > 0.5, 1, 0)
p_class <- factor(hd_or_nohd, levels = levels(test[["hd"]]))
InformationValue::confusionMatrix(p_class, test[["hd"]])


# 净重新分类指数(net reclassification improvement NRI)
# 在新旧指标，或者在不同模型上的重新分类的变化，似乎应该只适合两种情况下的对比使用
library(nricens)

fit3 <- glm(
  low ~ age+ftv+ht+lwt+ptl+smoke,
  data = mydata,
  family = binomial,
  x = TRUE
)

p.std <- fit2$fitted.values
p.new <- fit3$fitted.values

nribin(mdl.std = fit2,
       mdl.new = fit3,
       updown = 'category',
       cut = c(0.2, 0.4), # cut 为高低风险的临界值
       niter = 1000
       )

# 综合判别改善指数(integrated discrimination improvement IDI)
# PredictABEL::reclassification(
#   data = mydata,
#   cOutcome = '',
#   predrisk1 = p.std,
#   predrisk2 = p.new,
#   cutoff = c(0,0.2,0.4,1)
# )

# 以上两个指标的计算中NRI/IDI 需要的是两个模型间的比较
# 提取模型参数 ------------------------------------------------------------------

summary(fit)


# nomograph ----------------------------------------------------------------------
nom <- nomogram(fit = fit,
                fun = plogis,
                fun.at = c(.001, .01, .05,
                           seq(.1,.9, by = .1), .95, .99, .999),
                lp = F,
                funlabel = 'LowWeight rate'
                )
plot(nom)


# 校准曲线--------------------------------------
# 校准曲线用于评价模型校准度，即是实际发生概率和预测发生概率的散点图
# 用以判断预测发生的概率和实际发生概率的一致性, 罗辑回归  拟合优度检验
cal <- calibrate(fit = fit,
                 method = 'boot',
                 B = 1000
                 )
cal
plot(cal,xlim = c(0,1.0),ylim = c(0,1.0))

plot(p,
     add=F,
     conf.int=T,#95%CI（蓝色线）
     subtitles = F,#关闭副标题
     cex.subtitles=0.8,
     lwd=2,
     lty=1,
     errbar.col="blue",
     xlim=c(0.25,0.4),#调节x.y轴刻度范围
     ylim=c(0.25,0.4),
     xlab="列线图预测的5年OS",
     ylab="实际5年OS",
     col="red")


# decision curve analysis -------------------------------------------------

# 用以比对不同模型间的区分度等, 用以比对不同的模型优良度

library(rmda)

comp1 <- decision_curve(
  low ~age+ftv+ht,
  data = mydata,
  family = binomial(link = 'logit'),
  thresholds = seq(0, 1, by=0.01),
  confidence.intervals = 0.95,
  study.design = 'cohort'
)
comp2 <- decision_curve(
  low ~age+ftv+ht+lwt+ptl+smoke+ui+race1+race2,
  data = mydata,
  family = binomial(link = 'logit'),
  thresholds = seq(0, 1, by=0.01),
  confidence.intervals = 0.95,
  study.design = 'cohort'
)

models <- list(comp1, comp2)

plot_decision_curve(models,
                    curve.names = c('comp1', 'comp2'),
                    cost.benefit.axis = FALSE,
                    col = c('red', 'blue'),
                    confidence.intervals = FALSE,
                    standardize = TRUE
                    )
# 也可以只针对一个样本
plot_decision_curve(
  comp2
)

# 查看模型净收益率
summary(comp2, measure = 'sNB')

# 临床影响曲线(clinical impact curve)
plot_clinical_impact(
  comp1,
  population.size = 1000,
  cost.benefit.axis = TRUE,
  n.cost.benefits = 6,
  confidence.intervals = TRUE,
  col = c('red', 'blue')
)


# 用y叔叔的包, 据说碾压rmda
library(ggDCA)
d_res <- dca(
  fit, fit2
)

ggplot(d_res)

# 作用于新的数据集



# 外部验证 --------------------------------------------------------------------

# 外部验证是指用非建模时的数据，进行验证
library(ResourceSelection)
library(PredictABEL)

# 校准度评价,和rms包的calibrate函数应该是一致的
h1 <- hoslem.test(fit2$y,
                  fitted(fit2),
                  g=10
                  )
h1 # P值<0.05，模型拟合不良
cbind(h1$observed, h1$expected) # 生成Hosmer-Lemeshow检验列联表

pr.e <- predict(
  mod,
  exter.data,
  type=c('response')
)
hl.e <- hoslem.test(y, pr.e, g=10)
# 依据其的P值进行判断，如P值<0.05, 则说明模型拟合不良

# 使用校准曲线图 直观的评价模型
plotCalibration(data = exter.data,
                cOutcome = 2,# 结局变量所在的列
                predRisk = pr.e,
                groups = 10,
                rangeaxis = c(0,1)
                )

# 区分度评估
# 所谓区分度是指模型将二分类结果正确区分的能力，灵敏度特异性等指标
# ROC曲线和AUC值即是常用的指标

pr <- predict(mod, type=c('response'))
roccurve <- pROC::roc(y ~ pr)

plot.roc(roccurve,
         xlim = c(1, 0),
         ylim = c(0, 1)
         )

# 对于外部数据的ROC曲线
roccurve <- roc(y ~ pr.e)
plot.roc(roccurve)
auc(roccurve)
















