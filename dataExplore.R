# 处理数据的NA值，是常见的操作。在此做一个总结

library(mice)
library(ggmice)
library(VIM)
library(visdat)
library(DataExplorer)
library(janitor) # 较为方便整理数据的R包


# 数据探索 --------------------------------------------------------------------

data(cancer, package = "survival")
# colon 数据集，以下为html报告
DataExplorer::create_report(colon,
                            output_file = 'test_report',
                            output_dir = getwd()
                            )

# 亦可用report包得到文字描述性统计结果
report::report(colon)

# explore 包同样是个极好的交互数据探索性包
# DescrTab2 包可生成表格


# fill na -----------------------------------------------------------------
# mice包可以方便的统计NA值的数目
mice::md.pattern(colon)
# NA plot
visdat::vis_miss(colon)

# fill NA 需要依据具体情况做出具体的填充，比如均值，or 上下值



















