
# IC50的计算过程大概为DMSO组和药物组减去blank后，取药物组/DMSO组的比值。



# 四参数对于第一列是Blank的数据 -------------------------------------------------------------



all_sheets <- readxl::excel_sheets("~/OneDrive/kintor/Daily_Work/IC50_xiaodan/20221104  hPBMCs top40  13.xlsx")
conn_file <- readxl::read_excel("~/OneDrive/kintor/Daily_Work/IC50_xiaodan/concen_file.xlsx")



ic50_pipe_4p_1blank <- function(index,
                      file = "~/OneDrive/kintor/Daily_Work/IC50_xiaodan/20221104  hPBMCs top40  13.xlsx",
                      outpath = "~/OneDrive/kintor/Daily_Work/IC50_xiaodan/IC50_patch") {
  target <- all_sheets[[index]]
  print(target)

  if (length(str_split(target, pattern = " +")[[1]]) == 1) {
    message("Only one drug in this Sheet!!")
    t1 <- str_trim(target)


    all_raw_data <- readxl::read_excel(
      file,
      sheet = target
    ) %>%
      janitor::clean_names() %>%
      janitor::remove_empty(which = 'rows')

    all_raw_data %<>% mutate(group = rep(t1, each = 3)) %>%
      dplyr::select(-`temperature_c`)

    df1 <- all_raw_data %>%
      group_by(group) %>%
      mutate(
        group_mean = mean(`x1`),
        dmso_mean = mean(`x3`)
      ) %>%
      ungroup() %>%
      mutate(across(.cols = `x4`:`x11`, ~ (. - group_mean) / (dmso_mean - group_mean)))

    drug1_file <- df1 %>%
      # slice(1:3)
      slice_head(n = 3)

    cal_ic50(drug1_file, t1,
             outpath = outpath
    )

  } else {
    t1 <- str_split(target, pattern = " +")[[1]][1]
    t2 <- str_split(target, pattern = " +")[[1]][2]


    all_raw_data <- readxl::read_excel(
      file,
      sheet = target
    ) %>%
      janitor::clean_names() %>%
      janitor::remove_empty()

    all_raw_data %<>% mutate(group = rep(c(t1, t2), each = 3)) %>%
      dplyr::select(-`temperature_c`)

    df1 <- all_raw_data %>%
      group_by(group) %>%
      mutate(
        group_mean = mean(`x1`),
        dmso_mean = mean(`x3`)
      ) %>%
      ungroup() %>%
      mutate(across(.cols = `x4`:`x11`, ~ (. - group_mean) / (dmso_mean - group_mean)))

    drug1_file <- df1 %>%
      # slice(1:3)
      slice_head(n = 3)

    drug2_file <- df1 %>%
      # slice(1:3)
      slice_tail(n = 3)

    cal_ic50(drug1_file, t1,
             outpath = outpath
    )
    cal_ic50(drug2_file, t2,
             outpath = outpath
    )
  }

}


cal_ic50 <- function(drug_file, drug_name, outpath) {
  d_long <- drug_file %>%
    dplyr::select(x4:x11) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    pivot_longer(
      cols = -c(rowname),
      names_to = "rep",
      values_to = "response"
    ) %>%
    rename("con" = "rowname") %>%
    left_join(conn_file, by = c("con" = "id"))

  fitted_curve <- drc::drm(
    formula = response ~ drug_con,
    data = d_long,
    fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50"))
    # fct = LL.3()
  )

  E0 <- fitted_curve$coefficients[2]
  Einf <- fitted_curve$coefficients[3]
  EC50 <- fitted_curve$coefficients[4]
  H <- -fitted_curve$coefficients[1]
  # IC50 <- exp(log((((Einf - E0) / (50 - E0) ) - 1) * EC50**(-H))/ (-H))
  coefs <- setNames(
    fitted_curve$coefficients,
    c("hill", "min_value", "max_value", "ec_50")
  )

  IC50 <- with(
    as.list(coefs),
    exp(
      log(ec_50) + (1 / hill) * log(max_value / (max_value - 2 * min_value))
    )
  )

  efficacy_metrics <- data.frame(c(E0, Einf, EC50, H, IC50))
  rownames(efficacy_metrics) <- c("E0", "Einf", "EC50", "H", "IC50")
  colnames(efficacy_metrics) <- "Value"

  efficacy_metrics %>%
    rownames_to_column() %>%
    mutate(drug = drug_name) %>%
    write_excel_csv(file.path(outpath, glue::glue("{drug_name}_4P.csv")))

  plot_ic50(fitted_curve, d_long, outpath, drug_name)
}


plot_ic50 <- function(fitted_curve, d_long, outpath, drug_name) {
  newdata <- expand.grid(conc = exp(seq(log(0.01), log(100), length = 100)))

  pm <- predict(fitted_curve, newdata = newdata, interval = "confidence")
  # new data with predictions
  newdata$p <- pm[, 1]
  newdata$pmin <- pm[, 2]
  newdata$pmax <- pm[, 3]

  d_long$conc0 <- d_long$drug_con
  # d_p$conc0 <- log(d_p$drug_con)
  d_long$conc0[d_long$conc0 == 0] <- 0.001
  # plotting the curve
  p1 <- ggplot(d_long, aes(x = conc0, y = response)) +
    geom_ribbon(data = newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax), alpha = 0.2) +
    geom_line(data = newdata, aes(x = conc, y = p)) +
    scale_x_log10() +
    xlab(drug_name) +
    ylab("Response") +
    scale_colour_manual(
      values = c("#9FA3FE", "#00167B", "#9FA3FE")
    ) +
    ggnewscale::new_scale_colour() +
    geom_point(size = 3, color = "#9FA3FE") +
    scale_colour_prism(
      palette = "winter_bright",
    ) +
    scale_shape_prism() +
    theme_prism(palette = "winter_bright", base_size = 16)

  ggsave(
    filename = file.path(outpath, glue::glue("{drug_name}_4P.pdf")),
    width = 8,
    height = 8
  )
}



### start here
for (i in seq_along(head(all_sheets, n = 2))) {
  print(i)

  ic50_pipe(i)
}



# 四参数对于第二列是Blank数据的 ----------------------------------------------------------


all_sheets <- readxl::excel_sheets("~/OneDrive/kintor/Daily_Work/IC50_xiaodan/20221104  hPBMCs top40  13.xlsx")
conn_file <- readxl::read_excel("~/OneDrive/kintor/Daily_Work/IC50_xiaodan/concen_file.xlsx")


ic50_pipe_4p_2blank <- function(index,
                        file = "~/OneDrive/kintor/Daily_Work/IC50_xiaodan/20221104  hPBMCs top40  13.xlsx",
                        outpath = "~/OneDrive/kintor/Daily_Work/IC50_xiaodan/IC50_patch") {
  target <- all_sheets[[index]]
  print(target)

  if (length(str_split(target, pattern = " +")[[1]]) == 1) {
    message("Only one drug in this Sheet!!")

    t1 <- str_trim(target)

    all_raw_data <- readxl::read_excel(
      file,
      sheet = target
    ) %>%
      janitor::clean_names() %>%
      janitor::remove_empty(which = "rows")

    all_raw_data %<>% mutate(group = rep(t1, each = 3)) %>%
      dplyr::select(-`temperature_c`)

    df1 <- all_raw_data %>%
      group_by(group) %>%
      mutate(
        group_mean = mean(`x2`),
        dmso_mean = mean(`x3`)
      ) %>%
      ungroup() %>%
      mutate(across(.cols = `x4`:`x11`, ~ (. - group_mean) / (dmso_mean - group_mean)))

    drug1_file <- df1 %>%
      # slice(1:3)
      slice_head(n = 3)

    cal_ic50(drug1_file, t1,
             outpath = outpath
    )

  } else {

    t1 <- str_split(target, pattern = " +")[[1]][1]
    t2 <- str_split(target, pattern = " +")[[1]][2]


    all_raw_data <- readxl::read_excel(
      file,
      sheet = target
    ) %>%
      janitor::clean_names() %>%
      janitor::remove_empty(which = "rows")

    all_raw_data %<>% mutate(group = rep(c(t1, t2), each = 3)) %>%
      dplyr::select(-`temperature_c`)

    df1 <- all_raw_data %>%
      group_by(group) %>%
      mutate(
        group_mean = mean(`x2`),
        dmso_mean = mean(`x3`)
      ) %>%
      ungroup() %>%
      mutate(across(.cols = `x4`:`x11`, ~ (. - group_mean) / (dmso_mean - group_mean)))

    drug1_file <- df1 %>%
      # slice(1:3)
      slice_head(n = 3)

    drug2_file <- df1 %>%
      # slice(1:3)
      slice_tail(n = 3)

    cal_ic50(drug1_file, t1,
             outpath = outpath
    )
    cal_ic50(drug2_file, t2,
             outpath = outpath
    )
  }

}


cal_ic50 <- function(drug_file, drug_name, outpath) {
  d_long <- drug_file %>%
    dplyr::select(x4:x11) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    pivot_longer(
      cols = -c(rowname),
      names_to = "rep",
      values_to = "response"
    ) %>%
    rename("con" = "rowname") %>%
    left_join(conn_file, by = c("con" = "id"))

  fitted_curve <- drc::drm(
    formula = response ~ drug_con,
    data = d_long,
    fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50"))
    # fct = LL.3()
  )

  E0 <- fitted_curve$coefficients[2]
  Einf <- fitted_curve$coefficients[3]
  EC50 <- fitted_curve$coefficients[4]
  H <- -fitted_curve$coefficients[1]
  # IC50 <- exp(log((((Einf - E0) / (50 - E0) ) - 1) * EC50**(-H))/ (-H))
  coefs <- setNames(
    fitted_curve$coefficients,
    c("hill", "min_value", "max_value", "ec_50")
  )

  IC50 <- with(
    as.list(coefs),
    exp(
      log(ec_50) + (1 / hill) * log(max_value / (max_value - 2 * min_value))
    )
  )

  efficacy_metrics <- data.frame(c(E0, Einf, EC50, H, IC50))
  rownames(efficacy_metrics) <- c("E0", "Einf", "EC50", "H", "IC50")
  colnames(efficacy_metrics) <- "Value"

  efficacy_metrics %>%
    rownames_to_column() %>%
    mutate(drug = drug_name) %>%
    write_excel_csv(file.path(outpath, glue::glue("{drug_name}_4P.csv")))

  plot_ic50(fitted_curve, d_long, outpath, drug_name)
}



### start here
for (i in seq_along(all_sheets)) {
  print(i)

  ic50_pipe_2(i)
}




# 三参数 ---------------------------------------------------------------------

# Excel中的sheet名字
filepath1 = '~/OneDrive/kintor/Daily_Work/IC50_xiaodan/HK-2 20221111_wa.xlsx'
outpath1 = '~/OneDrive/kintor/Daily_Work/IC50_xiaodan/IC50_patch/wa_hk2_1111'

all_sheets <- readxl::excel_sheets(filepath1)
conn_file <- readxl::read_excel("~/OneDrive/kintor/Daily_Work/IC50_xiaodan/concen_file2.xlsx")


cal_ic50_3p <- function(drug_file, drug_name, outpath, cellname,
                        drug_columns = x4:x11
                        ) {

  drug_name = drug_name

  d_long <- drug_file %>%
    dplyr::select({{ drug_columns }}) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    pivot_longer(
      cols = -c(rowname),
      names_to = "rep",
      values_to = "response"
    ) %>%
    rename("con" = "rowname") %>%
    left_join(conn_file, by = c("con" = "id"))


  fitted_curve <- drc::drm(
    formula = response ~ drug_con,
    data = d_long,
    fct = LL.4(fixed = c(1, NA, NA, NA))
    # fct = LL.3()
  )

  E0_2 <- fitted_curve$coefficients[1]
  Einf_2 <- fitted_curve$coefficients[2]
  EC50_2 <- fitted_curve$coefficients[3]


  efficacy_metrics <- data.frame(c(E0_2, Einf_2, EC50_2))
  rownames(efficacy_metrics) <- c("E0", "Einf", "EC50")
  colnames(efficacy_metrics) <- "Value"

  efficacy_metrics2 <- broom::tidy(fitted_curve) %>%
    bind_cols(tibble(rowname = c('E0', 'Einf', 'EC50'))) %>%
    dplyr::relocate(rowname, everything()) %>%
    bind_cols(confint(fitted_curve))

  efficacy_metrics2 %>%
    mutate(drug = drug_name) %>%
    write_excel_csv(file.path(outpath, glue::glue("{drug_name}_3P.csv")))

  d_long %>% pivot_wider(id_cols = drug_con,
                         values_from = response,
                         names_from = rep
                         ) %>%
    write_excel_csv(file.path(outpath, glue::glue("{drug_name}_rawdata.csv")))

  plot_ic50_3p(fitted_curve, d_long, outpath, drug_name,efficacy_metrics2,cellname)

  plot_ic50_errorbar(d_long,efficacy_metrics2,drug_name,cellname,outpath)
}


plot_ic50_3p <- function(fitted_curve, d_long, outpath, drug_name, efficacy_metrics,
                         cellname) {

  ic50_value <- round(efficacy_metrics %>% filter(rowname == 'EC50') %>%
                        pull(estimate), 3)

  newdata <- expand.grid(conc = exp(seq(log(0.01), log(100), length = 100)))

  pm <- predict(fitted_curve, newdata = newdata, interval = "confidence")
  # new data with predictions
  newdata$p <- pm[, 1]
  newdata$pmin <- pm[, 2]
  newdata$pmax <- pm[, 3]

  d_long$conc0 <- d_long$drug_con
  # d_p$conc0 <- log(d_p$drug_con)
  d_long$conc0[d_long$conc0 == 0] <- 0.001
  # plotting the curve
  p1 <- ggplot(d_long, aes(x = conc0, y = response)) +
    geom_ribbon(data = newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax), alpha = 0.2) +
    geom_line(data = newdata, aes(x = conc, y = p)) +
    scale_x_log10() +
    xlab(glue::glue('{drug_name} {cellname}')) +
    ylab("Response") +
    ggplot2::annotate(
      geom = "text", x = 0.2, y = 0.2, size = 8,
      label = glue::glue("IC50 = {ic50_value}")
    ) +
    scale_colour_manual(
      values = c("#9FA3FE", "#00167B", "#9FA3FE")
    ) +
    ggnewscale::new_scale_colour() +
    geom_point(size = 3, color = "#9FA3FE") +
    scale_colour_prism(
      palette = "winter_bright",
    ) +
    scale_shape_prism() +
    theme_prism(palette = "winter_bright", base_size = 16)

  ggsave(plot = p1,
    filename = file.path(outpath, glue::glue("{drug_name}_3P.pdf")),
    width = 8,
    height = 8
  )
}



plot_ic50_errorbar <- function(d_long,efficacy_metrics,drug_name,cellname,outpath){

  data_1 <- aggregate(d_long["response"], list(conc = d_long$drug_con), FUN = mean)
  data_1$sd <- aggregate(d_long["response"], list(conc = d_long$drug_con), FUN = sd)$response
  ic50_value <- round(efficacy_metrics %>% filter(rowname == 'EC50') %>%
                        pull(estimate), 3)


  p2 <- ggplot(
    data = data_1,
    aes(x = conc, y = response)
  ) +
    geom_point(col = "red", size = 3, shape = 21, fill = "yellowgreen", position = "identity") +
    geom_errorbar(aes(ymin = response - sd, ymax = response + sd), width = .05, position = "identity") +
    geom_smooth(method = drm, col = "skyblue", method.args = list(fct = L.4()), se = F) +
    scale_x_log10() +
    ggplot2::annotate(
      geom = "text", x = 0.2, y = 0.2, size = 8,
      label = glue::glue("IC50 = {ic50_value}")
    ) +
    labs(
      title = glue::glue("{drug_name}  {cellname}"),
      x = "Log Drug Concentration(nM)", y = "Relative Cell viability",
      caption = 'Created by DTM'
    ) +
    theme(
      legend.position = "",
      panel.background = element_blank(),
      # panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, color = "darkblue", family = "sans", hjust = 0.5),
      axis.text.y = element_text(size = 12, color = "darkblue", family = "sans", vjust = 0.5, hjust = 0.5),
      axis.title = element_text(size = 16, color = "darkred", family = "sans"),
      axis.ticks = element_line(size = 1),
      axis.line = element_line(size = 1),
      axis.ticks.length = unit(3, "pt"),
      plot.title = element_text(size = 20, colour = "black", hjust = 0.5, vjust = -1, face = "bold")
    )

  ggsave(plot = p2,
    filename = file.path(outpath, glue::glue("{drug_name}_3P_errorbar.pdf")),
    width = 8,
    height = 8
  )
}




ic50_pipe_3p_1blank <- function(index,
                         file = filepath1,
                         outpath = outpath1,
                         cellname = 'hPBMC'
                         ) {

  target <- all_sheets[[index]]
  print(target)

  if (length(str_split(target, pattern = " +")[[1]]) == 1) {
    t1 <- str_trim(target)

    all_raw_data <- readxl::read_excel(
      file,
      sheet = target
    ) %>%
      janitor::clean_names() %>%
      janitor::remove_empty()

    all_raw_data %<>% mutate(group = rep(t1, each = 3)) %>%
      dplyr::select(-`temperature_c`)

    df1 <- all_raw_data %>%
      group_by(group) %>%
      mutate(
        group_mean = mean(`x1`),
        dmso_mean = mean(`x3`)
      ) %>%
      ungroup() %>%
      mutate(across(.cols = `x4`:`x11`, ~ (. - group_mean) / (dmso_mean - group_mean)))

    drug1_file <- df1 %>%
      # slice(1:3)
      slice_head(n = 3)

    cal_ic50_3p(drug_file = drug1_file, drug_name = t1,
      outpath = outpath, cellname = cellname
    )
  } else {
    t1 <- str_split(target, pattern = " +")[[1]][1]
    t2 <- str_split(target, pattern = " +")[[1]][2]


    all_raw_data <- readxl::read_excel(
      file,
      sheet = target
    ) %>%
      janitor::clean_names() %>%
      janitor::remove_empty(which = 'rows')

    all_raw_data %<>% mutate(group = rep(c(t1, t2), each = 3)) %>%
      dplyr::select(-`temperature_c`)

    df1 <- all_raw_data %>%
      group_by(group) %>%
      mutate(
        group_mean = mean(`x1`),
        dmso_mean = mean(`x3`)
      ) %>%
      ungroup() %>%
      mutate(across(.cols = `x4`:`x11`, ~ (. - group_mean) / (dmso_mean - group_mean)))

    drug1_file <- df1 %>%
      # slice(1:3)
      slice_head(n = 3)

    drug2_file <- df1 %>%
      # slice(1:3)
      slice_tail(n = 3)

    cal_ic50_3p(drug_file = drug1_file,drug_name = t1,
      outpath = outpath, cellname = cellname
    )
    cal_ic50_3p(drug_file = drug2_file, drug_name = t2,
      outpath = outpath, cellname = cellname
    )
  }
}


ic50_pipe_3p_2blank <- function(index,
                                file = filepath1,
                                outpath = outpath1,
                                cellname = 'hPBMC'
                                ) {
  target <- all_sheets[[index]]
  print(target)

  if (length(str_split(target, pattern = " +")[[1]]) == 1) {
    message("Only one drug in this Sheet!!")

    t1 <- str_trim(target)

    all_raw_data <- readxl::read_excel(
      file,
      sheet = target
    ) %>%
      janitor::clean_names() %>%
      janitor::remove_empty(which = "rows")

    all_raw_data %<>% mutate(group = rep(t1, each = 3)) %>%
      dplyr::select(-`temperature_c`)

    df1 <- all_raw_data %>%
      group_by(group) %>%
      mutate(
        group_mean = mean(`x2`),
        dmso_mean = mean(`x3`)
      ) %>%
      ungroup() %>%
      mutate(across(.cols = `x4`:`x11`, ~ (. - group_mean) / (dmso_mean - group_mean)))

    drug1_file <- df1 %>%
      # slice(1:3)
      slice_head(n = 3)

    cal_ic50_3p(drug1_file, t1,
      outpath = outpath,cellname = cellname
    )
  } else {
    t1 <- str_split(target, pattern = " +")[[1]][1]
    t2 <- str_split(target, pattern = " +")[[1]][2]


    all_raw_data <- readxl::read_excel(
      file,
      sheet = target
    ) %>%
      janitor::clean_names() %>%
      janitor::remove_empty(which = "rows")

    all_raw_data %<>% mutate(group = rep(c(t1, t2), each = 3)) %>%
      dplyr::select(-`temperature_c`)

    df1 <- all_raw_data %>%
      group_by(group) %>%
      mutate(
        group_mean = mean(`x2`),
        dmso_mean = mean(`x3`)
      ) %>%
      ungroup() %>%
      mutate(across(.cols = `x4`:`x11`, ~ (. - group_mean) / (dmso_mean - group_mean)))

    drug1_file <- df1 %>%
      # slice(1:3)
      slice_head(n = 3)

    drug2_file <- df1 %>%
      # slice(1:3)
      slice_tail(n = 3)

    cal_ic50_3p(drug1_file, t1,
      outpath = outpath,cellname = cellname
    )
    cal_ic50_3p(drug2_file, t2,
      outpath = outpath,cellname = cellname
    )
  }
}

### start here
for (i in seq_along(all_sheets)) {
  print(i)

  ic50_pipe_3p_2blank(i)
  # ic50_pipe_3p_1blank(i)
}



# 汇总EC50数据: ---------------------------------------------------------------


all_files <- list.files(
  path = "~/OneDrive/kintor/Daily_Work/IC50_xiaodan/IC50_patch/wa_hk2_1111/",
  # pattern = "*_4P.csv",
  pattern = "*_3P.csv",
  full.names = TRUE
)

f1 <- read_csv(all_files[[1]]) %>%
  filter(rowname == "EC50")

for (i in 2:length(all_files)) {
  fn <- read_csv(all_files[[i]]) %>%
    filter(rowname == "EC50")
  f1 <- f1 %>% bind_rows(fn)
}


f1 %>%
  mutate(label = ifelse((`p.value` >= 0.05 | `2.5 %` < 0), 'Warning', 'Looks Good')) %>%
  write_excel_csv("~/OneDrive/kintor/Daily_Work/IC50_xiaodan/wa_hk2_1111_3p.csv")


## 找到异常的值
f1 %>% filter(`p.value` >= 0.05 | `2.5 %` < 0)









