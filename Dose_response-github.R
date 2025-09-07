library(ggplot2)
library(dplyr)

# --- Read data (Excel with file.choose) ---
 library(readxl)
 df <- read_excel(file.choose())

 ### Imatinib only ####
 
 # Prepare summary for Imatinib
 df_summary <- df %>%
   filter(Drug == "Imatinib") %>%
   mutate(Dose = ifelse(Dose == "0", "0 µM", Dose)) %>%  # rename 0
   group_by(Cell_Type, Dose) %>%
   summarise(Mean = mean(Cell_Count),
             SEM = sd(Cell_Count)/sqrt(n()),
             .groups = "drop")
 
 # Set order of doses
 df_summary$Dose <- factor(df_summary$Dose,
                           levels = c("0 µM", "10µM", "50µM", "100µM"))
 
 # Plot
 ggplot(df_summary, aes(x = Dose, y = Mean, fill = Cell_Type)) +
   geom_col(position = position_dodge(width = 0.9)) +
   geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                 position = position_dodge(width = 0.9), width = 0.3) +
   scale_y_continuous(labels = scales::scientific) +
   labs(title = "Imatinib Dose Response (Day 7)",
        x = "Imatinib Dose",
        y = "Mean Cell Count ± SEM") +
   theme_bw(base_size = 14) +
   theme(legend.position = "bottom")
#------------------
library(rstatix) # for t_test and significance annotation
library(ggpubr)
 
 df_summary_imat <- df %>%
   filter(Drug == "Imatinib") %>%
   mutate(Dose = ifelse(Dose == "0", "0 µM", Dose)) %>%  # rename 0
   group_by(Cell_Type, Dose) %>%
   summarise(Mean = mean(Cell_Count),
             SEM = sd(Cell_Count)/sqrt(n()),
             .groups = "drop")
 
 # Subset for Imatinib
 df_imat <- df %>% filter(Drug == "Imatinib")
 
 # Perform t-tests at each dose
 stat_test_imat <- df_imat %>%
   group_by(Dose) %>%
   t_test(Cell_Count ~ Cell_Type) %>%
   adjust_pvalue(method = "BH") %>%
   add_significance("p.adj")
 
 stat_test_imat 
 
 #----alternative idea without ggpubr --------
 
 ggplot(df_summary_imat, aes(x = Dose, y = Mean, fill = Cell_Type)) +
   geom_col(position = position_dodge(width = 0.9)) +
   geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                 position = position_dodge(width = 0.9), width = 0.3) +
   scale_y_continuous(labels = scales::scientific) +
   labs(title = "Imatinib Dose Response (Day 7)",
        x = "Imatinib Dose",
        y = "Mean Cell Count ± SEM") +
   theme_bw(base_size = 14) +
   theme(legend.position = "bottom") +
   # Add significance stars manually
   annotate("text", x = 1, y = 2.3e5, label = "*") +
   annotate("text", x = 2, y = 2.5e5, label = "*") +
   annotate("text", x = 3, y = 1.8e5, label = "**") 
 
 
 
 ## Dasatinib only ## 
 
 # Prepare summary for Dasatinib
 df_summary_dasa <- df %>%
   filter(Drug == "Dasatinib") %>%
   mutate(Dose = ifelse(Dose == "0", "0 nM", Dose)) %>%  # rename 0
   group_by(Cell_Type, Dose) %>%
   summarise(Mean = mean(Cell_Count),
             SEM = sd(Cell_Count)/sqrt(n()),
             .groups = "drop")
 
 # Set order of doses
 df_summary_dasa$Dose <- factor(df_summary_dasa$Dose,
                                levels = c("0 nM", "10nM", "50nM", "100nM"))
 
 # Plot
 ggplot(df_summary_dasa, aes(x = Dose, y = Mean, fill = Cell_Type)) +
   geom_col(position = position_dodge(width = 0.9)) +
   geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                 position = position_dodge(width = 0.9), width = 0.3) +
   scale_y_continuous(labels = scales::scientific) +
   labs(title = "Dasatinib Dose Response (Day 7)",
        x = "Dasatinib Dose",
        y = "Mean Cell Count ± SEM") +
   theme_bw(base_size = 14) +
   theme(legend.position = "bottom")

 #----------
 library(rstatix) # for t_test and significance annotation
 library(ggpubr)
 
 
 # Subset for Dasatinib
 df_dasa <- df %>% filter(Drug == "Dasatinib")
 
 # Perform t-tests at each dose
 stat_test_dasa <- df_dasa %>%
   group_by(Dose) %>%
   t_test(Cell_Count ~ Cell_Type) %>%
   adjust_pvalue(method = "BH") %>%
   add_significance("p.adj")
 
 stat_test_dasa 
 
 
 # Update the y-position for stars (a bit above the max bar for each dose)
 stat_test_dasa <- stat_test_dasa %>%
   mutate(y.position = c(1.35e7, 1.3e7, 1.25e7))  # adjust if needed
 
 ggplot(df_summary_dasa, aes(x = Dose, y = Mean, fill = Cell_Type)) +
   geom_col(position = position_dodge(width = 0.9)) +
   geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                 position = position_dodge(width = 0.9), width = 0.3) +
   scale_y_continuous(labels = scales::scientific) +
   labs(title = "Dasatinib Dose Response (Day 7)",
        x = "Dasatinib Dose",
        y = "Mean Cell Count ± SEM") +
   theme_bw(base_size = 14) +
   theme(legend.position = "bottom") +
   stat_pvalue_manual(stat_test_dasa, label = "p.adj.signif")

 
 #----alternative idea without ggpubr --------
 
 ggplot(df_summary_dasa, aes(x = Dose, y = Mean, fill = Cell_Type)) +
   geom_col(position = position_dodge(width = 0.9)) +
   geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                 position = position_dodge(width = 0.9), width = 0.3) +
   scale_y_continuous(labels = scales::scientific) +
   labs(title = "Dasatinib Dose Response (Day 7)",
        x = "Dasatinib Dose",
        y = "Mean Cell Count ± SEM") +
   theme_bw(base_size = 14) +
   theme(legend.position = "bottom") +
   # Add significance stars manually
   annotate("text", x = 1, y = 1.3e7, label = "**") +
   annotate("text", x = 2, y = 1.25e7, label = "**") +
   annotate("text", x = 3, y = 1.2e7, label = "*") 
 