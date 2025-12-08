library(readr)
library(dplyr)
library(ggplot2)
library(survival)
library(AdequacyModel)
library(AICcmodavg)
library(zoo)
library(readxl)
library(plotly)
library(scales)

caminho = dirname(rstudioapi::getActiveDocumentContext()$path)
base_completa = read_excel(paste0(caminho, '/', "BASE_FINAL.xlsx"), 
                           sheet = "BASE_COMPLETA")
base_mod = read_excel(paste0(caminho, '/', "BASE_FINAL.xlsx"), 
                      sheet = "BASE_MOD") |>
  mutate(#TX_VETO = TX_VETO*100,
    TX_CONFLITO_MP = TX_CONFLITO_MP*100,
    TX_PARTC = TX_PARTC*100,
    TX_SUCESSO = TX_SUCESSO*100,
    TX_DOMINANCIA = TX_DOMINANCIA*100,
    TAM_PROP_BASE = TAM_PROP_BASE*100,
    TX_DISCIP_BASE = TX_DISCIP_BASE*100,
    TX_PROP_MNST = TX_PROP_MNST*100,
    POP_PRSD = POP_PRSD*100,
    DISCR_TX = DISCR_TX*100,
    ORC_TOTAL_TX = ORC_TOTAL_TX*100
  )

df_ger = base_mod |> filter(POL_ABV %in% c('pnab', 'pse') )
df_pnab = base_mod |> filter(POL_ABV == 'pnab')
df_pse = base_mod |> filter(POL_ABV == 'pse')


# Tempos - gráfico de dispersão ----

library(tidyverse)
library(reshape2)
library(lubridate)
library(dplyr)

cores <- c('Nenhuma' = '#C0C0C0', 'Ambas' = '#369a21', "Apenas PNAB" = "#2e27b0", "Apenas PSE" = '#FB1808')

df <- df_ger %>%
  mutate(
    POL_ABV = case_when(POL_ABV == "pnab" ~ "Apenas PNAB",
                        POL_ABV == "pse" ~ 'Apenas PSE',
                        T ~ POL_ABV)
  )

df_t = df |> dplyr::select(MUNICÍPIO, POL_ABV, t) |> pivot_wider(names_from = POL_ABV, values_from = t) |>
  rename(t_pol1 = 'Apenas PNAB', t_pol2 = 'Apenas PSE')
df_delta = df |> dplyr::select(MUNICÍPIO, POL_ABV, ADOPT) |> pivot_wider(names_from = POL_ABV, values_from = ADOPT) |>
  rename(delta_pol1 = 'Apenas PNAB', delta_pol2 = 'Apenas PSE')
df_ano = df |> dplyr::select(MUNICÍPIO, POL_ABV, ANO) |> pivot_wider(names_from = POL_ABV, values_from = ANO) |>
  rename(ano_pol1 = 'Apenas PNAB', ano_pol2 = 'Apenas PSE')

df = df_t |> left_join(df_delta) |> left_join(df_ano) # |> left_join(df_inicio) |> left_join(df_adopt)

df = df |> mutate(delta_pol1 = case_when(delta_pol1 == 0 ~ 'Não adotou',
                                         T ~ 'Adotou'),
                  delta_pol2 = case_when(delta_pol2 == 0 ~ 'Não adotou',
                                         T ~ 'Adotou')) |>
  mutate(
  adesao = case_when( (delta_pol1=='Adotou' & delta_pol2=='Adotou') ~ 'Ambas',
                      (delta_pol1=='Não adotou' & delta_pol2=='Não adotou') ~ 'Nenhuma',
                      (delta_pol1=='Adotou' & delta_pol2=='Não adotou') ~ 'Apenas PNAB',
                      TRUE ~ 'Apenas PSE')) |>
  mutate(adesao = factor(adesao, levels = c('Ambas', 'Apenas PNAB', 'Apenas PSE', 'Nenhuma')))

plot( df |> 
        ggplot(aes(x=t_pol1, y=t_pol2, color=adesao))+
        scale_color_manual(values = cores, name = "Adoção à política") +
        geom_point(size=2, alpha = 0.6) +
        scale_x_continuous(labels = label_number( big.mark = ".", decimal.mark = "," ) ) +
        scale_y_continuous(labels = label_number( big.mark = ".", decimal.mark = "," ) ) +
        labs(title = '', x = 'Tempos - PNAB', y = 'Tempos - PSE') +
        theme_minimal() )
#ggsave("imagens/graf_dispersao.pdf", width = 205, height = 121, units = "mm")

table(df$adesao); round(100*table(df$adesao)/sum(table(df$adesao)), 1)
cor(df$t_pol1, df$t_pol2)
cor(df|>filter(adesao=="Ambas")|>dplyr::select(t_pol1),
    df|>filter(adesao=="Ambas")|>dplyr::select(t_pol2) )



# Curva de sobrevivência - ambas ----

# Kaplan-Meier para PNAB
KMg_pnab <- survfit(data = df_pnab, Surv(time = t, event = ADOPT) ~ 1, conf.int = FALSE)
km_pnab <- data.frame(
  time = c(0, KMg_pnab$time),
  surv = c(1, KMg_pnab$surv),
  grupo = "PNAB"
)

# Kaplan-Meier para PSE
KMg_pse <- survfit(data = df_pse, Surv(time = t, event = ADOPT) ~ 1, conf.int = FALSE)
km_pse <- data.frame(
  time = c(0, KMg_pse$time),
  surv = c(1, KMg_pse$surv),
  grupo = "PSE"
)

# Combinar os dados
km_combined <- rbind(km_pnab, km_pse)

# Marcar censuras

# Criar dataframe com as censuras específicas que você quer marcar
censuras <- data.frame(
  time = c(max(df_pnab$t), max(df_pse$t)),
  grupo = c("PNAB", "PSE")
)

idx_pnab <- max(which(km_pnab$time <= max(df_pnab$t)))
surv_pnab <- km_pnab$surv[idx_pnab]

idx_pse <- max(which(km_pse$time <= max(df_pse$t)))
surv_pse <- km_pse$surv[idx_pse]

# Adicionar as sobrevivências ao dataframe
censuras$surv <- c(surv_pnab, surv_pse)

ggplot(km_combined, aes(x = time, y = surv, color = grupo)) +
  geom_step(size = 1) +
  # Censuras menores e mais discretas
  geom_point(data = censuras, shape = 3, size = 2, stroke = 0.8,
             show.legend = FALSE) +
  scale_color_manual(
    name = NULL, values = c("PNAB" = "steelblue", "PSE" = "#FB1808")
  ) +
  scale_y_continuous(limits = c(0,1), labels = number_format(decimal.mark = ",", big.mark = ".")) +
  scale_x_continuous(labels = number_format(decimal.mark = ",", big.mark = ".")) +
  labs(
    x = "Tempo", y = "Ŝ(t)"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "gray"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
  )
#ggsave("imagens/curvas_km.png", width = 205, height = 121, units = "mm", bg = "white")

table(df_pnab$ADOPT)
table(df_pnab$ADOPT)/sum(table(df_pnab$ADOPT))

table(df_pse$ADOPT)
table(df_pse$ADOPT)/sum(table(df_pse$ADOPT))



# Covariáveis ----

library(dplyr)
library(forcats)

df_ger = base_mod |> filter(POL_ABV %in% c('pnab', 'pse') )
df_ger <- df_ger %>%
  mutate(ADOPT = factor(ADOPT,
                        levels = c(1, 0),  # Valores originais
                        labels = c("Adotou", "Não adotou")) ) %>%
  mutate(POL_ABV = as.factor(POL_ABV) %>%
           fct_recode(
             "PNAB" = "pnab",
             "PSE" = "pse"
           ))

# Verificar os níveis do fator
levels(df_ger$ADOPT)


## PIB ----
ggplot(df_ger, aes(x = POL_ABV, y = PIB, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot do PIB por POL_ABRV e ADOPT",
       x = "POL_ABRV", y = "PIB per capita", fill = 'Adoção à política') +
  scale_y_continuous(labels = label_number(
    scale = 1/1000, prefix = "R$ ", suffix = " mil",
    big.mark = ".", decimal.mark = ","
  )) +
  scale_fill_manual(
    values = c("Adotou" = "#6A5ACD", "Não adotou" = "lightblue")
  ) +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_pib.pdf", width = 158, height = 93, units = "mm")

pnab1 = df_pnab|>filter(ADOPT==1)|>dplyr::select(PIB)|>unlist(); summary(pnab1); sd(pnab1)
pnab0 = df_pnab|>filter(ADOPT==0)|>dplyr::select(PIB)|>unlist(); summary(pnab0); sd(pnab0)
pse1 = df_pse|>filter(ADOPT==1)|>dplyr::select(PIB)|>unlist(); summary(pse1); sd(pse1)
pse0 = df_pse|>filter(ADOPT==0)|>dplyr::select(PIB)|>unlist(); summary(pse0); sd(pse0)


## ORC_TOTAL_EMPENHADO ----
# Avaliar valores e o ministério mencionado
ggplot(df_ger, aes(x = POL_ABV, y = ORC_TOTAL_EMPENHADO, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot do PIB por POL_ABRV e ADOPT",
       x = "POL_ABRV", y = "Orçamento total executado", fill = 'Adoção à política') +
  scale_y_continuous(labels = label_number(
    scale = 1/10^9, prefix = "R$ ", suffix = " bilhões",
    big.mark = ".", decimal.mark = ","
  )) +
  scale_fill_manual(
    values = c("Adotou" = "#6A5ACD", "Não adotou" = "lightblue")
  ) +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_orc_emp.pdf", width = 158, height = 93, units = "mm")

pnab1 = df_pnab|>filter(ADOPT==1)|>dplyr::select(ORC_TOTAL_EMPENHADO)|>unlist(); summary(pnab1)/(10^9); sd(pnab1)/(10^9)
pnab0 = df_pnab|>filter(ADOPT==0)|>dplyr::select(ORC_TOTAL_EMPENHADO)|>unlist(); summary(pnab0)/(10^9); sd(pnab0)/(10^9)
pse1 = df_pse|>filter(ADOPT==1)|>dplyr::select(ORC_TOTAL_EMPENHADO)|>unlist(); summary(pse1)/(10^9); sd(pse1)/(10^9)
pse0 = df_pse|>filter(ADOPT==0)|>dplyr::select(ORC_TOTAL_EMPENHADO)|>unlist(); summary(pse0)/(10^9); sd(pse0)/(10^9)



## TX_VETO ----
# Avaliar valores e o ministério mencionado
ggplot(df_ger, aes(x = POL_ABV, y = TX_VETO, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot do PIB por POL_ABRV e ADOPT",
       x = "POL_ABRV", y = "Taxa de veto", fill = 'Adoção à política') +
  scale_y_continuous(labels = label_number(
    scale = 1, prefix = "", suffix = "%",
    big.mark = ".", decimal.mark = ","
  )) +
  scale_fill_manual(
    values = c("Adotou" = "#6A5ACD", "Não adotou" = "lightblue")
  ) +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_tx_veto.pdf", width = 158, height = 93, units = "mm")

pnab1 = df_pnab|>filter(ADOPT==1)|>dplyr::select(TX_VETO)|>unlist(); summary(pnab1); sd(pnab1)
pnab0 = df_pnab|>filter(ADOPT==0)|>dplyr::select(TX_VETO)|>unlist(); summary(pnab0); sd(pnab0)
pse1 = df_pse|>filter(ADOPT==1)|>dplyr::select(TX_VETO)|>unlist(); summary(pse1); sd(pse1)
pse0 = df_pse|>filter(ADOPT==0)|>dplyr::select(TX_VETO)|>unlist(); summary(pse0); sd(pse0)



## DISCR_TX ----
# Avaliar valores e o ministério mencionado
ggplot(df_ger, aes(x = POL_ABV, y = DISCR_TX, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot do PIB por POL_ABRV e ADOPT",
       x = "POL_ABRV", y = "Taxa de execução do orçamento discricionário", fill = 'Adoção à política') +
  scale_y_continuous(labels = label_number(
    scale = 1, prefix = "", suffix = "%",
    big.mark = ".", decimal.mark = ","
  )) +
  scale_fill_manual(
    values = c("Adotou" = "#6A5ACD", "Não adotou" = "lightblue")
  ) +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_discr_tx.pdf", width = 158, height = 93, units = "mm")

pnab1 = df_pnab|>filter(ADOPT==1)|>dplyr::select(DISCR_TX)|>unlist(); summary(pnab1); sd(pnab1)
pnab0 = df_pnab|>filter(ADOPT==0)|>dplyr::select(DISCR_TX)|>unlist(); summary(pnab0); sd(pnab0)
pse1 = df_pse|>filter(ADOPT==1)|>dplyr::select(DISCR_TX)|>unlist(); summary(pse1); sd(pse1)
pse0 = df_pse|>filter(ADOPT==0)|>dplyr::select(DISCR_TX)|>unlist(); summary(pse0); sd(pse0)



## DISCR_DOTACAO ----
# Avaliar valores e o ministério mencionado
ggplot(df_ger, aes(x = POL_ABV, y = DISCR_DOTACAO, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot do PIB por POL_ABRV e ADOPT",
       x = "POL_ABRV", y = "Orçamento discricionário autorizado", fill = 'Adoção à política') +
  scale_y_continuous(labels = label_number(
    scale = 1/10^9, prefix = "R$ ", suffix = " bilhões",
    big.mark = ".", decimal.mark = ","
  )) +
  scale_fill_manual(
    values = c("Adotou" = "#6A5ACD", "Não adotou" = "lightblue")
  ) +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_discr_dotacao.pdf", width = 158, height = 93, units = "mm")

pnab1 = df_pnab|>filter(ADOPT==1)|>dplyr::select(DISCR_DOTACAO)|>unlist(); summary(pnab1)/(10^9); sd(pnab1)/(10^9)
pnab0 = df_pnab|>filter(ADOPT==0)|>dplyr::select(DISCR_DOTACAO)|>unlist(); summary(pnab0)/(10^9); sd(pnab0)/(10^9)
pse1 = df_pse|>filter(ADOPT==1)|>dplyr::select(DISCR_DOTACAO)|>unlist(); summary(pse1)/(10^9); sd(pse1)/(10^9)
pse0 = df_pse|>filter(ADOPT==0)|>dplyr::select(DISCR_DOTACAO)|>unlist(); summary(pse0)/(10^9); sd(pse0)/(10^9)



## SERVIDORES ----
# Avaliar valores e o ministério mencionado
ggplot(df_ger, aes(x = POL_ABV, y = SERVIDORES, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot do PIB por POL_ABRV e ADOPT",
       x = "POL_ABRV", y = "Número de servidores", fill = 'Adoção à política') +
  scale_y_continuous(labels = label_number(
    scale = 1, prefix = "", suffix = "",
    big.mark = ".", decimal.mark = ","
  )) +
  scale_fill_manual(
    values = c("Adotou" = "#6A5ACD", "Não adotou" = "lightblue")
  ) +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_servidores.pdf", width = 158, height = 93, units = "mm")

pnab1 = df_pnab|>filter(ADOPT==1)|>dplyr::select(SERVIDORES)|>unlist(); summary(pnab1); sd(pnab1)
pnab0 = df_pnab|>filter(ADOPT==0)|>dplyr::select(SERVIDORES)|>unlist(); summary(pnab0); sd(pnab0)
pse1 = df_pse|>filter(ADOPT==1)|>dplyr::select(SERVIDORES)|>unlist(); summary(pse1); sd(pse1)
pse0 = df_pse|>filter(ADOPT==0)|>dplyr::select(SERVIDORES)|>unlist(); summary(pse0); sd(pse0)

table(df_pnab$ANO, df_pnab$SERVIDORES)



## TX_SUCESSO ----
# Avaliar valores e o ministério mencionado
ggplot(df_ger, aes(x = POL_ABV, y = TX_SUCESSO, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot do PIB por POL_ABRV e ADOPT",
       x = "POL_ABRV", y = "Taxa de aprovação legislativa", fill = 'Adoção à política') +
  scale_y_continuous(labels = label_number(
    scale = 1, prefix = "", suffix = "%",
    big.mark = ".", decimal.mark = ","
  )) +
  scale_fill_manual(
    values = c("Adotou" = "#6A5ACD", "Não adotou" = "lightblue")
  ) +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_tx_sucesso.pdf", width = 158, height = 93, units = "mm")

pnab1 = df_pnab|>filter(ADOPT==1)|>dplyr::select(TX_SUCESSO)|>unlist(); summary(pnab1); sd(pnab1)
pnab0 = df_pnab|>filter(ADOPT==0)|>dplyr::select(TX_SUCESSO)|>unlist(); summary(pnab0); sd(pnab0)
pse1 = df_pse|>filter(ADOPT==1)|>dplyr::select(TX_SUCESSO)|>unlist(); summary(pse1); sd(pse1)
pse0 = df_pse|>filter(ADOPT==0)|>dplyr::select(TX_SUCESSO)|>unlist(); summary(pse0); sd(pse0)


## INDICE_MNST ----
# Avaliar valores e o ministério mencionado
ggplot(df_ger, aes(x = POL_ABV, y = INDICE_MNST, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot do PIB por POL_ABRV e ADOPT",
       x = "POL_ABRV", y = "Taxa ministerial", fill = 'Adoção à política') +
  scale_y_continuous(labels = label_number(
    scale = 1, prefix = "", suffix = "",
    big.mark = ".", decimal.mark = ","
  )) +
  scale_fill_manual(
    values = c("Adotou" = "#6A5ACD", "Não adotou" = "lightblue")
  ) +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_indice_mnst.pdf", width = 158, height = 93, units = "mm")

pnab1 = df_pnab|>filter(ADOPT==1)|>dplyr::select(INDICE_MNST)|>unlist(); summary(pnab1); sd(pnab1)
pnab0 = df_pnab|>filter(ADOPT==0)|>dplyr::select(INDICE_MNST)|>unlist(); summary(pnab0); sd(pnab0)
pse1 = df_pse|>filter(ADOPT==1)|>dplyr::select(INDICE_MNST)|>unlist(); summary(pse1); sd(pse1)
pse0 = df_pse|>filter(ADOPT==0)|>dplyr::select(INDICE_MNST)|>unlist(); summary(pse0); sd(pse0)

