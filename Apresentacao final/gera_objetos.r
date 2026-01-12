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
library(tidyr)


# Mapa ----

options(geobr.show_message = FALSE)
# Carregar pacotes necessários
library(ggplot2)
library(sf)
library(geobr)
library(dplyr)
library(ggrepel)  # Para melhor posicionamento dos labels

set.seed(123)

# 1. Obter dados dos biomas e remover "Zona Costeira"
biomas <- geobr::read_biomes(year = 2019) %>%
  filter(name_biome != "Sistema Costeiro")

# 2. Calcular centroides para posicionar os labels
centroides_biomas <- biomas %>%
  st_centroid() %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  as.data.frame()

# 3. Definir cores para cada bioma
cores_biomas <- c(
  "Amazônia" = "#1B9E77",
  "Caatinga" = "#D95F02",
  "Cerrado" = "#7570B3",
  "Mata Atlântica" = "#E7298A",
  "Pampa" = "#66A61E",
  "Pantanal" = "#E6AB02"
)

# 4. Criar o mapa com labels
graf1 = ggplot() +
  # Camada dos biomas
  geom_sf(data = biomas, aes(fill = name_biome), color = "white", size = 0.2, alpha = 0.8) +
  
  # Camada dos labels com repel para evitar sobreposição
  geom_label_repel(
    data = centroides_biomas,
    aes(x = lon, y = lat, label = name_biome),
    size = 3,
    fontface = "bold",
    box.padding = 0.5,
    min.segment.length = 0,
    segment.color = "gray40"
  ) +
  
  # Escala de cores
  scale_fill_manual(values = cores_biomas, name = "Biomas") +
  
  theme_void() +  # Tema mais limpo que theme_minimal()
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),  # Fundo transparente
    panel.background = element_rect(fill = "transparent", color = NA)  # Fundo transparente
  )

graf1

# ggsave("mapa_biomas_brasil.png", width = 5, height = 4, dpi = 300)



# 



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
names(base_mod)[names(base_mod) == "T_MEDIO_PERMAN_MINISTRO_P_MANDATO"] <- "T_MEDIO_PERMAN_MINISTRO"


df_ger = base_mod |> filter(POL_ABV %in% c('pnab', 'pse') )
df_pnab = base_mod |> filter(POL_ABV == 'pnab')
df_pse = base_mod |> filter(POL_ABV == 'pse')



# Curva de sobrevivência - exemplo PNAB ----

# Kaplan-Meier para PNAB
KMg_pnab <- survfit(data = df_pnab, Surv(time = t, event = ADOPT) ~ 1, conf.int = FALSE)
km_pnab <- data.frame(
  time = c(0, KMg_pnab$time),
  surv = c(1, KMg_pnab$surv),
  grupo = "PNAB"
)

# Criar dataframe com as censuras específicas que você quer marcar
censuras <- data.frame(
  time = max(df_pnab$t),
  grupo = "PNAB"
)

idx_pnab <- max(which(km_pnab$time <= max(df_pnab$t)))
surv_pnab <- km_pnab$surv[idx_pnab]

# Adicionar as sobrevivências ao dataframe
censuras$surv <- surv_pnab

ggplot(km_pnab, aes(x = time, y = surv, color = grupo)) +
  geom_step(size = 1) +
  # Censuras menores e mais discretas
  geom_point(data = censuras, shape = 3, size = 2, stroke = 0.8,
             show.legend = FALSE) +
  scale_color_manual(
    name = NULL, values = c("PNAB" = "steelblue")
  ) +
  scale_y_continuous(limits = c(0,1), labels = number_format(decimal.mark = ",", big.mark = ".")) +
  scale_x_continuous(labels = number_format(decimal.mark = ",", big.mark = ".")) +
  labs(
    title = "Curva de Sobrevivência - PNAB",
    x = "Tempo", 
    y = "Ŝ(t)"
  ) +
  theme_minimal() +
  theme(
    # Removida a legenda
    legend.position = "none",
    # Estilização do título
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 10)),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
  )

#ggsave("imagens/curva_km_pnab.png", width = 205, height = 121, units = "mm", bg = "white")



# Tempos - gráfico de dispersão ----

library(plotly)
library(htmlwidgets)

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

df <- df %>%
  left_join(df_pnab %>% select(MUNICÍPIO, SIGUF), 
            by = "MUNICÍPIO")

df = df |> mutate(label = paste0(MUNICÍPIO, ' - ', SIGUF,
                                 '\nPNAB: ', format(t_pol1, big.mark = '.', decimal.mark = ','),
                                 '\nPSE: ', format(t_pol2, big.mark = '.', decimal.mark = ','),
                                 '\nAdesão: ', adesao) )

# Certifique-se de ter os pacotes necessários
library(plotly)
library(dplyr)
library(scales)

# Gráfico
p <- df %>%
  plot_ly(
    x = ~t_pol1,
    y = ~t_pol2,
    color = ~adesao,
    colors = cores,
    type = "scatter",
    mode = "markers",
    marker = list(
      size = 8,
      opacity = 0.6,
      line = list(width = 0.5, color = "white")
    ),
    text = ~label,  # Apenas a label, sem X e Y
    hoverinfo = "text",  # Mostra apenas o texto personalizado
    hoverlabel = list(
      bgcolor = "white",
      font = list(size = 12, color = "black"),
      bordercolor = "gray"
    )
  ) %>%
  layout(
    title = list(
      text = "<b>Distribuição dos municípios por tempos até adoção das Políticas PNAB e PSE</b>",
      font = list(
        size = 20,
        family = "Arial, sans-serif",
        color = "black"
      ),
      x = 0.5,  # Centralizado
      xanchor = "center",
      y = 0.95,  # Posição vertical
      yanchor = "top"
    ),
    xaxis = list(
      title = "Tempos - PNAB",
      tickformat = ",.0f",  # Formato com separador de milhar
      tickformatstops = list(
        list(
          dtickrange = list(NULL, NULL),
          value = ".,0f"  # Força ponto como separador de milhar
        )
      ),
      separatethousands = TRUE,
      ticktext = ~sprintf("%.0f", t_pol1),  # Formata com ponto
      gridcolor = "lightgray",
      zerolinecolor = "lightgray"
    ),
    yaxis = list(
      title = "Tempos - PSE",
      tickformat = ",.0f",
      tickformatstops = list(
        list(
          dtickrange = list(NULL, NULL),
          value = ".,0f"  # Força ponto como separador de milhar
        )
      ),
      separatethousands = TRUE,
      ticktext = ~sprintf("%.0f", t_pol2),  # Formata com ponto
      gridcolor = "lightgray",
      zerolinecolor = "lightgray"
    ),
    hovermode = "closest",
    plot_bgcolor = "white",
    paper_bgcolor = "white",
    legend = list(
      title = list(text = "Adoção à política"),
      font = list(size = 11),
      orientation = "h",
      x = 0.5,
      y = -0.2,
      xanchor = "center"
    ),
    margin = list(l = 60, r = 40, t = 80, b = 80)  # Aumentei o top margin para acomodar o título
  )

# Mostrar o gráfico
p



# KM para ambas ----

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
    title = "Curvas de sobrevivência para a PNAB e PSE",
    x = "Tempo", y = "Ŝ(t)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "gray"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
  )
#ggsave("imagens/curvas_km_ambas.png", width = 205, height = 121, units = "mm", bg = "white")



# Ajuste paramétrico ----

## Ajuste PNAB ----

library(ggplot2)
library(survival)
library(survminer)

# Ajuste Kaplan-Meier
KMg <- survfit(data = df_pnab, Surv(time = t, event = ADOPT) ~ 1, conf.int = FALSE)

# Extrair dados do Kaplan-Meier
km_data <- data.frame(
  time = KMg$time,
  surv = KMg$surv
)

# Ajustes dos modelos paramétricos
ti <- df_pnab$t
deltai <- df_pnab$ADOPT

# Exponencial
MLE <- survreg(Surv(ti, deltai) ~ 1, dist = 'exponential')
lambda <- 1/exp(MLE$coefficients)

# Weibull
mwe <- survreg(Surv(ti, deltai) ~ 1, dist = 'weibull')
(alfa_weibull <- exp(mwe$coefficients[1]))
(gama_weibull <- 1/mwe$scale)

# Lognormal
mln <- survreg(Surv(ti, deltai) ~ 1, dist = 'lognormal')
mu <- mln$coefficients[1]
sigma <- mln$scale

# Loglogística
mll <- survreg(Surv(ti, deltai) ~ 1, dist = "loglogistic")
alfa_loglog <- exp(mll$coefficients[1])
gama_loglog <- 1/mll$scale

# Criar dados para as curvas teóricas
x_vals <- seq(0, max(ti), length.out = 1000)

curves_data <- data.frame(
  x = x_vals,
  Exponencial = exp(-lambda * x_vals),
  Weibull = exp(-(x_vals/alfa_weibull)^gama_weibull),
  Lognormal = pnorm((-log(x_vals) + mu)/sigma, mean = 0, sd = 1),
  Loglogística = 1/(1 + (x_vals/alfa_loglog)^gama_loglog)
)

km_data_completo <- data.frame(
  time = c(0, km_data$time),
  surv = c(1, km_data$surv)
)

ggplot() +
  # Curva Kaplan-Meier
  geom_step(data = km_data_completo, aes(x = time, y = surv, color = "Kaplan-Meier"), 
            size = .8) +
  # Curvas teóricas na ordem desejada
  geom_line(data = curves_data, aes(x = x, y = Exponencial, color = "Exponencial"), 
            size = .8) +
  geom_line(data = curves_data, aes(x = x, y = Weibull, color = "Weibull"), 
            size = .8) +
  geom_line(data = curves_data, aes(x = x, y = Loglogística, color = "Log-logística"), 
            size = .8) +
  geom_line(data = curves_data, aes(x = x, y = Lognormal, color = "Log-normal"), 
            size = .8) +
  # Escalas e labels
  scale_x_continuous(limits = c(0, max(ti)),
                     labels = number_format(decimal.mark = ",", big.mark = ".")) +
  scale_y_continuous(limits = c(0, 1),
                     labels = number_format(decimal.mark = ",", big.mark = ".")) +
  scale_color_manual(
    name = NULL,  # Remove o título
    breaks = c("Kaplan-Meier", "Exponencial", "Weibull", "Log-logística", "Log-normal"),
    values = c(
      "Kaplan-Meier" = "black",
      "Exponencial" = "#fcf458",
      "Weibull" = "#FB1808",
      "Log-logística" = "#3E9E0F",
      "Log-normal" = "steelblue"
    )
  ) +
  labs(
    title = 'Ajuste da distribuição de probabilidade para a PNAB',
    x = "Tempo",
    y = "Ŝ(t)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.85, 0.78),
    legend.background = element_rect(fill = "white", color = "gray"),
    legend.title = element_blank()
  )
#ggsave("imagens/ajuste_par_pnab.png", width = 205, height = 121, units = "mm", bg = "white")


## Ajuste PSE ----

library(ggplot2)
library(survival)
library(survminer)

# Ajuste Kaplan-Meier
KMg <- survfit(data = df_pse, Surv(time = t, event = ADOPT) ~ 1, conf.int = FALSE)

# Extrair dados do Kaplan-Meier
km_data <- data.frame(
  time = KMg$time,
  surv = KMg$surv
)

# Ajustes dos modelos paramétricos
ti <- df_pse$t
deltai <- df_pse$ADOPT

# Exponencial
MLE <- survreg(Surv(ti, deltai) ~ 1, dist = 'exponential')
lambda <- 1/exp(MLE$coefficients)

# Weibull
mwe <- survreg(Surv(ti, deltai) ~ 1, dist = 'weibull')
(alfa_weibull <- exp(mwe$coefficients[1]))
(gama_weibull <- 1/mwe$scale)

# Lognormal
mln <- survreg(Surv(ti, deltai) ~ 1, dist = 'lognormal')
mu <- mln$coefficients[1]
sigma <- mln$scale

# Loglogística
mll <- survreg(Surv(ti, deltai) ~ 1, dist = "loglogistic")
alfa_loglog <- exp(mll$coefficients[1])
gama_loglog <- 1/mll$scale

# Criar dados para as curvas teóricas
x_vals <- seq(0, max(ti), length.out = 1000)

curves_data <- data.frame(
  x = x_vals,
  Exponencial = exp(-lambda * x_vals),
  Weibull = exp(-(x_vals/alfa_weibull)^gama_weibull),
  Lognormal = pnorm((-log(x_vals) + mu)/sigma, mean = 0, sd = 1),
  Loglogística = 1/(1 + (x_vals/alfa_loglog)^gama_loglog)
)

km_data_completo <- data.frame(
  time = c(0, km_data$time),
  surv = c(1, km_data$surv)
)

ggplot() +
  # Curva Kaplan-Meier
  geom_step(data = km_data_completo, aes(x = time, y = surv, color = "Kaplan-Meier"), 
            size = .8) +
  # Curvas teóricas na ordem desejada
  geom_line(data = curves_data, aes(x = x, y = Exponencial, color = "Exponencial"), 
            size = .8) +
  geom_line(data = curves_data, aes(x = x, y = Weibull, color = "Weibull"), 
            size = .8) +
  geom_line(data = curves_data, aes(x = x, y = Loglogística, color = "Log-logística"), 
            size = .8) +
  geom_line(data = curves_data, aes(x = x, y = Lognormal, color = "Log-normal"), 
            size = .8) +
  # Escalas e labels
  scale_x_continuous(limits = c(0, max(ti)),
                     labels = number_format(decimal.mark = ",", big.mark = ".")) +
  scale_y_continuous(limits = c(0, 1),
                     labels = number_format(decimal.mark = ",", big.mark = ".")) +
  scale_color_manual(
    name = NULL,  # Remove o título
    breaks = c("Kaplan-Meier", "Exponencial", "Weibull", "Log-logística", "Log-normal"),
    values = c(
      "Kaplan-Meier" = "black",
      "Exponencial" = "#fcf458",
      "Weibull" = "#FB1808",
      "Log-logística" = "#3E9E0F",
      "Log-normal" = "steelblue"
    )
  ) +
  labs(
    title = 'Ajuste da distribuição de probabilidade para o PSE',
    x = "Tempo",
    y = "Ŝ(t)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.85, 0.78),
    legend.background = element_rect(fill = "white", color = "gray"),
    legend.title = element_blank()
  )
#ggsave("imagens/ajuste_par_pse.png", width = 205, height = 121, units = "mm", bg = "white")



# Covariáveis ----

library(patchwork)
library(ggplot2)
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


## Boxplot - ORC_TOTAL_EMPENHADO ----
# Avaliar valores e o ministério mencionado
ggplot(df_ger, aes(x = POL_ABV, y = ORC_TOTAL_EMPENHADO, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot do orçamento total executado",
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
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_orc_emp.png", width = 205, height = 121, units = "mm")
#ggsave("imagens/boxplot_orc_emp2.png", width = 158, height = 93, units = "mm", bg='white')


## TX_VETO ----
# Avaliar valores e o ministério mencionado
ggplot(df_ger, aes(x = POL_ABV, y = TX_VETO, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot da taxa de veto",
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
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_tx_veto.png", width = 158, height = 93, units = "mm", bg='white')



# Boxplots mod1 ----

## PIB ----
bp_pib = ggplot(df_ger, aes(x = POL_ABV, y = PIB, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot do PIB per capita",
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
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank()); bp_pib
#ggsave("imagens/boxplot_pib.png", width = 158, height = 93, units = "mm", bg='white')


## DISCR_DOTACAO ----
# Avaliar valores e o ministério mencionado
bp_discr_aut = ggplot(df_ger, aes(x = POL_ABV, y = DISCR_DOTACAO, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Orçamento discricionário autorizado",
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
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_discr_dotacao.pdf", width = 158, height = 93, units = "mm")


## DISCR_TX ----
# Avaliar valores e o ministério mencionado
bp_discr_tx = ggplot(df_ger, aes(x = POL_ABV, y = DISCR_TX, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Taxa de execução do orçamento discricionário",
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
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_discr_tx.pdf", width = 158, height = 93, units = "mm")


(bp_discr_aut + bp_discr_tx) + 
  plot_annotation(
    title = "Boxplots das variáveis explicativas",
    theme = theme(
      plot.title = element_text(hjust = 0.5),  # Reduzir tamanho se necessário
      plot.margin = margin(10, 0, 0, 10)  # cima, direita, baixo, esquerda
    )
  )
#ggsave("imagens/bp_mod1.png", width = 205, height = 121, units = "mm", bg = 'white')



# Boxplots mod1 ----

## PIB ----
bp_pib = ggplot(df_ger, aes(x = POL_ABV, y = PIB, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Boxplot do PIB per capita",
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
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank()); bp_pib
#ggsave("imagens/boxplot_pib.png", width = 158, height = 93, units = "mm", bg='white')


## DISCR_DOTACAO ----
# Avaliar valores e o ministério mencionado
bp_discr_aut = ggplot(df_ger, aes(x = POL_ABV, y = DISCR_DOTACAO, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Orçamento discricionário autorizado",
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
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_discr_dotacao.pdf", width = 158, height = 93, units = "mm")


## DISCR_TX ----
# Avaliar valores e o ministério mencionado
bp_discr_tx = ggplot(df_ger, aes(x = POL_ABV, y = DISCR_TX, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Taxa de execução do orçamento discricionário",
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
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank())
#ggsave("imagens/boxplot_discr_tx.pdf", width = 158, height = 93, units = "mm")



# Boxplots do modelo 4 ----

## SERVIDORES ----
# Avaliar valores e o ministério mencionado
bp_serv = ggplot(df_ger, aes(x = POL_ABV, y = SERVIDORES, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Número de servidores",
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
        plot.title = element_text(hjust = 0.5, size = 17),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank()); bp_serv
#ggsave("imagens/boxplot_servidores.pdf", width = 158, height = 93, units = "mm")


## TX_SUCESSO ----
# Avaliar valores e o ministério mencionado
bp_tx_sucesso = ggplot(df_ger, aes(x = POL_ABV, y = TX_SUCESSO, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Taxa de aprovação legislativa",
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
        plot.title = element_text(hjust = 0.5, size = 17),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank()); bp_tx_sucesso
#ggsave("imagens/boxplot_tx_sucesso.pdf", width = 158, height = 93, units = "mm")


## INDICE_MNST ----
# Avaliar valores e o ministério mencionado
bp_tx_mnst = ggplot(df_ger, aes(x = POL_ABV, y = INDICE_MNST, fill = as.factor(ADOPT))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black",
               position = position_dodge(width = 0.75) ) +
  facet_wrap(~ POL_ABV, scales = "free_x", strip.position = "bottom") +
  labs(title = "Taxa ministerial",
       x = "POL_ABRV", y = "Taxa ministerial", fill = 'Adoção à política') +
  scale_y_continuous(labels = label_number(
    scale = 100, prefix = "", suffix = "%",
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
        plot.title = element_text(hjust = 0.5, size = 17),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),  # Remove grade vertical principal
        panel.grid.minor.x = element_blank()) ; bp_tx_mnst
#ggsave("imagens/boxplot_indice_mnst.pdf", width = 158, height = 93, units = "mm")

(bp_serv + bp_tx_sucesso + bp_tx_mnst) + 
  plot_annotation(
    title = "Boxplots das variáveis explicativas",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 17),  # Reduzir tamanho se necessário
      plot.margin = margin(10, 0, 0, 10)  # cima, direita, baixo, esquerda
    )
  )
#ggsave("imagens/bp_mod4.png", width = 300, height = 121, units = "mm", bg = 'white')




# Resumo modelos ----

library(readxl)
library(dplyr)
avaliacao_modelos <- read_excel("avaliacao_modelos.xlsx", 
                                sheet = "Modelos bivariados", skip = 2)

avaliacao_modelos = avaliacao_modelos[avaliacao_modelos$`troca PNAB` != 'descarte', ]
avaliacao_modelos = avaliacao_modelos[avaliacao_modelos$`troca PSE` != 'descarte', ]

# Separar os dados
dados_na <- avaliacao_modelos |> filter(!is.na(label))
dados_validos <- avaliacao_modelos |> filter(is.na(label))

# Gráfico de dispersão interativo
disp_mod <- plot_ly() |>
  # Primeiro plotar os pontos com valores (gradiente normal)
  add_trace(
    data = dados_validos,
    x = ~BIC,
    y = ~theta,
    type = 'scatter',
    mode = 'markers',
    color = ~qtd_covs,
    colorscale = "Viridis",
    name = "Demais modelos",
    size = 0.8,
    text = ~paste0("BIC: ", BIC, "; phi: ", theta,
                   "\nPNAB: ", `Covariáveis PNAB`,
                   "\nPSE: ", `Covariáveis PSE`),
    hoverinfo = 'text'
  ) |>
  # Depois adicionar os pontos NA em vermelho por cima
  add_trace(
    data = dados_na,
    x = ~BIC,
    y = ~theta,
    type = 'scatter',
    mode = 'markers',
    marker = list(color = ~qtd_covs, colorscale = "Viridis", cmin = 0, cmax = 10),
    name = "Modelos destacados",
    showlegend = TRUE,
    size = 0.8,
    text = ~paste0(label, "\nBIC: ", BIC, "; phi: ", theta,
                   "\nPNAB: ", `Covariáveis PNAB`,
                   "\nPSE: ", `Covariáveis PSE`),
    hoverinfo = 'text'
  )

disp_mod <- layout(
  disp_mod,
  title = "Gráfico de dispersão dos modelos bivariados",
  xaxis = list(title = "BIC"),
  yaxis = list(title = "φ"),
  legend = list(
    orientation = "h",      # 'h' para horizontal (embaixo), 'v' para vertical
    x = 0.5,               # Posição horizontal: 0 (esquerda), 0.5 (centro), 1 (direita)
    y = -0.2,              # Posição vertical: negativo coloca abaixo do gráfico
    xanchor = "center",    # Ancora do ponto x: "left", "center", "right"
    yanchor = "top")       # Ancora do ponto y: "top", "middle", "bottom"
)

disp_mod




# Resíduos ----

## Modelo 3 - PNAB ----

variaveis = variaveis[variaveis != "NUM_ENTIDADES"]

df_pnab = base_mod |> filter(POL_ABV == 'pnab')
df_pnab$NUM_ENTIDADES = NULL

mod_pnab_teste = survreg(data = df_pnab, Surv(time = t, event = ADOPT) ~
                           scale(log(PIB))
                         + scale(log(ORC_TOTAL_EMPENHADO))
                         # scale(log(PIB)) + scale(log(ORC_TOTAL_EMPENHADO)) +
                         # scale(DISCR_TX) + scale(log(SERVIDORES))
                         ,
                         dist = 'weibull')
summary(mod_pnab_teste); c(AIC(mod_pnab_teste), AICc(mod_pnab_teste), BIC(mod_pnab_teste))
vif(mod_pnab_teste)
check_multicol_survreg(mod_pnab_teste)

library(ggplot2)
library(scales)
# Análise de resíduos
# 1. Resíduos de Cox-Snell
xitbeta <- mod_pnab_teste$linear.predictors
scale_weibull <- mod_pnab_teste$scale

# Calcular resíduos de Cox-Snell para Weibull
residuos_cs <- (df_pnab$t / exp(xitbeta))^(1/scale_weibull)

# 2. Curva de sobrevivência dos resíduos
surv_residuos <- Surv(time = residuos_cs, event = df_pnab$ADOPT)
res_fit <- survfit(surv_residuos ~ 1)

# 3. Preparar dados para ggplot
# Extrair dados do objeto survfit
km_residuos <- data.frame(
  time = res_fit$time,
  surv = res_fit$surv
)

# Criar dados para a curva de referência Exponencial(1)
x_ref <- seq(0, max(residuos_cs, na.rm = TRUE), length.out = 100)
ref_data <- data.frame(
  x = x_ref,
  Exponencial = exp(-x_ref)
)

# 4. Plot com ggplot2 usando a estética especificada
res_mod3_pnab = ggplot() +
  # Curva Kaplan-Meier dos resíduos
  geom_step(data = km_residuos, aes(x = time, y = surv, color = "Resíduos"), 
            size = .8) +
  # Curva de referência Exponencial(1)
  geom_line(data = ref_data, aes(x = x, y = Exponencial, color = "Exponencial Padrão"), 
            size = .8) +
  # Escalas e labels
  scale_x_continuous(
    limits = c(0, max(residuos_cs, na.rm = TRUE)),
    breaks = seq(0, max(residuos_cs, na.rm = TRUE), by = 0.5),
    labels = number_format(decimal.mark = ",", big.mark = ".")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = number_format(decimal.mark = ",", big.mark = ".")
  ) +
  scale_color_manual(
    name = NULL,  # Remove o título
    breaks = c("Resíduos", "Exponencial Padrão"),
    values = c(
      "Resíduos" = 'steelblue' , # "#FB1808",
      "Exponencial Padrão" = "black"
    )
  ) +
  labs( title = 'PNAB', 
    x = expression(paste(e[i], " - Resíduos de Cox-Snell")),
    y = expression(paste(Ŝ(e[i]),))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = c(0.85, 0.78),
    legend.background = element_rect(fill = "white", color = "gray"),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  ); res_mod3_pnab
#ggsave("imagens/residuos_pnab_mod3.png", width = 158, height = 93, units = "mm", bg = "white")



## Modelo 3 - PSE ----

df_pse = base_mod |> filter(POL_ABV == 'pse')
df_pse$NUM_ENTIDADES = NULL

mod_pse_teste = survreg(data = df_pse, Surv(time = t, event = ADOPT) ~
                          scale(log(PIB)) + scale(TX_VETO) +
                          scale(log(DISCR_DOTACAO))
                        # scale(log(PIB)) + scale(ORC_TOTAL_EMPENHADO) +
                        # scale(log(TX_SUCESSO)) + scale(log(DISCR_TX)) +
                        # scale(log(INDICE_MNST)) 
                        
                        ,
                        dist = 'weibull')
summary(mod_pse_teste); c(AIC(mod_pse_teste), AICc(mod_pse_teste), BIC(mod_pse_teste))
vif(mod_pse_teste)
check_multicol_survreg(mod_pse_teste)

# Análise de resíduos
# 1. Resíduos de Cox-Snell
xitbeta <- mod_pse_teste$linear.predictors
scale_weibull <- mod_pse_teste$scale

# Calcular resíduos de Cox-Snell para Weibull
residuos_cs <- (df_pse$t / exp(xitbeta))^(1/scale_weibull)

# 2. Curva de sobrevivência dos resíduos
surv_residuos <- Surv(time = residuos_cs, event = df_pse$ADOPT)
res_fit <- survfit(surv_residuos ~ 1)

# 3. Preparar dados para ggplot
# Extrair dados do objeto survfit
km_residuos <- data.frame(
  time = res_fit$time,
  surv = res_fit$surv
)

# Criar dados para a curva de referência Exponencial(1)
x_ref <- seq(0, max(residuos_cs, na.rm = TRUE), length.out = 100)
ref_data <- data.frame(
  x = x_ref,
  Exponencial = exp(-x_ref)
)

# 4. Plot com ggplot2 usando a estética especificada
res_mod3_pse = ggplot() +
  # Curva Kaplan-Meier dos resíduos
  geom_step(data = km_residuos, aes(x = time, y = surv, color = "Resíduos"), 
            size = .8) +
  # Curva de referência Exponencial(1)
  geom_line(data = ref_data, aes(x = x, y = Exponencial, color = "Exponencial Padrão"), 
            size = .8) +
  # Escalas e labels
  scale_x_continuous(
    limits = c(0, max(residuos_cs, na.rm = TRUE)),
    #breaks = seq(0, max(residuos_cs, na.rm = TRUE), by = 0.5),
    labels = number_format(decimal.mark = ",", big.mark = ".")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = number_format(decimal.mark = ",", big.mark = ".")
  ) +
  scale_color_manual(
    name = NULL,  # Remove o título
    breaks = c("Resíduos", "Exponencial Padrão"),
    values = c(
      "Resíduos" = 'steelblue' , # "#FB1808",
      "Exponencial Padrão" = "black"
    )
  ) +
  labs(title = 'PSE', 
    x = expression(paste(e[i], " - Resíduos de Cox-Snell")),
    y = expression(paste(Ŝ(e[i]),))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = c(0.85, 0.78),
    legend.background = element_rect(fill = "white", color = "gray"),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  ); res_mod3_pse
#ggsave("imagens/residuos_pse_mod3.png", width = 158, height = 93, units = "mm", bg = "white")

### Junção ----

# Com título geral
(res_mod3_pnab + res_mod3_pse ) + 
  plot_annotation(
    title = "Resíduos de Cox-Snell para as marginais do Modelo 3",
    theme = theme(
      plot.title = element_text(hjust = 0.5),  # Reduzir tamanho se necessário
      plot.margin = margin(10, 30, 0, 10)  # cima, direita, baixo, esquerda
    )
  )
#ggsave("imagens/res_mod3.png", width = 205, height = 121, units = "mm", bg = 'white')



## Modelo 4 - PNAB ----

mod_pnab_teste = survreg(data = df_pnab, Surv(time = t, event = ADOPT) ~
                           scale(log(PIB)) + scale(log(ORC_TOTAL_EMPENHADO)) +
                           scale(DISCR_TX) + scale(log(SERVIDORES))
                         ,
                         dist = 'weibull')
summary(mod_pnab_teste); c(AIC(mod_pnab_teste), AICc(mod_pnab_teste), BIC(mod_pnab_teste))
vif(mod_pnab_teste)
check_multicol_survreg(mod_pnab_teste)

library(ggplot2)
library(scales)
# Análise de resíduos
# 1. Resíduos de Cox-Snell
xitbeta <- mod_pnab_teste$linear.predictors
scale_weibull <- mod_pnab_teste$scale

# Calcular resíduos de Cox-Snell para Weibull
residuos_cs <- (df_pnab$t / exp(xitbeta))^(1/scale_weibull)

# 2. Curva de sobrevivência dos resíduos
surv_residuos <- Surv(time = residuos_cs, event = df_pnab$ADOPT)
res_fit <- survfit(surv_residuos ~ 1)

# 3. Preparar dados para ggplot
# Extrair dados do objeto survfit
km_residuos <- data.frame(
  time = res_fit$time,
  surv = res_fit$surv
)

# Criar dados para a curva de referência Exponencial(1)
x_ref <- seq(0, max(residuos_cs, na.rm = TRUE), length.out = 100)
ref_data <- data.frame(
  x = x_ref,
  Exponencial = exp(-x_ref)
)

# 4. Plot com ggplot2 usando a estética especificada
res_mod4_pnab = ggplot() +
  # Curva Kaplan-Meier dos resíduos
  geom_step(data = km_residuos, aes(x = time, y = surv, color = "Resíduos"), 
            size = .8) +
  # Curva de referência Exponencial(1)
  geom_line(data = ref_data, aes(x = x, y = Exponencial, color = "Exponencial Padrão"), 
            size = .8) +
  # Escalas e labels
  scale_x_continuous(
    limits = c(0, max(residuos_cs, na.rm = TRUE)),
    breaks = seq(0, max(residuos_cs, na.rm = TRUE), by = 1),
    labels = number_format(decimal.mark = ",", big.mark = ".")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = number_format(decimal.mark = ",", big.mark = ".")
  ) +
  scale_color_manual(
    name = NULL,  # Remove o título
    breaks = c("Resíduos", "Exponencial Padrão"),
    values = c(
      "Resíduos" = 'steelblue' , # "#FB1808",
      "Exponencial Padrão" = "black"
    )
  ) +
  labs( title = 'PNAB',
    x = expression(paste(e[i], " - Resíduos de Cox-Snell")),
    y = expression(paste(Ŝ(e[i]),))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = c(0.85, 0.78),
    legend.background = element_rect(fill = "white", color = "gray"),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  ); res_mod4_pnab
#ggsave("imagens/residuos_pnab_mod4.png", width = 158, height = 93, units = "mm", bg = "white")



## Modelo 4 - PSE ----

mod_pse_teste = survreg(data = df_pse, Surv(time = t, event = ADOPT) ~
                          scale(log(PIB)) + scale(ORC_TOTAL_EMPENHADO) +
                          scale(log(TX_SUCESSO)) + scale(log(DISCR_TX)) +
                          scale(log(INDICE_MNST))
                        
                        ,
                        dist = 'weibull')
summary(mod_pse_teste); c(AIC(mod_pse_teste), AICc(mod_pse_teste), BIC(mod_pse_teste))
vif(mod_pse_teste)
check_multicol_survreg(mod_pse_teste)

# Análise de resíduos
# 1. Resíduos de Cox-Snell
xitbeta <- mod_pse_teste$linear.predictors
scale_weibull <- mod_pse_teste$scale

# Calcular resíduos de Cox-Snell para Weibull
residuos_cs <- (df_pse$t / exp(xitbeta))^(1/scale_weibull)

# 2. Curva de sobrevivência dos resíduos
surv_residuos <- Surv(time = residuos_cs, event = df_pse$ADOPT)
res_fit <- survfit(surv_residuos ~ 1)

# 3. Preparar dados para ggplot
# Extrair dados do objeto survfit
km_residuos <- data.frame(
  time = res_fit$time,
  surv = res_fit$surv
)

# Criar dados para a curva de referência Exponencial(1)
x_ref <- seq(0, max(residuos_cs, na.rm = TRUE), length.out = 100)
ref_data <- data.frame(
  x = x_ref,
  Exponencial = exp(-x_ref)
)

# 4. Plot com ggplot2 usando a estética especificada
res_mod4_pse = ggplot() +
  # Curva Kaplan-Meier dos resíduos
  geom_step(data = km_residuos, aes(x = time, y = surv, color = "Resíduos"), 
            size = .8) +
  # Curva de referência Exponencial(1)
  geom_line(data = ref_data, aes(x = x, y = Exponencial, color = "Exponencial Padrão"), 
            size = .8) +
  # Escalas e labels
  scale_x_continuous(
    limits = c(0, max(residuos_cs, na.rm = TRUE)),
    breaks = seq(0, max(residuos_cs, na.rm = TRUE), by = 2),
    labels = number_format(decimal.mark = ",", big.mark = ".")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = number_format(decimal.mark = ",", big.mark = ".")
  ) +
  scale_color_manual(
    name = NULL,  # Remove o título
    breaks = c("Resíduos", "Exponencial Padrão"),
    values = c(
      "Resíduos" = 'steelblue' , # "#FB1808",
      "Exponencial Padrão" = "black"
    )
  ) +
  labs(title = 'PSE', 
       x = expression(paste(e[i], " - Resíduos de Cox-Snell")),
       y = expression(paste(Ŝ(e[i]),))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = c(0.85, 0.78),
    legend.background = element_rect(fill = "white", color = "gray"),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  ); res_mod4_pse
#ggsave("imagens/residuos_pse_mod4.png", width = 158, height = 93, units = "mm", bg = "white")


### Junção ----

# Com título geral
(res_mod4_pnab + res_mod4_pse ) + 
  plot_annotation(
    title = "Resíduos de Cox-Snell para as marginais do Modelo 4",
    theme = theme(
      plot.title = element_text(hjust = 0.5),  # Reduzir tamanho se necessário
      plot.margin = margin(10, 30, 0, 10)  # cima, direita, baixo, esquerda
    )
  )
#ggsave("imagens/res_mod4.png", width = 205, height = 121, units = "mm", bg = 'white')
