library(readr)
library(dplyr)
library(ggplot2)
library(survival)
library(AdequacyModel)
library(AICcmodavg)
library(zoo)
library(readxl)
library(car)
library(writexl)
library(numDeriv)
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
names(base_mod)[names(base_mod) == "T_MEDIO_PERMAN_MINISTRO_P_MANDATO"] <- "T_MEDIO_PERMAN_MINISTRO"

variaveis <- c("TX_VETO", "TX_CONFLITO_MP", "TX_PARTC", "TX_SUCESSO", "TX_DOMINANCIA", 
               "TX_FRAG_PARTD_NEP", "TAM_PROP_BASE", "TX_DISCIP_BASE", "TX_PROP_MNST", 
               "POP_PRSD", "T_MEDIO_PERMAN_MINISTRO", "CARREIRAS", 
               "NUM_ENTIDADES", "SERVIDORES", "SERVIDORES_RJU", "DISCR_DOTACAO", 
               "DISCR_EMPENHADO", "ORC_TOTAL_DOTACAO", "ORC_TOTAL_EMPENHADO", 
               "DISCR_TX", "ORC_TOTAL_TX", "TI_PROJ_LEI", "TI_EMPENHADO", 
               "INDICE_MNST", "PIB", "ANO_ELEITORAL_MUNIC", "ANO_ELEITORAL_NAC", 
               "MRG_PRSD", "IDEOL_PRSD", "IDEOL_PRF")

df_ger = base_mod |> filter(POL_ABV %in% c('pnab', 'pse') )
df_pnab = base_mod |> filter(POL_ABV == 'pnab')
df_pse = base_mod |> filter(POL_ABV == 'pse')

check_multicol_survreg <- function(model) {
  # Extrair matriz do modelo sem intercepto
  X <- model.matrix(model)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, -1, drop = FALSE]
  }
  # Calcular matriz de correlação
  cor_matrix <- cor(X)
  # 1. Verificar correlações altas entre pares
  high_cor <- which(abs(cor_matrix) > 0.7 & upper.tri(cor_matrix), arr.ind = TRUE)
  if (length(high_cor) > 0) {
    cat("Correlações altas (> 0.7):\n")
    for (i in 1:nrow(high_cor)) {
      row_idx <- high_cor[i, 1]
      col_idx <- high_cor[i, 2]
      cat(sprintf("  %s - %s: %.3f\n", 
                  colnames(X)[row_idx], colnames(X)[col_idx], 
                  cor_matrix[row_idx, col_idx]))
    }}
  # 2. Calcular VIF manualmente
  if (ncol(X) >= 2) {
    vif_values <- diag(solve(cor_matrix))
    names(vif_values) <- colnames(X)
    return(vif_values)
  }
  return(NULL)
}


calc_gradiente <- function(para) {
  numDeriv::grad(vero_cop_debug, para)
}



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
    x = "Tempo",
    y = "Ŝ(t)"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.78),
    legend.background = element_rect(fill = "white", color = "gray"),
    legend.title = element_blank()
  )
#ggsave("imagens/ajuste_par_pnab.png", width = 205, height = 121, units = "mm", bg = "white")

# Comparação
AIC(MLE); AIC(mwe); AIC(mln); AIC(mll)
AICc(MLE); AICc(mwe); AICc(mln); AICc(mll)
BIC(MLE); BIC(mwe); BIC(mln); BIC(mll)

(trv <- 2*(mwe$loglik[1] - MLE$loglik[1]))
1 - pchisq(trv, 1) # p-valor



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
    x = "Tempo",
    y = "Ŝ(t)"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.78),
    legend.background = element_rect(fill = "white", color = "gray"),
    legend.title = element_blank()
  )
#ggsave("imagens/ajuste_par_pse.png", width = 205, height = 121, units = "mm", bg = "white")

# Comparação
AIC(MLE); AIC(mwe); AIC(mln); AIC(mll)
AICc(MLE); AICc(mwe); AICc(mln); AICc(mll)
BIC(MLE); BIC(mwe); BIC(mln); BIC(mll)

(trv <- 2*(mwe$loglik[1] - MLE$loglik[1]))
1 - pchisq(trv, 1) # p-valor




# Gráfico de dispersão de modelos (avaliação de modelos) ----

library(readxl)
avaliacao_modelos <- read_excel("avaliacao_modelos.xlsx", 
                                sheet = "Modelos bivariados", skip = 2)

avaliacao_modelos = avaliacao_modelos[avaliacao_modelos$`troca PNAB` != 'descarte', ]
avaliacao_modelos = avaliacao_modelos[avaliacao_modelos$`troca PSE` != 'descarte', ]

avaliacao_modelos_2 = avaliacao_modelos |> filter(theta>0)

cor(avaliacao_modelos$BIC, avaliacao_modelos$theta)
cor(avaliacao_modelos_2$BIC, avaliacao_modelos_2$theta)

library(plotly)
library(dplyr)
# Gráfico de dispersão interativo
grafico <- plot_ly(
  data = avaliacao_modelos,
  x = ~BIC,
  y = ~theta,
  type = 'scatter',
  mode = 'markers',
  color = ~qtd_covs,
  # size = ~tamanho,
  text = ~paste0("BIC: ", BIC, "; theta: ", theta, "; tau: ", round(tau, 4),
                 "; qtd_covs: ", qtd_covs, "\nPNAB (id ", `ID PNAB`,
                 "); BIC: ", `BIC PNAB`, "; Troca sinal: ", `troca PNAB`, "\n",
                 `Covariáveis PNAB`, "\nPSE (id ", `ID PSE`, "); BIC: ", `BIC PSE`,
                 "; Troca sinal: ", `troca PSE`, "\n", `Covariáveis PSE`),
  hoverinfo = 'text'
) %>%
  layout(
    title = "Gráfico de Dispersão Bivariado Interativo",
    xaxis = list(title = "BIC"),
    yaxis = list(title = "theta (alpha)")
  ); grafico



# install.packages(c("ggplot2", "plotly"))

library(ggplot2)
library(plotly)
library(ggrepel)

gg <- ggplot(avaliacao_modelos, aes(x = BIC, y = theta, 
                                    text = paste0("BIC: ", BIC, "  theta: ", theta, 
                                                  "\nPNAB:\n", `Covariáveis PNAB`, 
                                                  "\nPSE:\n", `Covariáveis PSE`),
                                    color = qtd_covs)) +  # Cor baseada em qtd_covs
  geom_point(alpha = 0.7, size=2) +
  geom_label_repel(aes(x = BIC, y = theta, label = label),
                   box.padding = 1, point.padding = 0.3,
                   segment.color = 'grey50',
                   min.segment.length = 0, force = 10, seed = 1) +
  # Linha horizontal em y = 0
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) +
  scale_x_continuous(labels = number_format(decimal.mark = ",", big.mark = ".")) +
  scale_y_continuous(limits = c(-0.18, 1), breaks = seq(-0.2, 1, 0.2),
                     labels = number_format(decimal.mark = ",", big.mark = ".")) +
  scale_color_viridis_c(name = "Número de\ncovariáveis",
                        breaks = c(0, 2, 4, 6, 8, 10)) +  # Escala de cores contínua
  labs(x = "BIC", y = expression(phi)) +
  theme_minimal(); gg

#ggsave("imagens/dispersao_modelos.pdf", width = 205, height = 121, units = "mm")

# Converter para interativo
ggplotly(gg, tooltip = 'text')



# Superfícies ----

## funções de cópulas - Breno ----

t1 = df_pnab$t; delta1 = df_pnab$ADOPT
t2 = df_pse$t; delta2 = df_pse$ADOPT

num_grande = 1e+300

S_theta = function(t, alfa, gama, theta){
  res = exp( (t/alfa)^gama *theta )
  res = min(num_grande, res)
  return(res)
}
S_theta = function(t, alfa, gama, theta){
  exp( (t/alfa)^gama *theta )
}


log_S_theta = function(t, alfa, gama, theta){
  (t/alfa)^gama *theta
}

num_pequeno = 1e-300

S_cop = function(alfa1, gama1, alfa2, gama2, theta){
  res = case_when( (delta1 == 0 & delta2 == 0) ~ 
                     (S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1) ^(-1/theta),
                   T ~ 1 )
  res = ifelse(res == 0, num_pequeno, res)
  return(res)
}

log_S_cop = function(alfa1, gama1, alfa2, gama2, theta){
  res = case_when( (delta1 == 0 & delta2 == 0) ~ 
                     (-1/theta) * log(S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1),
                   T ~ 0 )
  return(res)
}


parc_t1 = function(alfa1, gama1, alfa2, gama2, theta){
  res = case_when( (delta1 == 1 & delta2 == 0) ~
                     gama1/alfa1 *(t1/alfa1)^(gama1-1) * S_theta(t1, alfa1, gama1, theta) *
                     (S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1) ^((-1/theta)-1),
                   T ~ 1 )
  res = ifelse(res == 0, num_pequeno, res)
  return(res)
}

log_parc_t1 = function(alfa1, gama1, alfa2, gama2, theta){
  res = case_when( (delta1 == 1 & delta2 == 0) ~
                     log(gama1) -log(alfa1) +(gama1-1)*log(t1/alfa1) +log_S_theta(t1, alfa1, gama1, theta) +
                     ((-1/theta)-1)*log(S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1),
                   T ~ 0 )
  return(res)
}


parc_t2 = function(alfa1, gama1, alfa2, gama2, theta){
  res = case_when( (delta1 == 0 & delta2 == 1) ~
                     gama2/alfa2 *(t2/alfa2)^(gama2-1) * S_theta(t2, alfa2, gama2, theta) *
                     (S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1) ^((-1/theta)-1),
                   T ~ 1 )
  res = ifelse(res == 0, num_pequeno, res)
}

log_parc_t2 = function(alfa1, gama1, alfa2, gama2, theta){
  res = case_when( (delta1 == 0 & delta2 == 1) ~
                     log(gama2) -log(alfa2) +(gama2-1)*log(t2/alfa2) +log_S_theta(t2, alfa2, gama2, theta) +
                     ((-1/theta)-1)*log(S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1),
                   T ~ 0 )
}


f_cop = function(alfa1, gama1, alfa2, gama2, theta){
  res = case_when( (delta1 == 1 & delta2 == 1) ~
                     (1+theta) *gama1/alfa1*(t1/alfa1)^(gama1-1) *gama2/alfa2*(t2/alfa2)^(gama2-1) *
                     S_theta(t1, alfa1, gama1, theta) * S_theta(t2, alfa2, gama2, theta) * 
                     (S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1) ^((-1/theta)-2),
                   T ~ 1 )
  res = ifelse(res == 0, num_pequeno, res)
  return(res)
}

log_f_cop = function(alfa1, gama1, alfa2, gama2, theta){
  res = case_when( (delta1 == 1 & delta2 == 1) ~
                     log(1+theta) +log(gama1) -log(alfa1) +(gama1-1)*log(t1/alfa1) +log(gama2) -log(alfa2) +(gama2-1)*log(t2/alfa2) +
                     log_S_theta(t1, alfa1, gama1, theta) +log_S_theta(t2, alfa2, gama2, theta) + 
                     ((-1/theta)-2)*log(S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1),
                   T ~ 0 )
  return(res)
}


## Estimativa pontual ----

S_estim = function(t, cov, b0, b1, gama){
  alfa = exp(b0 + b1*cov)
  exp( -(t/alfa)^gama )
}

S_cop_estim = function(t1, t2, cov1, cov2, b01, b11, gama1, b02, b12, gama2, theta){
  alfa1 = exp(b01 + b11*cov1); alfa2 = exp(b02 + b12*cov2)
  (S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1) ^(-1/theta)
}

f_cop_estim = function(t1, t2, cov1, cov2, b01, b11, gama1, b02, b12, gama2, theta){
  alfa1 = exp(b01 + b11*cov1); alfa2 = exp(b02 + b12*cov2)
  
  (1+theta) *gama1/alfa1*(t1/alfa1)^(gama1-1) *gama2/alfa2*(t2/alfa2)^(gama2-1) *
    S_theta(t1, alfa1, gama1, theta) * S_theta(t2, alfa2, gama2, theta) *
    (S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1) ^((-1/theta)-2)
}

h_cop_estim = function(t1, t2, cov1, cov2, b01, b11, gama1, b02, b12, gama2, theta){
  f_cop_estim(t1, t2, cov1, cov2, b01, b11, gama1, b02, b12, gama2, theta) /
    S_cop_estim(t1, t2, cov1, cov2, b01, b11, gama1, b02, b12, gama2, theta)
}


## Modelo 1 ----

ajust1 = survreg(data = df_pnab, Surv(time = t, event = ADOPT) ~
                   scale(log(PIB)) + scale(log(ORC_TOTAL_EMPENHADO)) 
                 + scale(DISCR_TX)
                 ,
                 dist = 'weibull')
summary(ajust1); c(AIC(ajust1), AICc(ajust1), BIC(ajust1))
vif(ajust1); check_multicol_survreg(ajust1)

ajust2 = survreg(data = df_pse, Surv(time = t, event = ADOPT) ~
                   scale(log(PIB)) + scale(TX_VETO) + 
                   scale(log(DISCR_DOTACAO))
                 ,
                 dist = 'weibull')
summary(ajust2); c(AIC(ajust2), AICc(ajust2), BIC(ajust2))
vif(ajust2); check_multicol_survreg(ajust2)

(valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, .15))

vero_cop_debug = function(para) {
  cat("Parâmetros:", round(para, 6), "\n")
  
  qtd_par = length(para)
  if ( length(ajust1$coefficients) == length(ajust2$coefficients) ){
    if (qtd_par==7){ # 1 cov
      b01 = para[1]; b11 = para[2]; gama1 = para[3];
      b02 = para[4]; b12 = para[5]; gama2 = para[6];
      theta = para[7]
    } else if(qtd_par==9){ # 2 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; gama1 = para[4];
      b02 = para[5]; b12 = para[6]; b22 = para[7]; gama2 = para[8];
      theta = para[9]
    } else if(qtd_par==11){ # 3 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; gama1 = para[5];
      b02 = para[6]; b12 = para[7]; b22 = para[8]; b32 = para[9]; gama2 = para[10];
      theta = para[11]
    } else if(qtd_par==13){ # 4 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; b41 = para[5]; gama1 = para[6];
      b02 = para[7]; b12 = para[8]; b22 = para[9]; b32 = para[10]; b42 = para[11]; gama2 = para[12];
      theta = para[13]
    } else if(qtd_par==15){ # 5 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; b41 = para[5]; b51 = para[6]; gama1 = para[7];
      b02 = para[8]; b12 = para[9]; b22 = para[10]; b32 = para[11]; b42 = para[12]; b52 = para[13]; gama2 = para[14];
      theta = para[15]
    }
  } else { # nums de covs diferentes
    b01 = para[1]; b11 = para[2]; b21 = para[3]; gama1 = para[4];
    b02 = para[5]; b12 = para[6]; b22 = para[7]; b32 = para[8]; gama2 = para[9];
    theta = para[10]
  }
  
  
  if (gama1 <= 0 | gama2 <= 0) {
    cat("Gama não positivo:", gama1, gama2, "\n")
    return(1e10)
  }
  
  alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) +
                b31 * scale(df_pnab$DISCR_TX) ) #+ b41 * scale(log(df_pnab$SERVIDORES)) )
  # alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) +
  #               b31 * scale(df_pnab$DISCR_TX) + b41 * scale(log(df_pnab$SERVIDORES)) )
  #               #b51 * df_pnab$ANO_ELEITORAL_MUNIC)
  # alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) )
  # alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) )
  alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(df_pse$TX_VETO) +
                b32 * scale(log(df_pse$DISCR_DOTACAO)) ) #+ b42 * scale(df_pse$DISCR_TX) +
  #b52 * scale(log(df_pse$INDICE_MNST)) )
  # alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(log(df_pse$ORC_TOTAL_EMPENHADO)) +
  #               b32 * scale(log(df_pse$TX_SUCESSO)) + b42 * scale(log(df_pse$DISCR_TX)) +
  #               b52 * scale(log(df_pse$INDICE_MNST)) )
  # alfa = exp( b0 + b1*scale(df_pse$TX_VETO) + b2*df_pse$CARREIRAS +
  #               b3*scale(log(df_pse$DISCR_DOTACAO))+ b4*scale(df_pse$PIB))
  
  if ( any(!is.finite(alfa1)) || any(alfa1 <= 0) || any(!is.finite(alfa2)) || any(alfa2 <= 0) ) {
    cat("Alfa inválido\n")
    return(1e10)
  }
  
  log_vero_S = log_S_cop(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t1 = log_parc_t1(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t2 = log_parc_t2(alfa1, gama1, alfa2, gama2, theta)
  log_vero_f = log_f_cop(alfa1, gama1, alfa2, gama2, theta)
  
  result <- -sum(log_vero_S, log_vero_parc_t1, log_vero_parc_t2, log_vero_f)
  
  # cat("log_vero_S", sum(log_vero_S), '\n')
  # cat("log_vero_parc_t1", sum(log_vero_parc_t1), '\n')
  # cat("log_vero_parc_t1", sum(log_vero_parc_t2), '\n')
  # cat("log_vero_f", sum(log_vero_f), '\n')
  
  cat("Resultado:", ifelse(is.finite(result), round(result, 4), "Não finito"), "\n\n")
  
  if (!is.finite(result)) {
    return(1e10)
  }
  
  return(result)
}

valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, 0.1)
# nlminb é geralmente mais robusto que optim para problemas com bounds
MLE_nlminb <- nlminb(start = valorin1,
                     objective = vero_cop_debug,
                     # lower = c(rep(-Inf, 10), 0.00248),
                     # upper = c(100, 10, 100, 100, 10, 100, 4, 2),
                     control = list(iter.max = 2000, eval.max = 2000))
tau = MLE_nlminb[["par"]][[length(valorin1)]] / (2+MLE_nlminb[["par"]][[length(valorin1)]])
log_lik = -MLE_nlminb$objective; k = length(MLE_nlminb$par); n_observacoes = nrow(df_pnab)
aic = 2 * k - 2 * log_lik; aicc = aic + (2 * k * (k + 1)) / (n_observacoes - k - 1); bic = log(n_observacoes) * k - 2 * log_lik
print(MLE_nlminb); paste0('tau: ', round(tau, 4) )
paste0("AIC: ", round(aic, 3), "; AICc: ", round(aicc, 3), "; BIC: ", round(bic, 3), "; Vero: ", round(log_lik, 3) )

S_estim = function(t, cov1, cov2, cov3, b0, b1, b2, b3, gama){
  alfa = exp(b0 + b1*cov1 + b2*cov2 + b3*cov3)
  exp( -(t/alfa)^gama )
}
S_cop_estim = function(t1, t2, cov11, cov21, cov31, cov12, cov22, cov32,
                       b01, b11, b21, b31, gama1, b02, b12, b22, b32, gama2, theta){
  alfa1 = exp(b01 + b11*cov11 + b21*cov21 + b31*cov31); alfa2 = exp(b02 + b12*cov12 + b22*cov22 + b32*cov32)
  (S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1) ^(-1/theta)
}

# Definir a função f(x,y) com os parâmetros fixos
f <- function(x, y) {
  S_cop_estim(x, y, 0, 0, 0, 0, 0, 0, 
              MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
              MLE_nlminb$par[4], MLE_nlminb$par[5], MLE_nlminb$par[6],
              MLE_nlminb$par[7], MLE_nlminb$par[8], MLE_nlminb$par[9], 
              MLE_nlminb$par[10], MLE_nlminb$par[11])
}

f2 <- function(x, y) {
  S_estim(x, 0, 0, 0,
          MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
          MLE_nlminb$par[4], MLE_nlminb$par[5]) *
    S_estim(y, 0, 0, 0,
            MLE_nlminb$par[6], MLE_nlminb$par[7], MLE_nlminb$par[8], 
            MLE_nlminb$par[9], MLE_nlminb$par[10])
}

# Criar sequências de valores de 0 a 6000
x <- seq(0, 2800, length.out = 30)
y <- seq(0, 3500, length.out = 30)
# Calcular a matriz z
z <- outer(x, y, f)
z2 <- outer(x, y, f2)


### Gráfico para exportar ----

library(lattice)
# Preparar dados
df_surface <- expand.grid(x = x, y = y)
df_surface$z <- as.vector(z)


# Função para formatar números
formatar_numero <- function(x) {
  format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE, digits = 2)
}

# Gerar valores para eixos e legenda
x_vals <- pretty(df_surface$x)
y_vals <- pretty(df_surface$y)
z_vals <- pretty(df_surface$z)
colorkey_vals <- pretty(range(df_surface$z), n = 5)


# pdf("imagens/plano_mod1.pdf", 
#     width = 10, height = 8,  # Polegadas
#     pointsize = 12)

wireframe(z ~ x * y, data = df_surface,
          xlab = list("Tempos da PNAB", rot = -55, cex = 1), 
          ylab = list("Tempos do PSE", rot = 20, cex = 1),
          zlab = list(expression(paste("S(", t[1], ",", t[2], ")")), rot = 90, cex = 1),
          drape = TRUE,
          colorkey = list(
            title = expression(paste("S(", t[1], ",", t[2], ")")),
            width = 0.8,
            height = 0.6,
            at = colorkey_vals,
            labels = list(
              at = colorkey_vals,
              labels = formatar_numero(colorkey_vals)
            )
          ),
          scales = list(
            arrows = FALSE,
            x = list(at = x_vals, labels = formatar_numero(x_vals)),
            y = list(at = y_vals, labels = formatar_numero(y_vals)),
            z = list(at = z_vals, labels = formatar_numero(z_vals))
          ),
          screen = list(z = -60, x = -60),
          shade = FALSE,
          #col.regions = hcl.colors(100, "Viridis"),
          par.settings = list(axis.line = list(col = "transparent")))

# dev.off()


## Modelo 2 ----

ajust1 = survreg(data = df_pnab, Surv(time = t, event = ADOPT) ~
                   scale(log(PIB)) + scale(log(ORC_TOTAL_EMPENHADO)) 
                 + scale(DISCR_TX)
                 ,
                 dist = 'weibull')
summary(ajust1); c(AIC(ajust1), AICc(ajust1), BIC(ajust1))
vif(ajust1); check_multicol_survreg(ajust1)

ajust2 = survreg(data = df_pse, Surv(time = t, event = ADOPT) ~
                   scale(log(PIB)) + scale(TX_VETO) + 
                   scale(DISCR_TX)
                 ,
                 dist = 'weibull')
summary(ajust2); c(AIC(ajust2), AICc(ajust2), BIC(ajust2))
vif(ajust2); check_multicol_survreg(ajust2)

(valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, .15))

vero_cop_debug = function(para) {
  cat("Parâmetros:", round(para, 6), "\n")
  
  qtd_par = length(para)
  if ( length(ajust1$coefficients) == length(ajust2$coefficients) ){
    if (qtd_par==7){ # 1 cov
      b01 = para[1]; b11 = para[2]; gama1 = para[3];
      b02 = para[4]; b12 = para[5]; gama2 = para[6];
      theta = para[7]
    } else if(qtd_par==9){ # 2 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; gama1 = para[4];
      b02 = para[5]; b12 = para[6]; b22 = para[7]; gama2 = para[8];
      theta = para[9]
    } else if(qtd_par==11){ # 3 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; gama1 = para[5];
      b02 = para[6]; b12 = para[7]; b22 = para[8]; b32 = para[9]; gama2 = para[10];
      theta = para[11]
    } else if(qtd_par==13){ # 4 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; b41 = para[5]; gama1 = para[6];
      b02 = para[7]; b12 = para[8]; b22 = para[9]; b32 = para[10]; b42 = para[11]; gama2 = para[12];
      theta = para[13]
    } else if(qtd_par==15){ # 5 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; b41 = para[5]; b51 = para[6]; gama1 = para[7];
      b02 = para[8]; b12 = para[9]; b22 = para[10]; b32 = para[11]; b42 = para[12]; b52 = para[13]; gama2 = para[14];
      theta = para[15]
    }
  } else { # nums de covs diferentes
    b01 = para[1]; b11 = para[2]; b21 = para[3]; gama1 = para[4];
    b02 = para[5]; b12 = para[6]; b22 = para[7]; b32 = para[8]; gama2 = para[9];
    theta = para[10]
  }
  
  
  if (gama1 <= 0 | gama2 <= 0) {
    cat("Gama não positivo:", gama1, gama2, "\n")
    return(1e10)
  }
  
  alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) +
                b31 * scale(df_pnab$DISCR_TX) ) #+ b41 * scale(log(df_pnab$SERVIDORES)) )
  # alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) +
  #               b31 * scale(df_pnab$DISCR_TX) + b41 * scale(log(df_pnab$SERVIDORES)) )
  #               #b51 * df_pnab$ANO_ELEITORAL_MUNIC)
  # alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) )
  # alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) )
  alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(df_pse$TX_VETO) +
                b32 * scale(df_pse$DISCR_TX) ) #+ b42 * scale(df_pse$DISCR_TX) +
  #b52 * scale(log(df_pse$INDICE_MNST)) )
  # alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(log(df_pse$ORC_TOTAL_EMPENHADO)) +
  #               b32 * scale(log(df_pse$TX_SUCESSO)) + b42 * scale(log(df_pse$DISCR_TX)) +
  #               b52 * scale(log(df_pse$INDICE_MNST)) )
  # alfa = exp( b0 + b1*scale(df_pse$TX_VETO) + b2*df_pse$CARREIRAS +
  #               b3*scale(log(df_pse$DISCR_DOTACAO))+ b4*scale(df_pse$PIB))
  
  if ( any(!is.finite(alfa1)) || any(alfa1 <= 0) || any(!is.finite(alfa2)) || any(alfa2 <= 0) ) {
    cat("Alfa inválido\n")
    return(1e10)
  }
  
  log_vero_S = log_S_cop(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t1 = log_parc_t1(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t2 = log_parc_t2(alfa1, gama1, alfa2, gama2, theta)
  log_vero_f = log_f_cop(alfa1, gama1, alfa2, gama2, theta)
  
  result <- -sum(log_vero_S, log_vero_parc_t1, log_vero_parc_t2, log_vero_f)
  
  # cat("log_vero_S", sum(log_vero_S), '\n')
  # cat("log_vero_parc_t1", sum(log_vero_parc_t1), '\n')
  # cat("log_vero_parc_t1", sum(log_vero_parc_t2), '\n')
  # cat("log_vero_f", sum(log_vero_f), '\n')
  
  cat("Resultado:", ifelse(is.finite(result), round(result, 4), "Não finito"), "\n\n")
  
  if (!is.finite(result)) {
    return(1e10)
  }
  
  return(result)
}

valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, 0.1)
# nlminb é geralmente mais robusto que optim para problemas com bounds
MLE_nlminb <- nlminb(start = valorin1,
                     objective = vero_cop_debug,
                     # lower = c(rep(-Inf, 10), 0.00248),
                     # upper = c(100, 10, 100, 100, 10, 100, 4, 2),
                     control = list(iter.max = 2000, eval.max = 2000))
tau = MLE_nlminb[["par"]][[length(valorin1)]] / (2+MLE_nlminb[["par"]][[length(valorin1)]])
log_lik = -MLE_nlminb$objective; k = length(MLE_nlminb$par); n_observacoes = nrow(df_pnab)
aic = 2 * k - 2 * log_lik; aicc = aic + (2 * k * (k + 1)) / (n_observacoes - k - 1); bic = log(n_observacoes) * k - 2 * log_lik
print(MLE_nlminb); paste0('tau: ', round(tau, 4) )
paste0("AIC: ", round(aic, 3), "; AICc: ", round(aicc, 3), "; BIC: ", round(bic, 3), "; Vero: ", round(log_lik, 3) )


S_estim = function(t, cov1, cov2, cov3, b0, b1, b2, b3, gama){
  alfa = exp(b0 + b1*cov1 + b2*cov2 + b3*cov3)
  exp( -(t/alfa)^gama )
}
S_cop_estim = function(t1, t2, cov11, cov21, cov31, cov12, cov22, cov32,
                       b01, b11, b21, b31, gama1, b02, b12, b22, b32, gama2, theta){
  alfa1 = exp(b01 + b11*cov11 + b21*cov21 + b31*cov31); alfa2 = exp(b02 + b12*cov12 + b22*cov22 + b32*cov32)
  (S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1) ^(-1/theta)
}
# Definir a função f(x,y) com os parâmetros fixos
f <- function(x, y) {
  S_cop_estim(x, y, 0, 0, 0, 0, 0, 0, 
              MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
              MLE_nlminb$par[4], MLE_nlminb$par[5], MLE_nlminb$par[6],
              MLE_nlminb$par[7], MLE_nlminb$par[8], MLE_nlminb$par[9], 
              MLE_nlminb$par[10], MLE_nlminb$par[11])
}

f2 <- function(x, y) {
  S_estim(x, 0, 0, 0,
          MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
          MLE_nlminb$par[4], MLE_nlminb$par[5]) *
    S_estim(y, 0, 0, 0,
            MLE_nlminb$par[6], MLE_nlminb$par[7], MLE_nlminb$par[8], 
            MLE_nlminb$par[9], MLE_nlminb$par[10])
}

# Criar sequências de valores de 0 a 6000
x <- seq(0, 2800, length.out = 30)
y <- seq(0, 3500, length.out = 30)
# Calcular a matriz z
z <- outer(x, y, f)
z2 <- outer(x, y, f2)

### Gráfico para exportar ----

library(lattice)
# Preparar dados
df_surface <- expand.grid(x = x, y = y)
df_surface$z <- as.vector(z)


# Função para formatar números
formatar_numero <- function(x) {
  format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE, digits = 2)
}

# Gerar valores para eixos e legenda
x_vals <- pretty(df_surface$x)
y_vals <- pretty(df_surface$y)
z_vals <- pretty(df_surface$z)
colorkey_vals <- pretty(range(df_surface$z), n = 5)


# pdf("imagens/plano_mod2.pdf", 
#     width = 10, height = 8,  # Polegadas
#     pointsize = 12)

wireframe(z ~ x * y, data = df_surface,
          xlab = list("Tempos da PNAB", rot = -55, cex = 1), 
          ylab = list("Tempos do PSE", rot = 20, cex = 1),
          zlab = list(expression(paste("S(", t[1], ",", t[2], ")")), rot = 90, cex = 1),
          drape = TRUE,
          colorkey = list(
            title = expression(paste("S(", t[1], ",", t[2], ")")),
            width = 0.8,
            height = 0.6,
            at = colorkey_vals,
            labels = list(
              at = colorkey_vals,
              labels = formatar_numero(colorkey_vals)
            )
          ),
          scales = list(
            arrows = FALSE,
            x = list(at = x_vals, labels = formatar_numero(x_vals)),
            y = list(at = y_vals, labels = formatar_numero(y_vals)),
            z = list(at = z_vals, labels = formatar_numero(z_vals))
          ),
          screen = list(z = -60, x = -60),
          shade = FALSE,
          #col.regions = hcl.colors(100, "Viridis"),
          par.settings = list(axis.line = list(col = "transparent")))

# dev.off()



## Modelo 3 ----

ajust1 = survreg(data = df_pnab, Surv(time = t, event = ADOPT) ~
                   scale(log(PIB)) + scale(log(ORC_TOTAL_EMPENHADO))
                 ,
                 dist = 'weibull')
summary(ajust1); c(AIC(ajust1), AICc(ajust1), BIC(ajust1))
vif(ajust1); check_multicol_survreg(ajust1)

ajust2 = survreg(data = df_pse, Surv(time = t, event = ADOPT) ~
                   scale(log(PIB)) + scale(TX_VETO) +
                   scale(log(DISCR_DOTACAO))
                 
                 ,
                 dist = 'weibull')
summary(ajust2); c(AIC(ajust2), AICc(ajust2), BIC(ajust2))
vif(ajust2); check_multicol_survreg(ajust2)

(valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, .15))

vero_cop_debug = function(para) {
  cat("Parâmetros:", round(para, 6), "\n")
  
  qtd_par = length(para)
  if ( length(ajust1$coefficients) == length(ajust2$coefficients) ){
    if (qtd_par==7){ # 1 cov
      b01 = para[1]; b11 = para[2]; gama1 = para[3];
      b02 = para[4]; b12 = para[5]; gama2 = para[6];
      theta = para[7]
    } else if(qtd_par==9){ # 2 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; gama1 = para[4];
      b02 = para[5]; b12 = para[6]; b22 = para[7]; gama2 = para[8];
      theta = para[9]
    } else if(qtd_par==11){ # 3 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; gama1 = para[5];
      b02 = para[6]; b12 = para[7]; b22 = para[8]; b32 = para[9]; gama2 = para[10];
      theta = para[11]
    } else if(qtd_par==13){ # 4 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; b41 = para[5]; gama1 = para[6];
      b02 = para[7]; b12 = para[8]; b22 = para[9]; b32 = para[10]; b42 = para[11]; gama2 = para[12];
      theta = para[13]
    } else if(qtd_par==15){ # 5 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; b41 = para[5]; b51 = para[6]; gama1 = para[7];
      b02 = para[8]; b12 = para[9]; b22 = para[10]; b32 = para[11]; b42 = para[12]; b52 = para[13]; gama2 = para[14];
      theta = para[15]
    }
  } else { # nums de covs diferentes
    b01 = para[1]; b11 = para[2]; b21 = para[3]; gama1 = para[4];
    b02 = para[5]; b12 = para[6]; b22 = para[7]; b32 = para[8]; gama2 = para[9];
    theta = para[10]
  }
  
  
  if (gama1 <= 0 | gama2 <= 0) {
    cat("Gama não positivo:", gama1, gama2, "\n")
    return(1e10)
  }
  
  alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) ) #+
  #b31 * scale(log(df_pnab$DISCR_DOTACAO)) ) #+ b41 * scale(log(df_pnab$SERVIDORES)) )
  # alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) +
  #               b31 * scale(df_pnab$DISCR_TX) + b41 * scale(log(df_pnab$SERVIDORES)) )
  #               #b51 * df_pnab$ANO_ELEITORAL_MUNIC)
  # alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) )
  # alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) )
  alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(df_pse$TX_VETO) +
                b32 * scale(log(df_pse$DISCR_DOTACAO)) ) #+ b42 * scale(df_pse$DISCR_TX) +
  #b52 * scale(log(df_pse$INDICE_MNST)) )
  # alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(log(df_pse$ORC_TOTAL_EMPENHADO)) +
  #               b32 * scale(log(df_pse$TX_SUCESSO)) + b42 * scale(log(df_pse$DISCR_TX)) +
  #               b52 * scale(log(df_pse$INDICE_MNST)) )
  # alfa = exp( b0 + b1*scale(df_pse$TX_VETO) + b2*df_pse$CARREIRAS +
  #               b3*scale(log(df_pse$DISCR_DOTACAO))+ b4*scale(df_pse$PIB))
  
  if ( any(!is.finite(alfa1)) || any(alfa1 <= 0) || any(!is.finite(alfa2)) || any(alfa2 <= 0) ) {
    cat("Alfa inválido\n")
    return(1e10)
  }
  
  log_vero_S = log_S_cop(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t1 = log_parc_t1(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t2 = log_parc_t2(alfa1, gama1, alfa2, gama2, theta)
  log_vero_f = log_f_cop(alfa1, gama1, alfa2, gama2, theta)
  
  result <- -sum(log_vero_S, log_vero_parc_t1, log_vero_parc_t2, log_vero_f)
  
  # cat("log_vero_S", sum(log_vero_S), '\n')
  # cat("log_vero_parc_t1", sum(log_vero_parc_t1), '\n')
  # cat("log_vero_parc_t1", sum(log_vero_parc_t2), '\n')
  # cat("log_vero_f", sum(log_vero_f), '\n')
  
  cat("Resultado:", ifelse(is.finite(result), round(result, 4), "Não finito"), "\n\n")
  
  if (!is.finite(result)) {
    return(1e10)
  }
  
  return(result)
}

valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, 0.1)
# nlminb é geralmente mais robusto que optim para problemas com bounds
MLE_nlminb <- nlminb(start = valorin1,
                     objective = vero_cop_debug,
                     # lower = c(rep(-Inf, 10), 0.00248),
                     # upper = c(100, 10, 100, 100, 10, 100, 4, 2),
                     control = list(iter.max = 2000, eval.max = 2000))
tau = MLE_nlminb[["par"]][[length(valorin1)]] / (2+MLE_nlminb[["par"]][[length(valorin1)]])
log_lik = -MLE_nlminb$objective; k = length(MLE_nlminb$par); n_observacoes = nrow(df_pnab)
aic = 2 * k - 2 * log_lik; aicc = aic + (2 * k * (k + 1)) / (n_observacoes - k - 1); bic = log(n_observacoes) * k - 2 * log_lik
print(MLE_nlminb); paste0('tau: ', round(tau, 4) )
paste0("AIC: ", round(aic, 3), "; AICc: ", round(aicc, 3), "; BIC: ", round(bic, 3), "; Vero: ", round(log_lik, 3) )


S_estim1 = function(t, cov1, cov2, b0, b1, b2, gama){
  alfa = exp(b0 + b1*cov1 + b2*cov2)
  exp( -(t/alfa)^gama )
}
S_estim2 = function(t, cov1, cov2, cov3, b0, b1, b2, b3, gama){
  alfa = exp(b0 + b1*cov1 + b2*cov2 + b3*cov3)
  exp( -(t/alfa)^gama )
}
S_cop_estim = function(t1, t2, cov11, cov21, cov12, cov22, cov32,
                       b01, b11, b21, gama1, b02, b12, b22, b32, gama2, theta){
  alfa1 = exp(b01 + b11*cov11 + b21*cov21); alfa2 = exp(b02 + b12*cov12 + b22*cov22 + b32*cov32)
  (S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1) ^(-1/theta)
}
# Definir a função f(x,y) com os parâmetros fixos
f <- function(x, y) {
  S_cop_estim(x, y, 0, 0, 0, 0, 0, 
              MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
              MLE_nlminb$par[4], MLE_nlminb$par[5], MLE_nlminb$par[6],
              MLE_nlminb$par[7], MLE_nlminb$par[8], MLE_nlminb$par[9], 
              MLE_nlminb$par[10])
}

f2 <- function(x, y) {
  S_estim1(x, 0, 0,
           MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
           MLE_nlminb$par[4]) *
    S_estim2(y, 0, 0, 0,
             MLE_nlminb$par[5], MLE_nlminb$par[6], MLE_nlminb$par[7], 
             MLE_nlminb$par[8], MLE_nlminb$par[9])
}

# Criar sequências de valores de 0 a 6000
x <- seq(0, 2800, length.out = 30)
y <- seq(0, 3500, length.out = 30)
# Calcular a matriz z
z <- outer(x, y, f)
z2 <- outer(x, y, f2)

### Gráfico para exportar ----

library(lattice)
# Preparar dados
df_surface <- expand.grid(x = x, y = y)
df_surface$z <- as.vector(z)


# Função para formatar números
formatar_numero <- function(x) {
  format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE, digits = 2)
}

# Gerar valores para eixos e legenda
x_vals <- pretty(df_surface$x)
y_vals <- pretty(df_surface$y)
z_vals <- pretty(df_surface$z)
colorkey_vals <- pretty(range(df_surface$z), n = 5)


# pdf("imagens/plano_mod3.pdf", 
#     width = 10, height = 8,  # Polegadas
#     pointsize = 12)

wireframe(z ~ x * y, data = df_surface,
          xlab = list("Tempos da PNAB", rot = -55, cex = 1), 
          ylab = list("Tempos do PSE", rot = 20, cex = 1),
          zlab = list(expression(paste("S(", t[1], ",", t[2], ")")), rot = 90, cex = 1),
          drape = TRUE,
          colorkey = list(
            title = expression(paste("S(", t[1], ",", t[2], ")")),
            width = 0.8,
            height = 0.6,
            at = colorkey_vals,
            labels = list(
              at = colorkey_vals,
              labels = formatar_numero(colorkey_vals)
            )
          ),
          scales = list(
            arrows = FALSE,
            x = list(at = x_vals, labels = formatar_numero(x_vals)),
            y = list(at = y_vals, labels = formatar_numero(y_vals)),
            z = list(at = z_vals, labels = formatar_numero(z_vals))
          ),
          screen = list(z = -60, x = -60),
          shade = FALSE,
          #col.regions = hcl.colors(100, "Viridis"),
          par.settings = list(axis.line = list(col = "transparent")))

# dev.off()



## Modelo 4 ----

ajust1 = survreg(data = df_pnab, Surv(time = t, event = ADOPT) ~
                   scale(log(PIB)) + scale(log(ORC_TOTAL_EMPENHADO)) +
                   scale(DISCR_TX) + scale(log(SERVIDORES))
                 
                 
                 ,
                 dist = 'weibull')
summary(ajust1); c(AIC(ajust1), AICc(ajust1), BIC(ajust1))
vif(ajust1); check_multicol_survreg(ajust1)

ajust2 = survreg(data = df_pse, Surv(time = t, event = ADOPT) ~
                   scale(log(PIB)) + scale(ORC_TOTAL_EMPENHADO) +
                   scale(log(TX_SUCESSO)) + scale(log(DISCR_TX)) +
                   scale(log(INDICE_MNST))
                 
                 
                 ,
                 dist = 'weibull')
summary(ajust2); c(AIC(ajust2), AICc(ajust2), BIC(ajust2))
vif(ajust2); check_multicol_survreg(ajust2)

(valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, .15))

vero_cop_debug = function(para) {
  cat("Parâmetros:", round(para, 6), "\n")
  
  qtd_par = length(para)
  if ( length(ajust1$coefficients) == length(ajust2$coefficients) ){
    if (qtd_par==7){ # 1 cov
      b01 = para[1]; b11 = para[2]; gama1 = para[3];
      b02 = para[4]; b12 = para[5]; gama2 = para[6];
      theta = para[7]
    } else if(qtd_par==9){ # 2 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; gama1 = para[4];
      b02 = para[5]; b12 = para[6]; b22 = para[7]; gama2 = para[8];
      theta = para[9]
    } else if(qtd_par==11){ # 3 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; gama1 = para[5];
      b02 = para[6]; b12 = para[7]; b22 = para[8]; b32 = para[9]; gama2 = para[10];
      theta = para[11]
    } else if(qtd_par==13){ # 4 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; b41 = para[5]; gama1 = para[6];
      b02 = para[7]; b12 = para[8]; b22 = para[9]; b32 = para[10]; b42 = para[11]; gama2 = para[12];
      theta = para[13]
    } else if(qtd_par==15){ # 5 covs
      b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; b41 = para[5]; b51 = para[6]; gama1 = para[7];
      b02 = para[8]; b12 = para[9]; b22 = para[10]; b32 = para[11]; b42 = para[12]; b52 = para[13]; gama2 = para[14];
      theta = para[15]
    }
  } else { # nums de covs diferentes
    b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; b41 = para[5]; gama1 = para[6];
    b02 = para[7]; b12 = para[8]; b22 = para[9]; b32 = para[10]; b42 = para[11]; b52 = para[12]; gama2 = para[13];
    theta = para[14]
  }
  
  
  if (gama1 <= 0 | gama2 <= 0) {
    cat("Gama não positivo:", gama1, gama2, "\n")
    return(1e10)
  }
  
  alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) +
                b31 * scale(df_pnab$DISCR_TX) + b41 * scale(log(df_pnab$SERVIDORES)) )
  # alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) +
  #               b31 * scale(df_pnab$DISCR_TX) + b41 * scale(log(df_pnab$SERVIDORES)) )
  #               #b51 * df_pnab$ANO_ELEITORAL_MUNIC)
  # alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) )
  # alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) )
  alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(df_pse$ORC_TOTAL_EMPENHADO) +
                b32 * scale(log(df_pse$TX_SUCESSO)) + b42 * scale(log(df_pse$DISCR_TX)) +
                b52 * scale(log(df_pse$INDICE_MNST)) )
  # alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(log(df_pse$ORC_TOTAL_EMPENHADO)) +
  #               b32 * scale(log(df_pse$TX_SUCESSO)) + b42 * scale(log(df_pse$DISCR_TX)) +
  #               b52 * scale(log(df_pse$INDICE_MNST)) )
  # alfa = exp( b0 + b1*scale(df_pse$TX_VETO) + b2*df_pse$CARREIRAS +
  #               b3*scale(log(df_pse$DISCR_DOTACAO))+ b4*scale(df_pse$PIB))
  
  if ( any(!is.finite(alfa1)) || any(alfa1 <= 0) || any(!is.finite(alfa2)) || any(alfa2 <= 0) ) {
    cat("Alfa inválido\n")
    return(1e10)
  }
  
  log_vero_S = log_S_cop(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t1 = log_parc_t1(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t2 = log_parc_t2(alfa1, gama1, alfa2, gama2, theta)
  log_vero_f = log_f_cop(alfa1, gama1, alfa2, gama2, theta)
  
  result <- -sum(log_vero_S, log_vero_parc_t1, log_vero_parc_t2, log_vero_f)
  
  # cat("log_vero_S", sum(log_vero_S), '\n')
  # cat("log_vero_parc_t1", sum(log_vero_parc_t1), '\n')
  # cat("log_vero_parc_t1", sum(log_vero_parc_t2), '\n')
  # cat("log_vero_f", sum(log_vero_f), '\n')
  
  cat("Resultado:", ifelse(is.finite(result), round(result, 4), "Não finito"), "\n\n")
  
  if (!is.finite(result)) {
    return(1e10)
  }
  
  return(result)
}

valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, 0.1)
# nlminb é geralmente mais robusto que optim para problemas com bounds
MLE_nlminb <- nlminb(start = valorin1,
                     objective = vero_cop_debug,
                     # lower = c(rep(-Inf, 10), 0.00248),
                     # upper = c(100, 10, 100, 100, 10, 100, 4, 2),
                     control = list(iter.max = 2000, eval.max = 2000))
tau = MLE_nlminb[["par"]][[length(valorin1)]] / (2+MLE_nlminb[["par"]][[length(valorin1)]])
log_lik = -MLE_nlminb$objective; k = length(MLE_nlminb$par); n_observacoes = nrow(df_pnab)
aic = 2 * k - 2 * log_lik; aicc = aic + (2 * k * (k + 1)) / (n_observacoes - k - 1); bic = log(n_observacoes) * k - 2 * log_lik
print(MLE_nlminb); paste0('tau: ', round(tau, 4) )
paste0("AIC: ", round(aic, 3), "; AICc: ", round(aicc, 3), "; BIC: ", round(bic, 3), "; Vero: ", round(log_lik, 3) )


S_estim1 = function(t, cov1, cov2, cov3, cov4, b0, b1, b2, b3, b4, gama){
  alfa = exp(b0 + b1*cov1 + b2*cov2 + b3*cov3 + b4*cov4)
  exp( -(t/alfa)^gama )
}
S_estim2 = function(t, cov1, cov2, cov3, cov4, cov5, b0, b1, b2, b3, b4, b5, gama){
  alfa = exp(b0 + b1*cov1 + b2*cov2 + b3*cov3 + b4*cov4 + b5*cov5)
  exp( -(t/alfa)^gama )
}
S_cop_estim = function(t1, t2, cov11, cov21, cov31, cov41, cov12, cov22, cov32, cov42, cov52,
                       b01, b11, b21, b31, b41, gama1, b02, b12, b22, b32, b42, b52, gama2, theta){
  alfa1 = exp(b01 + b11*cov11 + b21*cov21 + b31*cov31 + b41*cov41)
  alfa2 = exp(b02 + b12*cov12 + b22*cov22 + b32*cov32 + b42*cov42 + b52*cov52)
  (S_theta(t1, alfa1, gama1, theta) + S_theta(t2, alfa2, gama2, theta) - 1) ^(-1/theta)
}
# Definir a função f(x,y) com os parâmetros fixos
f <- function(x, y) {
  S_cop_estim(x, y, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
              MLE_nlminb$par[4], MLE_nlminb$par[5], MLE_nlminb$par[6],
              MLE_nlminb$par[7], MLE_nlminb$par[8], MLE_nlminb$par[9], 
              MLE_nlminb$par[10], MLE_nlminb$par[11], MLE_nlminb$par[12],
              MLE_nlminb$par[13], MLE_nlminb$par[14])
}

f2 <- function(x, y) {
  S_estim1(x, 0, 0, 0, 0,
           MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
           MLE_nlminb$par[4], MLE_nlminb$par[5], MLE_nlminb$par[6]) *
    S_estim2(y, 0, 0, 0, 0, 0,
             MLE_nlminb$par[7], MLE_nlminb$par[8], MLE_nlminb$par[9], 
             MLE_nlminb$par[10], MLE_nlminb$par[11], MLE_nlminb$par[12], 
             MLE_nlminb$par[13])
}

# Criar sequências de valores de 0 a 6000
x <- seq(0, 2800, length.out = 30)
y <- seq(0, 2500, length.out = 30)
# Calcular a matriz z
z <- outer(x, y, f); z[is.na(z)] <- 0
z2 <- outer(x, y, f2)


### Gráfico para exportar ----

library(lattice)
# Preparar dados
df_surface <- expand.grid(x = x, y = y)
df_surface$z <- as.vector(z)


# Função para formatar números
formatar_numero <- function(x) {
  format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE, digits = 2)
}

# Gerar valores para eixos e legenda
x_vals <- pretty(df_surface$x)
y_vals <- pretty(df_surface$y)
z_vals <- pretty(df_surface$z)
colorkey_vals <- pretty(range(df_surface$z), n = 5)


# pdf("imagens/plano_mod4.pdf", 
#     width = 10, height = 8,  # Polegadas
#     pointsize = 12)

wireframe(z ~ x * y, data = df_surface,
          xlab = list("Tempos da PNAB", rot = 30, cex = 1), 
          ylab = list("Tempos do PSE", rot = -12, cex = 1),
          zlab = list(expression(paste("S(", t[1], ",", t[2], ")")), rot = 90, cex = 1),
          drape = TRUE,
          colorkey = list(
            title = expression(paste("S(", t[1], ",", t[2], ")")),
            width = 0.8,
            height = 0.6,
            at = colorkey_vals,
            labels = list(
              at = colorkey_vals,
              labels = formatar_numero(colorkey_vals)
            )
          ),
          scales = list(
            arrows = FALSE,
            x = list(at = x_vals, labels = formatar_numero(x_vals)),
            y = list(at = y_vals, labels = formatar_numero(y_vals)),
            z = list(at = z_vals, labels = formatar_numero(z_vals))
          ),
          screen = list(z = -120, x = -60),
          #screen = list(z = 120, x = -60),  # z mudou para 120
          shade = FALSE,
          #col.regions = hcl.colors(100, "Viridis"),
          par.settings = list(axis.line = list(col = "transparent")))

# dev.off()



# Residuos ----

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
ggplot() +
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
  labs(
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
  )
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
ggplot() +
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
  labs(
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
  )
#ggsave("imagens/residuos_pse_mod3.png", width = 158, height = 93, units = "mm", bg = "white")



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
ggplot() +
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
  labs(
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
  )
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
ggplot() +
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
  labs(
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
  )
#ggsave("imagens/residuos_pse_mod4.png", width = 158, height = 93, units = "mm", bg = "white")

