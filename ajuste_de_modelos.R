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
library(msm)

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


# PNAB ----

variaveis = variaveis[variaveis != "NUM_ENTIDADES"]

df_pnab = base_mod |> filter(POL_ABV == 'pnab')
df_pnab$NUM_ENTIDADES = NULL

mod_pnab_teste = survreg(data = df_pnab, Surv(time = t, event = ADOPT) ~
                            scale(log(PIB))
                         + scale(log(ORC_TOTAL_EMPENHADO))
                          # scale(log(PIB)) + scale(log(ORC_TOTAL_EMPENHADO)) +
                          #  scale(DISCR_TX) + scale(log(SERVIDORES))
                         # + ANO_ELEITORAL_MUNIC
                         # + IDEOL_PRF
                           #scale(PIB) + scale(log(ORC_TOTAL_EMPENHADO)) #+
                           #scale(log(DISCR_TX)) + scale(log(SERVIDORES))
                         
                         
                         
                         
                         
                         ,
                         dist = 'weibull')
summary(mod_pnab_teste); c(AIC(mod_pnab_teste), AICc(mod_pnab_teste), BIC(mod_pnab_teste))
vif(mod_pnab_teste)
check_multicol_survreg(mod_pnab_teste)

sumario <- summary(mod_pnab_teste)
# Extrair a tabela de coeficientes
sumario$table

# Aplicar o método delta
# 1/scale = exp(-log(scale))
se_inv_scale <- deltamethod(
  ~ exp(-x1),
  mean = log(mod_pnab_teste$scale),
  cov = vcov(mod_pnab_teste)["Log(scale)", "Log(scale)"]
)
# Resultado
cat("Forma da Weibull (1/scale):", 1/mod_pnab_teste$scale, "\n")
cat("Erro padrão:", se_inv_scale, "\n")

cov1 = scale(log(df_pnab$DISCR_TX)); cov2 = scale(df_pnab$ANO_ELEITORAL_NAC)
cor(cov1, cov2); cor(cov1, cov2, method ="kendall")

# Análise de resíduos
# 1. Resíduos de Cox-Snell
xitbeta <- mod_pnab_teste$linear.predictors; scale_weibull <- mod_pnab_teste$scale
# Calcular resíduos de Cox-Snell para Weibull
residuos_cs <- (df_pnab$t / exp(xitbeta))^(1/scale_weibull)
# 2. Curva de sobrevivência dos resíduos
surv_residuos <- Surv(time = residuos_cs, event = df_pnab$ADOPT)
res_fit <- survfit(surv_residuos ~ 1)
# 3. Plot principal
plot(res_fit, conf.int = FALSE, 
     ylab = "S(ei) - Função de Sobrevivência", 
     xlab = "ei - Resíduos de Cox-Snell",
     main = "Ajuste do Modelo Weibull - Resíduos Cox-Snell")
# Linha de referência Exponencial(1)
curve(exp(-x), from = 0, to = max(residuos_cs, na.rm = TRUE), 
      col = "red", lwd = 2, add = TRUE)
legend("topright", 
       legend = c("Resíduos Observados (KM)", "Exp(1) - Referência"),
       col = c("black", "red"), lwd = 2, bty = "n")
# Curva de sobrevivência Kaplan-Meier dos resíduos
surv_residuos <- Surv(time = residuos_cs, event = df_pnab$ADOPT)
km_residuos <- survfit(surv_residuos ~ 1)
# Gerar o gráfico comparativo
plot(exp(-km_residuos$time), km_residuos$surv,
     type = "l", lwd = 2,
     xlab = expression(exp(-hat(e[i]))),
     ylab = expression(widehat(S[KM])(hat(e[i]))),
     main = "Verificação do Ajuste Weibull\nResíduos de Cox-Snell",
     xlim = c(0, 1), ylim = c(0, 1), col = "blue")
# Adicionar linha de referência (y = x)
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
# Legenda
legend("bottomright", 
       legend = c("Curva Observada (KM)", "Linha de Referência (y = x)"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n"); grid()



# PSE ----

df_pse = base_mod |> filter(POL_ABV == 'pse')
df_pse$NUM_ENTIDADES = NULL

mod_pse_teste = survreg(data = df_pse, Surv(time = t, event = ADOPT) ~
                          # scale(log(PIB)) + scale(TX_VETO) +
                          # scale(log(DISCR_DOTACAO))
                          scale(log(PIB)) + scale(ORC_TOTAL_EMPENHADO) +
                          scale(log(TX_SUCESSO)) + scale(log(DISCR_TX)) +
                          scale(log(INDICE_MNST))
                        
                        ,
                         dist = 'weibull')
summary(mod_pse_teste); c(AIC(mod_pse_teste), AICc(mod_pse_teste), BIC(mod_pse_teste))
vif(mod_pse_teste)
check_multicol_survreg(mod_pse_teste)

sumario <- summary(mod_pse_teste)
# Extrair a tabela de coeficientes
sumario$table

# Aplicar o método delta
# 1/scale = exp(-log(scale))
se_inv_scale <- deltamethod(
  ~ exp(-x1),
  mean = log(mod_pse_teste$scale),
  cov = vcov(mod_pse_teste)["Log(scale)", "Log(scale)"]
)
# Resultado
cat("Forma da Weibull (1/scale):", 1/mod_pse_teste$scale, "\n")
cat("Erro padrão:", se_inv_scale, "\n")


# Análise de resíduos
# 1. Resíduos de Cox-Snell
xitbeta <- mod_pse_teste$linear.predictors; scale_weibull <- mod_pse_teste$scale
# Calcular resíduos de Cox-Snell para Weibull
residuos_cs <- (df_pse$t / exp(xitbeta))^(1/scale_weibull)
# 2. Curva de sobrevivência dos resíduos
surv_residuos <- Surv(time = residuos_cs, event = df_pse$ADOPT)
res_fit <- survfit(surv_residuos ~ 1)
# 3. Plot principal
plot(res_fit, conf.int = FALSE, 
     ylab = "S(ei) - Função de Sobrevivência", 
     xlab = "ei - Resíduos de Cox-Snell",
     main = "Ajuste do Modelo Weibull - Resíduos Cox-Snell")
# Linha de referência Exponencial(1)
curve(exp(-x), from = 0, to = max(residuos_cs, na.rm = TRUE), 
      col = "red", lwd = 2, add = TRUE)
 legend("topright", 
       legend = c("Resíduos Observados (KM)", "Exp(1) - Referência"),
       col = c("black", "red"), lwd = 2, bty = "n")
# Curva de sobrevivência Kaplan-Meier dos resíduos
surv_residuos <- Surv(time = residuos_cs, event = df_pse$ADOPT)
km_residuos <- survfit(surv_residuos ~ 1)
# Gerar o gráfico comparativo
plot(exp(-km_residuos$time), km_residuos$surv,
     type = "l", lwd = 2,
     xlab = expression(exp(-hat(e[i]))),
     ylab = expression(widehat(S[KM])(hat(e[i]))),
     main = "Verificação do Ajuste Weibull\nResíduos de Cox-Snell",
     xlim = c(0, 1), ylim = c(0, 1), col = "blue")
# Adicionar linha de referência (y = x)
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
# Legenda
legend("bottomright", 
       legend = c("Curva Observada (KM)", "Linha de Referência (y = x)"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n"); grid()



# cópulas - ajuste de Breno ----

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



## cópulas sem covariáveis ----

vero_cop_debug = function(para) {
  cat("Parâmetros:", round(para, 6), "\n")
  
  # b01 = para[1]; b11 = para[2]; gama1 = para[3]; b02 = para[4]; b12 = para[5]; gama2 = para[6]; theta = para[7]
  alfa1 = para[1]; gama1 = para[2]; alfa2 = para[3]; gama2 = para[4]; theta = para[5]
  
  if (gama1 <= 0 | gama2 <= 0) {
    cat("Gama não positivo:", gama1, gama2, "\n")
    return(1e10)
  }
  
  #alfa1 = exp(b01 + b11 * log(df_pnab$ORC_TOTAL_EMPENHADO) )
  #alfa1 = exp(b01 + b11 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO) ) )
  #alfa2 = exp(b02 + b12 * scale(df_pse$TX_VETO) )
  #alfa2 = exp(b02 + b12 * scale(log(df_pse$ORC_TOTAL_EMPENHADO) ) )
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


# nlminb é geralmente mais robusto que optim para problemas com bounds
MLE_nlminb <- nlminb(start = 
                       c(2100, 1, 2100, 1, .8), # sem covs 
                     objective = vero_cop_debug,
                     control = list(iter.max = 200, eval.max = 2000))
print(MLE_nlminb)


# cópulas com covariáveis ----

## função para teste de modelos ----

# Análise de correlação entre covariáveis
cov1pib = scale(df_pnab$PIB); cov2pib = scale(df_pse$PIB)
cov1 = scale(df_pnab$TX_VETO); cov2 = scale(df_pse$TX_VETO)
cov1 = scale(log(df_pnab$ORC_TOTAL_EMPENHADO)); cov2 = scale(log(df_pse$ORC_TOTAL_EMPENHADO))

cov1 = df_pnab$TX_PROP_MNST; cov2 = log(df_pse$TAM_PROP_BASE)

cor(cov1, cov2); cor(cov1, cov2, method ="kendall")

# Ajusta-se o modelo simples para as mesmas covariáveis do modelo conjunto.
# A ideia é usar estes parâmetros como chutes
ajust1 = survreg(data = df_pnab, Surv(time = t, event = ADOPT) ~
                 scale(log(PIB)) + scale(log(ORC_TOTAL_EMPENHADO)) #+
                 #scale(DISCR_TX) + scale(log(SERVIDORES))
                 ,
        dist = 'weibull')
summary(ajust1); c(AIC(ajust1), AICc(ajust1), BIC(ajust1))
vif(ajust1); check_multicol_survreg(ajust1)

ajust2 = survreg(data = df_pse, Surv(time = t, event = ADOPT) ~
                   scale(log(PIB)) + scale(TX_VETO) +
                   scale(DISCR_TX) #+ scale(log(DISCR_TX)) +
                   #scale(log(INDICE_MNST))
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
                #b31 * scale(df_pnab$DISCR_TX) + b41 * scale(log(df_pnab$SERVIDORES)) )
  # alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) +
  #               b31 * scale(df_pnab$DISCR_TX) + b41 * scale(log(df_pnab$SERVIDORES)) )
  #               #b51 * df_pnab$ANO_ELEITORAL_MUNIC)
  # alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) )
  # alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) )
  alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(df_pse$TX_VETO) +
                b32 * scale(df_pse$DISCR_TX) ) #+ b42 * scale(log(df_pse$DISCR_TX)) +
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

valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, .15)
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

# p-valores
vero<-MLE_nlminb$objective; parametro<-MLE_nlminb$par        
hessi <- numDeriv::hessian(func = vero_cop_debug, x = MLE_nlminb$par)
invR<-solve(hessi); variancia<-diag(invR); epp<-sqrt(variancia); estz=parametro/epp
( pvalor<-2*(1-pnorm(abs(estz))) )

# Em vez de solve(), usar pseudoinversa
library(MASS)
invR <- ginv(hessi); variancia <- diag(invR); epp <- sqrt(abs(variancia)); estz=parametro/epp

# gradiente
gradiente <- calc_gradiente(MLE_nlminb$par)
cat("Gradiente no ótimo:\n")
print(round(gradiente, 6))
cat("Norma do gradiente:", norm(gradiente, "2"), "\n")


## comparação com surv.copula ----
cov1 <- matrix(c(scale(log(df_pnab$PIB)), scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) ),
                 ncol = 2)
cov2 <- matrix(c(scale(log(df_pse$PIB)), scale(df_pse$TX_VETO),
                 scale(df_pse$DISCR_TX) ), ncol = 3)

library(Copula.surv)
(mod_clayton = Weib.reg.Clayton(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_frank = Weib.reg.Frank(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_gumbel = Weib.reg.Gumbel(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_bb1 = Weib.reg.BB1(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))


# pacote Copula.surv ----
library(Copula.surv)

t1 = df_pnab$t; delta1 = df_pnab$ADOPT
t2 = df_pse$t; delta2 = df_pse$ADOPT

## Uma covariável ----
cov1 = scale(log(df_pnab$ORC_TOTAL_EMPENHADO)); cov2 = scale(log(df_pse$ORC_TOTAL_EMPENHADO))
#cov1 = scale(df_pnab$TX_VETO); cov2 = scale(df_pse$TX_VETO)

## Mais covs
cov1 <- matrix(c(scale(log(df_pnab$PIB)), scale(log(df_pnab$ORC_TOTAL_EMPENHADO)),
                 scale(df_pnab$DISCR_TX)), ncol = 3)
cov2 <- matrix(c(scale(log(df_pse$PIB)), scale(df_pse$TX_VETO),
                 scale(df_pse$DISCR_TX) ), ncol = 3)

# colnames(cov1) <- c("a", "b")
# print(cov1)

(mod_clayton = Weib.reg.Clayton(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_frank = Weib.reg.Frank(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_gumbel = Weib.reg.Gumbel(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_bb1 = Weib.reg.BB1(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))

mod_clayton[["convergence"]][["ML"]]
(AIC_clayton = -2*mod_clayton[["convergence"]][["ML"]] + 2*mod_clayton[["convergence"]][["DF"]] )



## Modelo 1 ----

ajust1 = survreg(data = df_pnab, Surv(time = t, event = ADOPT) ~
                   scale(log(PIB)) + scale(log(ORC_TOTAL_EMPENHADO)) +
                   scale(DISCR_TX)
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
    b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; b41 = para[5]; gama1 = para[6];
    b02 = para[7]; b12 = para[8]; b22 = para[9]; b32 = para[10]; b42 = para[11]; b52 = para[12]; gama2 = para[13];
    theta = para[14]
  }
  
  if (gama1 <= 0 | gama2 <= 0) {
    cat("Gama não positivo:", gama1, gama2, "\n")
    return(1e10)
  }
  
  alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) +
                b31 * scale(df_pnab$DISCR_TX) )
  alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(df_pse$TX_VETO) +
                b32 * scale(log(df_pse$DISCR_DOTACAO)) )
  if ( any(!is.finite(alfa1)) || any(alfa1 <= 0) || any(!is.finite(alfa2)) || any(alfa2 <= 0) ) {
    cat("Alfa inválido\n")
    return(1e10)
  }
  
  log_vero_S = log_S_cop(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t1 = log_parc_t1(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t2 = log_parc_t2(alfa1, gama1, alfa2, gama2, theta)
  log_vero_f = log_f_cop(alfa1, gama1, alfa2, gama2, theta)
  
  result <- -sum(log_vero_S, log_vero_parc_t1, log_vero_parc_t2, log_vero_f)
  
  cat("Resultado:", ifelse(is.finite(result), round(result, 4), "Não finito"), "\n\n")
  
  if (!is.finite(result)) {
    return(1e10)
  }
  return(result)
}

valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, .8)
# nlminb é geralmente mais robusto que optim para problemas com bounds
MLE_nlminb <- nlminb(start = valorin1,
                     #lower = c(rep(-Inf, 10), 0.9028),
                     objective = vero_cop_debug,
                     control = list(iter.max = 2000, eval.max = 2000))
tau = MLE_nlminb[["par"]][[length(valorin1)]] / (2+MLE_nlminb[["par"]][[length(valorin1)]])
log_lik = -MLE_nlminb$objective; k = length(MLE_nlminb$par); n_observacoes = nrow(df_pnab)
aic = 2 * k - 2 * log_lik; aicc = aic + (2 * k * (k + 1)) / (n_observacoes - k - 1); bic = log(n_observacoes) * k - 2 * log_lik
print(MLE_nlminb); paste0('tau: ', round(tau, 4) )
paste0("AIC: ", round(aic, 3), "; AICc: ", round(aicc, 3), "; BIC: ", round(bic, 3), "; Vero: ", round(log_lik, 3) )

# p-valores
vero<-MLE_nlminb$objective; parametro<-MLE_nlminb$par        
hessi <- numDeriv::hessian(func = vero_cop_debug, x = MLE_nlminb$par)
invR<-solve(hessi); variancia<-diag(invR); epp<-sqrt(variancia); estz=parametro/epp
gsub('\\.', ',', round(epp, 4))
( pvalor<-2*(1-pnorm(abs(estz))) )

pvalor_formatado <- format(pvalor, scientific = TRUE, digits = )
print(pvalor_formatado)


# Em vez de solve(), usar pseudoinversa
library(MASS)
invR <- ginv(hessi); variancia <- diag(invR); epp <- sqrt(abs(variancia)); estz=parametro/epp

# gradiente
gradiente <- calc_gradiente(MLE_nlminb$par)
cat("Gradiente no ótimo:\n")
print(round(gradiente, 6))
cat("Norma do gradiente:", norm(gradiente, "2"), "\n")

## surv.copulas
cov1 <- matrix(c(scale(log(df_pnab$PIB)), scale(log(df_pnab$ORC_TOTAL_EMPENHADO)),
                 scale(df_pnab$DISCR_TX) ), ncol = 3)
cov2 <- matrix(c(scale(log(df_pse$PIB)), scale(df_pse$TX_VETO),
                 scale(log(df_pse$DISCR_DOTACAO)) ), ncol = 3)
library(Copula.surv)
(mod_clayton = Weib.reg.Clayton(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_frank = Weib.reg.Frank(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_gumbel = Weib.reg.Gumbel(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_bb1 = Weib.reg.BB1(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))



## Modelo 2 ----

ajust1 = survreg(data = df_pnab, Surv(time = t, event = ADOPT) ~
                   scale(log(PIB)) + scale(log(ORC_TOTAL_EMPENHADO)) +
                   scale(DISCR_TX)
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
    b01 = para[1]; b11 = para[2]; b21 = para[3]; b31 = para[4]; b41 = para[5]; gama1 = para[6];
    b02 = para[7]; b12 = para[8]; b22 = para[9]; b32 = para[10]; b42 = para[11]; b52 = para[12]; gama2 = para[13];
    theta = para[14]
  }
  
  if (gama1 <= 0 | gama2 <= 0) {
    cat("Gama não positivo:", gama1, gama2, "\n")
    return(1e10)
  }
  
  alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) +
                b31 * scale(df_pnab$DISCR_TX) )
  alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(df_pse$TX_VETO) +
                b32 * scale(df_pse$DISCR_TX) )
  if ( any(!is.finite(alfa1)) || any(alfa1 <= 0) || any(!is.finite(alfa2)) || any(alfa2 <= 0) ) {
    cat("Alfa inválido\n")
    return(1e10)
  }
  
  log_vero_S = log_S_cop(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t1 = log_parc_t1(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t2 = log_parc_t2(alfa1, gama1, alfa2, gama2, theta)
  log_vero_f = log_f_cop(alfa1, gama1, alfa2, gama2, theta)
  
  result <- -sum(log_vero_S, log_vero_parc_t1, log_vero_parc_t2, log_vero_f)
  
  cat("Resultado:", ifelse(is.finite(result), round(result, 4), "Não finito"), "\n\n")
  
  if (!is.finite(result)) {
    return(1e10)
  }
  return(result)
}

valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, 0.2)
# nlminb é geralmente mais robusto que optim para problemas com bounds
MLE_nlminb <- nlminb(start = valorin1,
                     #lower = c(rep(-Inf, 10), 0.9028),
                     objective = vero_cop_debug,
                     control = list(iter.max = 2000, eval.max = 2000))
tau = MLE_nlminb[["par"]][[length(valorin1)]] / (2+MLE_nlminb[["par"]][[length(valorin1)]])
log_lik = -MLE_nlminb$objective; k = length(MLE_nlminb$par); n_observacoes = nrow(df_pnab)
aic = 2 * k - 2 * log_lik; aicc = aic + (2 * k * (k + 1)) / (n_observacoes - k - 1); bic = log(n_observacoes) * k - 2 * log_lik
print(MLE_nlminb); paste0('tau: ', round(tau, 4) )
paste0("AIC: ", round(aic, 3), "; AICc: ", round(aicc, 3), "; BIC: ", round(bic, 3), "; Vero: ", round(log_lik, 3) )

# p-valores
vero<-MLE_nlminb$objective; parametro<-MLE_nlminb$par        
hessi <- numDeriv::hessian(func = vero_cop_debug, x = MLE_nlminb$par)
invR<-solve(hessi); variancia<-diag(invR); epp<-sqrt(variancia); estz=parametro/epp
( pvalor<-2*(1-pnorm(abs(estz))) )

gsub('\\.', ',', round(MLE_nlminb$par, 4))
gsub('\\.', ',', round(epp, 4))
gsub('\\.', ',', round(pvalor, 4))


# Em vez de solve(), usar pseudoinversa
library(MASS)
invR <- ginv(hessi); variancia <- diag(invR); epp <- sqrt(abs(variancia)); estz=parametro/epp

# gradiente
gradiente <- calc_gradiente(MLE_nlminb$par)
cat("Gradiente no ótimo:\n")
print(round(gradiente, 6))
cat("Norma do gradiente:", norm(gradiente, "2"), "\n")

## surv.copulas
cov1 <- matrix(c(scale(log(df_pnab$PIB)), scale(log(df_pnab$ORC_TOTAL_EMPENHADO)),
                 scale(df_pnab$DISCR_TX) ), ncol = 3)
cov2 <- matrix(c(scale(log(df_pse$PIB)), scale(df_pse$TX_VETO),
                 scale(df_pse$DISCR_TX) ), ncol = 3)
library(Copula.surv)
(mod_clayton = Weib.reg.Clayton(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_frank = Weib.reg.Frank(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_gumbel = Weib.reg.Gumbel(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_bb1 = Weib.reg.BB1(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))



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
  
  alfa1 = exp(b01 + b11 * scale(log(df_pnab$PIB)) + b21 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) ) 
  alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(df_pse$TX_VETO) +
                b32 * scale(log(df_pse$DISCR_DOTACAO)) ) 
  
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

valorin1<-c(ajust1$coefficients, 1/ajust1$scale, ajust2$coefficients, 1/ajust2$scale, .2)
# nlminb é geralmente mais robusto que optim para problemas com bounds
MLE_nlminb <- nlminb(start = valorin1,
                     #lower = c(rep(-Inf, 10), 0.9028),
                     objective = vero_cop_debug,
                     control = list(iter.max = 2000, eval.max = 2000))
tau = MLE_nlminb[["par"]][[length(valorin1)]] / (2+MLE_nlminb[["par"]][[length(valorin1)]])
log_lik = -MLE_nlminb$objective; k = length(MLE_nlminb$par); n_observacoes = nrow(df_pnab)
aic = 2 * k - 2 * log_lik; aicc = aic + (2 * k * (k + 1)) / (n_observacoes - k - 1); bic = log(n_observacoes) * k - 2 * log_lik
print(MLE_nlminb); paste0('tau: ', round(tau, 4) )
paste0("AIC: ", round(aic, 3), "; AICc: ", round(aicc, 3), "; BIC: ", round(bic, 3), "; Vero: ", round(log_lik, 3) )

# p-valores
vero<-MLE_nlminb$objective; parametro<-MLE_nlminb$par        
hessi <- numDeriv::hessian(func = vero_cop_debug, x = MLE_nlminb$par)
invR<-solve(hessi); variancia<-diag(invR); epp<-sqrt(variancia); estz=parametro/epp
( pvalor<-2*(1-pnorm(abs(estz))) )

gsub('\\.', ',', round(MLE_nlminb$par, 4))
gsub('\\.', ',', round(epp, 4))
gsub('\\.', ',', round(pvalor, 4))

# Em vez de solve(), usar pseudoinversa
library(MASS)
invR <- ginv(hessi); variancia <- diag(invR); epp <- sqrt(abs(variancia)); estz=parametro/epp

# gradiente
gradiente <- calc_gradiente(MLE_nlminb$par)
cat("Gradiente no ótimo:\n")
print(round(gradiente, 6))
cat("Norma do gradiente:", norm(gradiente, "2"), "\n")

## surv.copulas
cov1 <- matrix(c(scale(log(df_pnab$PIB)), scale(log(df_pnab$ORC_TOTAL_EMPENHADO)) ), 
               ncol = 2)
cov2 <- matrix(c(scale(log(df_pse$PIB)), scale(df_pse$TX_VETO) ),
                ncol = 2)
library(Copula.surv)
(mod_clayton = Weib.reg.Clayton(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_frank = Weib.reg.Frank(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_gumbel = Weib.reg.Gumbel(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_bb1 = Weib.reg.BB1(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))



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
                   scale(log(INDICE_MNST))                 ,
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
  alfa2 = exp(b02 + b12 * scale(log(df_pse$PIB)) + b22 * scale(df_pse$ORC_TOTAL_EMPENHADO) +
                b32 * scale(log(df_pse$TX_SUCESSO)) + b42 * scale(log(df_pse$DISCR_TX)) +
                b52 * scale(log(df_pse$INDICE_MNST)) )
  if ( any(!is.finite(alfa1)) || any(alfa1 <= 0) || any(!is.finite(alfa2)) || any(alfa2 <= 0) ) {
    cat("Alfa inválido\n")
    return(1e10)
  }
  
  log_vero_S = log_S_cop(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t1 = log_parc_t1(alfa1, gama1, alfa2, gama2, theta)
  log_vero_parc_t2 = log_parc_t2(alfa1, gama1, alfa2, gama2, theta)
  log_vero_f = log_f_cop(alfa1, gama1, alfa2, gama2, theta)
  
  result <- -sum(log_vero_S, log_vero_parc_t1, log_vero_parc_t2, log_vero_f)
  
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
                     control = list(iter.max = 2000, eval.max = 2000))
tau = MLE_nlminb[["par"]][[length(valorin1)]] / (2+MLE_nlminb[["par"]][[length(valorin1)]])
log_lik = -MLE_nlminb$objective; k = length(MLE_nlminb$par); n_observacoes = nrow(df_pnab)
aic = 2 * k - 2 * log_lik; aicc = aic + (2 * k * (k + 1)) / (n_observacoes - k - 1); bic = log(n_observacoes) * k - 2 * log_lik
print(MLE_nlminb); paste0('tau: ', round(tau, 4) )
paste0("AIC: ", round(aic, 3), "; AICc: ", round(aicc, 3), "; BIC: ", round(bic, 3), "; Vero: ", round(log_lik, 3) )



# p-valores
vero<-MLE_nlminb$objective; parametro<-MLE_nlminb$par        
hessi <- numDeriv::hessian(func = vero_cop_debug, x = MLE_nlminb$par)
invR<-solve(hessi); variancia<-diag(invR); epp<-sqrt(variancia); estz=parametro/epp
( pvalor<-2*(1-pnorm(abs(estz))) )

# Em vez de solve(), usar pseudoinversa
library(MASS)
invR <- ginv(hessi); variancia <- diag(invR); epp <- sqrt(abs(variancia)); estz=parametro/epp

gsub('\\.', ',', round(MLE_nlminb$par, 4))
gsub('\\.', ',', round(epp, 4))
gsub('\\.', ',', round(pvalor, 4))

# gradiente
gradiente <- calc_gradiente(MLE_nlminb$par)
cat("Gradiente no ótimo:\n")
print(round(gradiente, 6))
cat("Norma do gradiente:", norm(gradiente, "2"), "\n")

## surv.copulas
cov1 <- matrix(c(scale(log(df_pnab$PIB)), scale(log(df_pnab$ORC_TOTAL_EMPENHADO)),
                 scale(df_pnab$DISCR_TX), scale(log(df_pnab$SERVIDORES)) ), ncol = 4)
cov2 <- matrix(c(scale(log(df_pse$PIB)), scale(df_pse$ORC_TOTAL_EMPENHADO),
                 scale(log(df_pse$TX_SUCESSO)), scale(log(df_pse$DISCR_TX)),
                 scale(log(df_pse$INDICE_MNST)) ), ncol = 5)
library(Copula.surv)
(mod_clayton = Weib.reg.Clayton(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_frank = Weib.reg.Frank(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_gumbel = Weib.reg.Gumbel(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))
(mod_bb1 = Weib.reg.BB1(t1, t2, delta1, delta2, cov1, cov2, convergence.par=F))




# Estimativa pontual - mod log(orc) com log(orc) ----

vero_cop_debug = function(para) {
  cat("Parâmetros:", round(para, 6), "\n")
  
  b01 = para[1]; b11 = para[2]; gama1 = para[3]; b02 = para[4]; b12 = para[5]; gama2 = para[6]; theta = para[7]
  #alfa1 = para[1]; gama1 = para[2]; alfa2 = para[3]; gama2 = para[4]; theta = para[5]
  
  if (gama1 <= 0 | gama2 <= 0) {
    cat("Gama não positivo:", gama1, gama2, "\n")
    return(1e10)
  }
  
  #alfa1 = exp(b01 + b11 * log(df_pnab$ORC_TOTAL_EMPENHADO) )
  alfa1 = exp(b01 + b11 * scale(log(df_pnab$ORC_TOTAL_EMPENHADO) ) )
  #alfa1 = exp(b01 + b11 * scale(df_pnab$TX_VETO) )
  #alfa2 = exp(b02 + b12 * scale(df_pse$TX_VETO) )
  alfa2 = exp(b02 + b12 * scale(log(df_pse$ORC_TOTAL_EMPENHADO) ) )
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


# nlminb é geralmente mais robusto que optim para problemas com bounds
MLE_nlminb <- nlminb(start = #c(b01, b11, gama1, b02, b12, gama2, theta),
                       #c(-147.366100, 6.073400, 2.415459, 7.706900, 0.954600, 3.546099, 0.847800),
                       c(7, 1, 2.415459, 7.706900, 0.954600, 3.546099, 0.85),
                     #c(2100, 1, 2100, 1, .8), # sem covs 
                     objective = vero_cop_debug,
                     #lower = c(-15, 0, -0.005, -1, 0.0000001, 1, .2),
                     # upper = c(100, 10, 100, 100, 10, 100, 4, 2),
                     control = list(iter.max = 200, eval.max = 2000))

print(MLE_nlminb)


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


S_cop_estim(500, 500, median(scale(log(df_pnab$ORC_TOTAL_EMPENHADO))),
            median(scale(log(df_pse$ORC_TOTAL_EMPENHADO))), MLE_nlminb$par[1],
            MLE_nlminb$par[2], MLE_nlminb$par[3], MLE_nlminb$par[4],
            MLE_nlminb$par[5], MLE_nlminb$par[6], MLE_nlminb$par[7])

f_cop_estim(500, 500, median(scale(log(df_pnab$ORC_TOTAL_EMPENHADO))),
            median(scale(log(df_pse$ORC_TOTAL_EMPENHADO))), MLE_nlminb$par[1],
            MLE_nlminb$par[2], MLE_nlminb$par[3], MLE_nlminb$par[4],
            MLE_nlminb$par[5], MLE_nlminb$par[6], MLE_nlminb$par[7])

h_cop_estim(500, 500, median(scale(log(df_pnab$ORC_TOTAL_EMPENHADO))),
            median(scale(log(df_pse$ORC_TOTAL_EMPENHADO))), MLE_nlminb$par[1],
            MLE_nlminb$par[2], MLE_nlminb$par[3], MLE_nlminb$par[4],
            MLE_nlminb$par[5], MLE_nlminb$par[6], MLE_nlminb$par[7])

# Comparação modelo bivariado com univariado
S_cop_estim(500, 3000, median(scale(log(df_pnab$ORC_TOTAL_EMPENHADO))),
            median(scale(log(df_pse$TX_VETO))), MLE_nlminb$par[1],
            MLE_nlminb$par[2], MLE_nlminb$par[3], MLE_nlminb$par[4],
            MLE_nlminb$par[5], MLE_nlminb$par[6], MLE_nlminb$par[7])
S_estim(500, median(scale(log(df_pnab$ORC_TOTAL_EMPENHADO))),
            MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3]) *
  S_estim(3000, median(scale(log(df_pse$TX_VETO))),
          MLE_nlminb$par[4], MLE_nlminb$par[5], MLE_nlminb$par[6])



# Gráficos ----

## S ----
# Definir a função f(x,y) com os parâmetros fixos
f <- function(x, y) {
  S_cop_estim(x, y, # 0, 0,
              median(scale(log(df_pnab$ORC_TOTAL_EMPENHADO))),
              median(scale(log(df_pse$ORC_TOTAL_EMPENHADO))), 
              MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
              MLE_nlminb$par[4], MLE_nlminb$par[5], MLE_nlminb$par[6], 
              MLE_nlminb$par[7])
}
# Criar sequências de valores de 0 a 6000
x <- seq(0, 1100, length.out = 30)
y <- seq(0, 6000, length.out = 30)
# Calcular a matriz z
z <- outer(x, y, f)


library(lattice)
# Preparar dados
df_surface <- expand.grid(x = x, y = y)
df_surface$z <- as.vector(z)
wireframe(z ~ x * y, data = df_surface,
          main = expression(paste("Superfície - S"[cop], "(x,y)")),
          xlab = "x", ylab = "y", zlab = "S(x,y)",
          drape = TRUE,  # Adiciona cores
          colorkey = TRUE,
          screen = list(z = -60, x = -60))


library(plotly)
plot_ly(x = ~x, y = ~y, z = ~t(z), type = "surface") %>%
  layout(
    title = expression(paste("Superfície 3D Interativa - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "S(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )


## f ----
# Definir a função f(x,y) com os parâmetros fixos
f <- function(x, y) {
  f_cop_estim(x, y, # 0, 0,
              median(scale(log(df_pnab$ORC_TOTAL_EMPENHADO))),
              median(scale(log(df_pse$ORC_TOTAL_EMPENHADO))), 
              MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
              MLE_nlminb$par[4], MLE_nlminb$par[5], MLE_nlminb$par[6], 
              MLE_nlminb$par[7])
}
# Criar sequências de valores de 0 a 6000
x <- seq(0, 1300, length.out = 60)
y <- seq(0, 7000, length.out = 60)
# Calcular a matriz z
z <- outer(x, y, f)


library(lattice)
# Preparar dados
df_surface <- expand.grid(x = x, y = y)
df_surface$z <- as.vector(z)
wireframe(z ~ x * y, data = df_surface,
          main = expression(paste("Superfície - S"[cop], "(x,y)")),
          xlab = "x", ylab = "y", zlab = "f(x,y)",
          drape = TRUE,  # Adiciona cores
          colorkey = TRUE,
          screen = list(z = -60, x = -60))


library(plotly)
plot_ly(x = ~x, y = ~y, z = ~t(z), type = "surface") %>%
  layout(
    title = expression(paste("Superfície 3D Interativa - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "f(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )


library(ggplot2)
library(reshape2)
# Converter dados para formato longo (necessário para ggplot2)
df <- melt(z)
names(df) <- c("x_index", "y_index", "value")
df$x <- x[df$x_index]
df$y <- y[df$y_index]
# Gráfico de contorno com ggplot2
ggplot(df, aes(x = x, y = y, z = value)) +
  geom_contour_filled() +
  labs(title = "Gráfico de Contorno - f(x,y)",
       x = "X", y = "Y", fill = "f(x,y)") +
  theme_minimal()


## h ----
# Definir a função h(x,y) com os parâmetros fixos
f <- function(x, y) {
  h_cop_estim(x, y, # 0, 0,
              median(scale(log(df_pnab$ORC_TOTAL_EMPENHADO))),
              median(scale(log(df_pse$ORC_TOTAL_EMPENHADO))),
              MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
              MLE_nlminb$par[4], MLE_nlminb$par[5], MLE_nlminb$par[6], 
              MLE_nlminb$par[7])
}
# Criar sequências de valores de 0 a 6000
x <- seq(0, 1300, length.out = 60)
y <- seq(0, 7000, length.out = 60)
# x <- seq(3500, 12000, length.out = 120)
# y <- seq(12000, 25000, length.out = 120)
# Calcular a matriz z
z <- outer(x, y, f)


library(lattice)
# Preparar dados
df_surface <- expand.grid(x = x, y = y)
df_surface$z <- as.vector(z)
wireframe(z ~ x * y, data = df_surface,
          main = expression(paste("Superfície - h"[cop], "(x,y)")),
          xlab = "x", ylab = "y", zlab = "h(x,y)",
          drape = TRUE,  # Adiciona cores
          colorkey = TRUE,
          screen = list(z = -60, x = -60))


library(plotly)
plot_ly(x = ~x, y = ~y, z = ~t(z), type = "surface") %>%
  layout(
    title = expression(paste("Superfície 3D Interativa - h"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "f(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )

## Duas superfícies para comparação - teste ----

# Definir a função f(x,y) com os parâmetros fixos
f <- function(x, y) {
  S_cop_estim(x, y, # 0, 0,
              median(scale(log(df_pnab$ORC_TOTAL_EMPENHADO))),
              median(scale(log(df_pse$ORC_TOTAL_EMPENHADO))), 
              MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
              MLE_nlminb$par[4], MLE_nlminb$par[5], MLE_nlminb$par[6], 
              MLE_nlminb$par[7])
}

f2 <- function(x, y) {
  S_cop_estim(x, y, 0, 0,
              # median(scale(log(df_pnab$ORC_TOTAL_EMPENHADO))),
              # median(scale(log(df_pse$ORC_TOTAL_EMPENHADO))), 
              MLE_nlminb$par[1], MLE_nlminb$par[2], MLE_nlminb$par[3], 
              MLE_nlminb$par[4], MLE_nlminb$par[5], MLE_nlminb$par[6], 
              MLE_nlminb$par[7])
}

# Criar sequências de valores de 0 a 6000
x <- seq(0, 2000, length.out = 60)
y <- seq(0, 6000, length.out = 60)
# Calcular a matriz z
z <- outer(x, y, f)
z2 <- outer(x, y, f2)


library(plotly)
plot_ly() %>%
  add_surface(x = ~x, y = ~y, z = ~t(z), 
              colorscale = list(c(0, 1), c("#1a074a", "#caf4fc")),
              name = "Superfície 1") %>%
  add_surface(x = ~x, y = ~y, z = ~t(z2), 
              colorscale = list(c(0, .7, 1), c("#871313", "#f7a516", "#fcf477")),
              name = "Superfície 2") %>%
  layout(
    title = expression(paste("Superfícies 3D - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "S(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )

z_dif = z-z2
z_zero <- matrix(0, nrow = length(x), ncol = length(y))

plot_ly() %>%
  add_surface(x = ~x, y = ~y, z = ~t(z_dif), 
              colorscale = "Blues",
              opacity = 0.8,
              name = "S(x,y)") %>%
  add_surface(x = ~x, y = ~y, z = ~z_zero, 
              opacity = 0.4,
              colorscale = list(c(0, 1), c("gray", "gray")),
              showscale = FALSE,
              name = "Plano z = 0") %>%
  layout(
    title = expression(paste("Superfície 3D com Plano de Referência - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "S(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )




# Duas superfícies para comparação - concreto ----

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


library(plotly)
plot_ly() %>%
  add_surface(x = ~x, y = ~y, z = ~t(z), 
              colorscale = list(c(0, 1), c("#1a074a", "#caf4fc")),
              name = "Superfície 1") %>%
  add_surface(x = ~x, y = ~y, z = ~t(z2), 
              colorscale = list(c(0, .7, 1), c("#871313", "#f7a516", "#fcf477")),
              name = "Superfície 2") %>%
  layout(
    title = expression(paste("Superfícies 3D - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "S(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )

z_dif = z-z2
z_zero <- matrix(0, nrow = length(x), ncol = length(y))

plot_ly(x = ~x, y = ~y, z = ~t(z_dif), type = "surface") %>%
  add_surface(x = ~x, y = ~y, z = ~z_zero, 
              opacity = 0.4,
              colorscale = list(c(0, 1), c("gray", "gray")),
              showscale = FALSE,
              name = "Plano z = 0") %>%
  layout(
    title = expression(paste("Superfície 3D com Plano de Referência - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "S(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )


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


pdf("imagens/plano_mod1.pdf", 
    width = 10, height = 8,  # Polegadas
    pointsize = 12)

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

dev.off()


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


library(plotly)
plot_ly() %>%
  add_surface(x = ~x, y = ~y, z = ~t(z), 
              colorscale = list(c(0, 1), c("#1a074a", "#caf4fc")),
              name = "Superfície 1") %>%
  add_surface(x = ~x, y = ~y, z = ~t(z2), 
              colorscale = list(c(0, .7, 1), c("#871313", "#f7a516", "#fcf477")),
              name = "Superfície 2") %>%
  layout(
    title = expression(paste("Superfícies 3D - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "S(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )

z_dif = z-z2
z_zero <- matrix(0, nrow = length(x), ncol = length(y))


plot_ly(x = ~x, y = ~y, z = ~t(z_dif), type = "surface") %>%
  add_surface(x = ~x, y = ~y, z = ~z_zero, 
              opacity = 0.4,
              colorscale = list(c(0, 1), c("gray", "gray")),
              showscale = FALSE,
              name = "Plano z = 0") %>%
  layout(
    title = expression(paste("Superfície 3D com Plano de Referência - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "S(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )

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


pdf("imagens/plano_mod2.pdf", 
    width = 10, height = 8,  # Polegadas
    pointsize = 12)

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

dev.off()



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


library(plotly)
plot_ly() %>%
  add_surface(x = ~x, y = ~y, z = ~t(z), 
              colorscale = list(c(0, 1), c("#1a074a", "#caf4fc")),
              name = "Superfície 1") %>%
  add_surface(x = ~x, y = ~y, z = ~t(z2), 
              colorscale = list(c(0, .7, 1), c("#871313", "#f7a516", "#fcf477")),
              name = "Superfície 2") %>%
  layout(
    title = expression(paste("Superfícies 3D - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "S(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )

z_dif = z-z2
z_zero <- matrix(0, nrow = length(x), ncol = length(y))


plot_ly(x = ~x, y = ~y, z = ~t(z_dif), type = "surface") %>%
  add_surface(x = ~x, y = ~y, z = ~z_zero, 
              opacity = 0.4,
              colorscale = list(c(0, 1), c("gray", "gray")),
              showscale = FALSE,
              name = "Plano z = 0") %>%
  layout(
    title = expression(paste("Superfície 3D com Plano de Referência - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "S(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )

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


pdf("imagens/plano_mod3.pdf", 
    width = 10, height = 8,  # Polegadas
    pointsize = 12)

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

dev.off()



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


library(plotly)
plot_ly() %>%
  add_surface(x = ~x, y = ~y, z = ~t(z), 
              colorscale = list(c(0, 1), c("#1a074a", "#caf4fc")),
              name = "Superfície 1") %>%
  # add_surface(x = ~x, y = ~y, z = ~t(z2), 
  #             colorscale = list(c(0, .7, 1), c("#871313", "#f7a516", "#fcf477")),
  #             name = "Superfície 2") %>%
  layout(
    title = expression(paste("Superfícies 3D - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "S(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )

z_dif = z-z2
z_zero <- matrix(0, nrow = length(x), ncol = length(y))


plot_ly(x = ~x, y = ~y, z = ~t(z_dif), type = "surface") %>%
  # add_surface(x = ~x, y = ~y, z = ~z_zero, 
  #             opacity = 0.4,
  #             colorscale = list(c(0, 1), c("gray", "gray")),
  #             showscale = FALSE,
  #             name = "Plano z = 0") %>%
  layout(
    title = expression(paste("Superfície 3D com Plano de Referência - S"[cop], "(x,y)")),
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"), 
      zaxis = list(title = "S(x,y)"),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    )
  )

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
