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

caminho = dirname(rstudioapi::getActiveDocumentContext()$path)
base_completa = read_excel(paste0(caminho, '/', "BASE_FINAL.xlsx"), 
                           sheet = "BASE_COMPLETA")
base_mod = read_excel(paste0(caminho, '/', "BASE_FINAL.xlsx"), 
                      sheet = "BASE_MOD")
names(base_mod)[names(base_mod) == "T_MEDIO_PERMAN_MINISTRO_P_MANDATO"] <- "T_MEDIO_PERMAN_MINISTRO"

variaveis <- c("TX_VETO", "TX_CONFLITO_MP", "TX_PARTC", "TX_SUCESSO", "TX_DOMINANCIA", 
               "TX_FRAG_PARTD_NEP", "TAM_PROP_BASE", "TX_DISCIP_BASE", "TX_PROP_MNST", 
               "POP_PRSD", "T_MEDIO_PERMAN_MINISTRO", "CARREIRAS", 
               "NUM_ENTIDADES", "SERVIDORES", "SERVIDORES_RJU", "DISCR_DOTACAO", 
               "DISCR_EMPENHADO", "ORC_TOTAL_DOTACAO", "ORC_TOTAL_EMPENHADO", 
               "DISCR_TX", "ORC_TOTAL_TX", "TI_PROJ_LEI", "TI_EMPENHADO", 
               "INDICE_MNST", "PIB", "ANO_ELEITORAL_MUNIC", "ANO_ELEITORAL_NAC", 
               "MRG_PRSD", "IDEOL_PRSD", "IDEOL_PRF")

variaveis <- c("t", "ADOPT", "TX_VETO", "TX_SUCESSO", "SERVIDORES", "DISCR_DOTACAO", 
               "ORC_TOTAL_EMPENHADO", "DISCR_TX", "INDICE_MNST", "PIB")


# Geral - sem filtro pro política ----
# Selecionar apenas as variáveis de interesse do banco de dados
banco_selecionado <- base_mod[variaveis]

# PNAB ----
variaveis = variaveis[variaveis != "NUM_ENTIDADES"]

df_pnab = base_mod |> filter(POL_ABV == 'pnab')
#df_pnab$NUM_ENTIDADES = NULL

df_pnab = df_pnab[variaveis]

# Função para aplicar log com tratamento de valores não-positivos
aplicar_log_seguro <- function(x) {
  if (is.numeric(x)) {
    # Verificar se há valores não-positivos
    if (any(x <= 0, na.rm = TRUE)) {
      # Encontrar o menor valor positivo
      min_positivo <- min(x[x > 0], na.rm = TRUE)
      # Adicionar constante para tornar todos positivos
      constante <- ifelse(min_positivo <= 1, 1 - min(x, na.rm = TRUE) + 0.01, 0)
      return(log(x + constante))
    } else {
      return(log(x))
    }
  } else {
    return(x) # Retorna original se não for numérico
  }
}

# Identificar colunas numéricas (excluindo variáveis de tempo e evento)
colunas_numericas <- names(df_pnab)[sapply(df_pnab, is.numeric)]

# Remover colunas que não devem receber log (como tempo, evento, IDs, etc.)
colunas_para_excluir <- c("t", "ADOPT", "ID", "id", "codigo", "ano", "year") # ajuste conforme necessário
colunas_para_log <- setdiff(colunas_numericas, colunas_para_excluir)

# Aplicar log a todas as colunas numéricas selecionadas
for (coluna in colunas_para_log) {
  nova_coluna <- paste0("log(", coluna, ")")
  df_pnab[[nova_coluna]] <- aplicar_log_seguro(df_pnab[[coluna]])
  
  # Verificar se a transformação foi aplicada corretamente
  if (any(is.infinite(df_pnab[[nova_coluna]]), na.rm = TRUE)) {
    warning("Valores infinitos encontrados em ", nova_coluna)
  }
  if (any(is.na(df_pnab[[nova_coluna]]), na.rm = TRUE)) {
    warning("NAs encontrados em ", nova_coluna)
  }
}


vars = names(df_pnab)


# Criar dataframe para armazenar os critérios de informação
resumo_modelos_pnab <- data.frame(
  Variavel = character(),
  Beta = numeric(),
  Erro_Padrao = numeric(),  # NOVO CAMPO
  P_valor = numeric(),
  BIC = numeric(),
  stringsAsFactors = FALSE
)

# Ajustar modelos e armazenar resultados
for (variavel in vars) {
  tryCatch({
    # Verificar se a variável existe no dataframe
    if (!variavel %in% names(df_pnab)) {
      cat("Variável", variavel, "não encontrada no dataframe. Pulando.\n")
      next
    }
    
    formula_modelo <- as.formula(paste("Surv(time = t, event = ADOPT) ~ scale(", variavel, ")"))
    modelo <- survreg(data = df_pnab, formula = formula_modelo, dist = 'weibull')
    
    # Extrair sumário do modelo
    sumario <- summary(modelo)
    
    # Extrair tabela de coeficientes
    tabela_coef <- sumario$table
    
    # Encontrar a linha correspondente à variável (não ao intercepto)
    linha_variavel <- which(!rownames(tabela_coef) == "(Intercept)")[1]
    
    # Extrair o coeficiente da variável
    coeficiente <- coef(modelo)[linha_variavel]
    
    # EXTRAIR O ERRO PADRÃO DA VARIÁVEL
    erro_padrao <- tabela_coef[linha_variavel, 2]  # Segunda coluna é o Std. Error
    
    # Extrair o p-valor da variável
    p_valor <- tabela_coef[linha_variavel, 4]
    
    # Determinar o sinal do beta
    sinal <- ifelse(coeficiente > 0, "Positivo", "Negativo")
    
    # Adicionar ao resumo
    resumo_modelos_pnab <- rbind(resumo_modelos_pnab, data.frame(
      Variavel = variavel,
      Beta = round(coeficiente, 4),
      Erro_Padrao = round(erro_padrao, 4),  # NOVO CAMPO
      P_valor = round(p_valor, 4),
      BIC = BIC(modelo),
      stringsAsFactors = FALSE
    ))
    
    
  }, error = function(e) {
    cat("Modelo para", variavel, "não pôde ser ajustado:", e$message, "\n")
  })
}

# Ordenar por BIC
resumo_modelos_pnab <- resumo_modelos_pnab[order(resumo_modelos_pnab$BIC), ]

# Visualizar o resultado
print(resumo_modelos_pnab)

#write_xlsx(resumo_modelos_pnab, "resumo_modelos_pnab_final.xlsx")

## gera matriz de corr ----

#write_xlsx(df_pnab, "pnab_vars_mod.xlsx")

pnab_vars_mod_att <- read_excel("pnab_vars_mod_att.xlsx")

matriz_cor <- cor(pnab_vars_mod_att[,3:18], use = "complete.obs")
library(corrplot)

ordem_manual <- c("PIB per capita", "log(PIB per capita)",
                  "Orc total executado", "log(Orc total executado)",
                  "Taxa de veto", "log(Taxa de veto)", 
                  "Orc discr autorizado", "log(Orc discr autorizado)",
                  "Tx de exec do orc discr", "log(Tx de exec do orc discr)", 
                  "N de servidores", "log(N de servidores)",
                  "Taxa de aprov leg", "log(Taxa de aprov leg)", 
                  "Taxa ministerial", "log(Taxa ministerial)" )

matriz_cor <- matriz_cor[ordem_manual, ordem_manual]

# PDF com engine Cairo (melhor para fontes e gráficos vetoriais)
# cairo_pdf("imagens/matriz_correlacao_pnab.pdf", 
#           width = 10, 
#           height = 8)

corrplot(matriz_cor, 
         method = "color",
         type = "upper",
         order = "original",
         addCoef.col = "black",
         tl.col = "black",
         tl.cex = 0.8,
         number.cex = 0.7)

#dev.off()



## exportar tabela com maiores corrs ----

tabela_cor_sem_duplicatas <- as.data.frame(as.table(matriz_cor)) %>%
  rename(Var1 = 1, Var2 = 2, Correlacao = 3) %>%
  filter(
    abs(Correlacao) >= 0.7,
    as.character(Var1) != as.character(Var2)  # Remove correlação consigo mesma
  ) %>%
  # Remove duplicatas onde (Var1, Var2) = (Var2, Var1)
  rowwise() %>%
  mutate(par_ordenado = paste(sort(c(Var1, Var2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(par_ordenado, .keep_all = TRUE) %>%
  select(-par_ordenado) %>%
  # Ordenar e formatar
  arrange(desc(abs(Correlacao))) %>%
  mutate(Correlacao = round(Correlacao, 3),
         Tipo = ifelse(Correlacao > 0, "Positiva", "Negativa"))

#write_xlsx(tabela_cor_sem_duplicatas, "filtro_corr_fortes_pnab.xlsx")

library(readxl)
filtro_corr_fortes_pnab_final <- read_excel("filtro_corr_fortes_pnab_final.xlsx")

library(xtable)
# Converter para LaTeX
print(xtable(filtro_corr_fortes_pnab_final, caption = "Tabela exemplo"), 
      include.rownames = FALSE,
      caption.placement = "top")

## exporta tabela dos ajustes ----

library(xtable)
library(dplyr)
library(stringr)

# Carregar os dados
resumo_modelos_pnab_final2 <- read_excel("resumo_modelos_pnab_final2.xlsx")

# Função para formatar números com condições específicas
formatar_numero_br <- function(x, nome_coluna, decimais = 4) {
  if(is.numeric(x)) {
    # Condição especial para coluna BIC - apenas 1 casa decimal
    if(nome_coluna == "BIC") {
      decimais <- 1
    }
    
    # Arredondar
    x_arredondado <- round(x, decimais)
    
    # Verificar se o valor arredondado é zero e substituir
    if(abs(x_arredondado) < 0.0001) {
      return("$<$0,0001")
    } else {
      # Formatar com separadores brasileiros
      return(format(x_arredondado, 
                    nsmall = decimais, 
                    decimal.mark = ",", 
                    big.mark = ".", 
                    scientific = FALSE))
    }
  } else {
    return(x)
  }
}

# Aplicar formatação às colunas numéricas
resumo_formatado <- resumo_modelos_pnab_final2

# Aplicar a função para cada coluna numérica
for(coluna in names(resumo_formatado)) {
  if(is.numeric(resumo_formatado[[coluna]])) {
    resumo_formatado[[coluna]] <- sapply(resumo_formatado[[coluna]], 
                                         function(x) formatar_numero_br(x, coluna))
  }
}

# Converter para LaTeX
print(xtable(resumo_formatado, caption = "Resultados dos Modelos de Sobrevivência"), 
      include.rownames = FALSE,
      caption.placement = "top",
      sanitize.text.function = identity)





# PSE ----

variaveis = variaveis[variaveis != "NUM_ENTIDADES"]

df_pse = base_mod |> filter(POL_ABV == 'pse')
#df_pse$NUM_ENTIDADES = NULL

df_pse = df_pse[variaveis]

# Identificar colunas numéricas (excluindo variáveis de tempo e evento)
colunas_numericas <- names(df_pse)[sapply(df_pse, is.numeric)]

# Remover colunas que não devem receber log (como tempo, evento, IDs, etc.)
colunas_para_excluir <- c("t", "ADOPT", "ID", "id", "codigo", "ano", "year") # ajuste conforme necessário
colunas_para_log <- setdiff(colunas_numericas, colunas_para_excluir)

# Aplicar log a todas as colunas numéricas selecionadas
for (coluna in colunas_para_log) {
  nova_coluna <- paste0("log(", coluna, ")")
  df_pse[[nova_coluna]] <- aplicar_log_seguro(df_pse[[coluna]])
  
  # Verificar se a transformação foi aplicada corretamente
  if (any(is.infinite(df_pse[[nova_coluna]]), na.rm = TRUE)) {
    warning("Valores infinitos encontrados em ", nova_coluna)
  }
  if (any(is.na(df_pse[[nova_coluna]]), na.rm = TRUE)) {
    warning("NAs encontrados em ", nova_coluna)
  }
}


vars = names(df_pse)


# Criar dataframe para armazenar os critérios de informação
resumo_modelos_pse <- data.frame(
  Variavel = character(),
  Beta = numeric(),
  Erro_Padrao = numeric(),  # NOVO CAMPO
  P_valor = numeric(),
  BIC = numeric(),
  stringsAsFactors = FALSE
)

# Ajustar modelos e armazenar resultados
for (variavel in vars) {
  tryCatch({
    # Verificar se a variável existe no dataframe
    if (!variavel %in% names(df_pse)) {
      cat("Variável", variavel, "não encontrada no dataframe. Pulando.\n")
      next
    }
    
    formula_modelo <- as.formula(paste("Surv(time = t, event = ADOPT) ~ scale(", variavel, ")"))
    modelo <- survreg(data = df_pse, formula = formula_modelo, dist = 'weibull')
    
    # Extrair sumário do modelo
    sumario <- summary(modelo)
    
    # Extrair tabela de coeficientes
    tabela_coef <- sumario$table
    
    # Encontrar a linha correspondente à variável (não ao intercepto)
    linha_variavel <- which(!rownames(tabela_coef) == "(Intercept)")[1]
    
    # Extrair o coeficiente da variável
    coeficiente <- coef(modelo)[linha_variavel]
    
    # EXTRAIR O ERRO PADRÃO DA VARIÁVEL
    erro_padrao <- tabela_coef[linha_variavel, 2]  # Segunda coluna é o Std. Error
    
    # Extrair o p-valor da variável
    p_valor <- tabela_coef[linha_variavel, 4]
    
    # Determinar o sinal do beta
    sinal <- ifelse(coeficiente > 0, "Positivo", "Negativo")
    
    # Adicionar ao resumo
    resumo_modelos_pse <- rbind(resumo_modelos_pse, data.frame(
      Variavel = variavel,
      Beta = round(coeficiente, 4),
      Erro_Padrao = round(erro_padrao, 4),  # NOVO CAMPO
      P_valor = round(p_valor, 4),
      BIC = BIC(modelo),
      stringsAsFactors = FALSE
    ))
    
    
  }, error = function(e) {
    cat("Modelo para", variavel, "não pôde ser ajustado:", e$message, "\n")
  })
}

# Ordenar por BIC
resumo_modelos_pse <- resumo_modelos_pse[order(resumo_modelos_pse$BIC), ]

# Visualizar o resultado
print(resumo_modelos_pse)

#write_xlsx(resumo_modelos_pse, "resumo_modelos_pse_final.xlsx")


## gera matriz de corr ----

#write_xlsx(df_pse, "pse_vars_mod.xlsx")

pse_vars_mod_att <- read_excel("pse_vars_mod_att.xlsx")

matriz_cor <- cor(pse_vars_mod_att[,3:18], use = "complete.obs")
library(corrplot)

ordem_manual <- c("PIB per capita", "log(PIB per capita)",
                  "Orc total executado", "log(Orc total executado)",
                  "Taxa de veto", "log(Taxa de veto)", 
                  "Orc discr autorizado", "log(Orc discr autorizado)",
                  "Tx de exec do orc discr", "log(Tx de exec do orc discr)", 
                  "N de servidores", "log(N de servidores)",
                  "Taxa de aprov leg", "log(Taxa de aprov leg)", 
                  "Taxa ministerial", "log(Taxa ministerial)" )

matriz_cor <- matriz_cor[ordem_manual, ordem_manual]

# PDF com engine Cairo (melhor para fontes e gráficos vetoriais)
# cairo_pdf("imagens/matriz_correlacao_pse.pdf", 
#           width = 10, 
#           height = 8)

corrplot(matriz_cor, 
         method = "color",
         type = "upper",
         order = 'original',
         addCoef.col = "black",
         tl.col = "black",
         tl.cex = 0.8,
         number.cex = 0.7)

# dev.off()

tabela_cor_sem_duplicatas <- as.data.frame(as.table(matriz_cor)) %>%
  rename(Var1 = 1, Var2 = 2, Correlacao = 3) %>%
  filter(
    abs(Correlacao) >= 0.7,
    as.character(Var1) != as.character(Var2)  # Remove correlação consigo mesma
  ) %>%
  # Remove duplicatas onde (Var1, Var2) = (Var2, Var1)
  rowwise() %>%
  mutate(par_ordenado = paste(sort(c(Var1, Var2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(par_ordenado, .keep_all = TRUE) %>%
  select(-par_ordenado) %>%
  # Ordenar e formatar
  arrange(desc(abs(Correlacao))) %>%
  mutate(Correlacao = round(Correlacao, 3),
         Tipo = ifelse(Correlacao > 0, "Positiva", "Negativa"))

#write_xlsx(tabela_cor_sem_duplicatas, "filtro_corr_fortes_pse.xlsx")

library(readxl)
filtro_corr_fortes_pse_final <- read_excel("filtro_corr_fortes_pse_final.xlsx")

library(xtable)
# Converter para LaTeX
print(xtable(filtro_corr_fortes_pse_final, caption = "Tabela exemplo"), 
      include.rownames = FALSE,
      caption.placement = "top")



## exporta tabela dos ajustes ----

library(xtable)
library(dplyr)
library(stringr)

# Carregar os dados
resumo_modelos_pse_final2 <- read_excel("resumo_modelos_pse_final2.xlsx")

# Função para formatar números com condições específicas
formatar_numero_br <- function(x, nome_coluna, decimais = 4) {
  if(is.numeric(x)) {
    # Condição especial para coluna BIC - apenas 1 casa decimal
    if(nome_coluna == "BIC") {
      decimais <- 1
    }
    
    # Arredondar
    x_arredondado <- round(x, decimais)
    
    # Verificar se o valor arredondado é zero e substituir
    if(abs(x_arredondado) < 0.0001) {
      return("$<$0,0001")
    } else {
      # Formatar com separadores brasileiros
      return(format(x_arredondado, 
                    nsmall = decimais, 
                    decimal.mark = ",", 
                    big.mark = ".", 
                    scientific = FALSE))
    }
  } else {
    return(x)
  }
}

# Aplicar formatação às colunas numéricas
resumo_formatado <- resumo_modelos_pse_final2

# Aplicar a função para cada coluna numérica
for(coluna in names(resumo_formatado)) {
  if(is.numeric(resumo_formatado[[coluna]])) {
    resumo_formatado[[coluna]] <- sapply(resumo_formatado[[coluna]], 
                                         function(x) formatar_numero_br(x, coluna))
  }
}

# Converter para LaTeX
print(xtable(resumo_formatado, caption = "Resultados dos Modelos de Sobrevivência"), 
      include.rownames = FALSE,
      caption.placement = "top",
      sanitize.text.function = identity)
