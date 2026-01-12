library(ggplot2)

# Gráfico da Weibull ----
alfa <- 1
gama_values <- c(1, 2.5, 0.5)
colors <- c("#3E9E0F", "#FB1808", "steelblue")
labels <- c(
  expression(Weibull(alpha==1~","~gamma==1)),
  expression(Weibull(alpha==1~","~gamma==2.5)),
  expression(Weibull(alpha==1~","~gamma==0.5))
)

# Criar um data frame com as curvas
x <- seq(0, 0.8, length.out = 100)
df <- data.frame(x = rep(x, length(gama_values)),
                 gama = rep(gama_values, each = length(x)),
                 h = unlist(lapply(gama_values, function(g) (g/alfa^g) * x^(g-1))),
                 group = rep(as.character(labels), each = length(x)))

# Criar o gráfico
ggplot(df, aes(x = x, y = h, color = group)) +
  geom_line(size = 1) +
  scale_color_manual(values = setNames(colors, as.character(labels)),
                     labels = function(x) parse(text = x)) +
  scale_x_continuous(labels = number_format(decimal.mark = ",", big.mark = ".")) +
  scale_y_continuous(labels = number_format(decimal.mark = ",", big.mark = ".")) +
  labs(x = "t",
       y = "h(t)",
       color = "Distribuição") +
  ylim(0, 4) +
  theme_minimal() +
  theme(
    legend.position = c(0.75, 0.85),
    legend.text = element_text(size = 9),
    legend.title = element_blank(),
    legend.background = element_rect(
      color = "gray",
      fill = "white",
      linetype = ,
      size = 0.5
    )#,
    #legend.key.height = unit(0.8, "lines"),  # Ajusta altura dos itens
    #legend.margin = margin(0, 15, 0, 10)  # Aumenta margens internas
  )
#ggsave("imagens/ref_h_wei.pdf", width = 158, height = 93, units = "mm")


# Gráfico da log-normal ----
mu <- 0
sigma <- 1

# Criar um data frame com as curvas
x <- seq(0, 4, length.out = 100)
df <- data.frame(x = x, g='1',
                 h = dlnorm(x, mean = mu, sd = sigma) / # f(t)
                   pnorm((-log(x)+mu)/sigma, mean=0, sd=1) )

ggplot(df, aes(x = x, y = h, color = g)) +
  geom_line(size = 1) +  # Cor da curva alterada para vermelho
  labs(x = "t",
       y = "h(t)",
       color = "Legenda") +  # Título da legenda
  scale_x_continuous(labels = number_format(decimal.mark = ",", big.mark = ".")) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(decimal.mark = ",", big.mark = ".")) +
  #ylim(0, 1) +
  xlim(0, 4) +
  theme_minimal() +
  theme(
    legend.position = c(0.75, 0.85),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.background = element_rect(
      color = "gray",
      fill = "white",
      linetype = "solid",
      size = 0.5
    )#,
    #legend.margin = margin(0, 23, 0, 5)
  ) +
  scale_color_manual(values = "steelblue",  # Define a cor da legenda para corresponder à curva
                     labels = expression(log-normal(mu==0~","~sigma==1)))  # Texto da legenda
#ggsave("imagens/ref_h_lnorm.pdf", width = 158, height = 93, units = "mm")


# Gráfico da log-logística ----
alfa <- 1
gama <- 3

# Criar um data frame com as curvas
x <- seq(0, 4, length.out = 100)
df <- data.frame(x = x, g='1',
                 h = gama*(x/alfa)^(gama-1) / (alfa*(1+(x/alfa)^gama)) )

ggplot(df, aes(x = x, y = h, color = g)) +
  geom_line(size = 1) +  # Cor da curva alterada para vermelho
  labs(x = "t",
       y = "h(t)") +  # Título da legenda
  scale_y_continuous(limits = c(0, 2), labels = number_format(decimal.mark = ",", big.mark = ".")) +
  xlim(0, 4) +
  theme_minimal() +
  theme(
    legend.position = c(0.75, 0.85),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.background = element_rect(
      color = "gray",
      fill = "white",
      linetype = "solid",
      size = 0.5
    )#,
    #legend.margin = margin(0, 23, 0, 5)
  ) +
  scale_color_manual(values = "steelblue",  # Define a cor da legenda para corresponder à curva
                     labels = 
                     expression(log-logística(alpha==1~","~gamma==3)))  # Texto da legenda
#ggsave("imagens/ref_h_llogis.pdf", width = 158, height = 93, units = "mm")

