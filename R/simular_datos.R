# simular_datos.R ---------------------------------------------------------
# Genera un dataset FICTICIO con la misma estructura que el del estudio.
# Los datos reales no se publican: son de pacientes.
#
# Qué reproduce (de forma aproximada):
#   - n = 2409 pacientes.
#   - La prevalencia de sensibilización (IgE específica >= 0,35 kUA/L) a cada
#     componente del panel: 13 nsLTP y el extracto de Cannabis sativa (Can s).
#   - La co-sensibilización: Can s 3 casi nunca aparece aislado; va con las
#     LTP de alimentos y pólenes. Par j 2 va bastante por libre.
#   - Los formatos de laboratorio ya convertidos a número:
#       IgE específica "< 0,10"  -> 0,05 (la mitad del límite de detección)
#       IgE total      "< 20"    -> 10;  "> 2500" -> 2500
# Las cifras no serán idénticas a las publicadas.
#
# Uso: source("R/simular_datos.R") crea data/cannabis_ltp_simulado.csv
# -------------------------------------------------------------------------

set.seed(2026)
n <- 2409

# Prevalencia objetivo (%) y dependencia de la sensibilización LTP común
# (alto = muy co-sensibilizado con el resto; bajo = va por libre)
components <- data.frame(
  var      = c("pru_p_3", "mal_d_3", "zea_m_14", "ara_h_9", "vit_v_1", "pla_a_3",
               "art_v_3", "cor_a_8", "par_j_2", "jug_r_3", "can_s_3", "act_d_2",
               "tri_a_14", "can_s"),
  prev     = c(11.7, 10.0, 8.8, 8.1, 6.8, 6.8, 6.7, 6.3, 5.3, 5.1, 3.7, 1.7, 1.2, 2.2) / 100,
  loading  = c(3.2, 3.4, 3.2, 3.2, 3.4, 3.2, 3.0, 3.4, 0.1, 3.2, 3.6, 1.6, 2.2, 3.0)
)

z <- rnorm(n)  # sensibilización LTP latente de cada paciente

ltp <- data.frame(id = sprintf("P%04d", seq_len(n)))
for (i in seq_len(nrow(components))) {
  lam <- components$loading[i]
  # Umbral de la variable latente que da la prevalencia buscada
  b <- qnorm(components$prev[i]) * sqrt(lam^2 + 1)
  positive  <- lam * z + b + rnorm(n) > 0
  level_pos <- 0.35 + exp(rnorm(n, log(0.8) + 0.4 * pmax(z, 0), 1.0))
  level_neg <- ifelse(runif(n) < 0.85, 0.05, runif(n, 0.10, 0.34))
  ltp[[components$var[i]]] <- round(ifelse(positive, level_pos, level_neg), 2)
}

# IgE total: sube con el número de componentes positivos; censurada en 20 y 2500
n_pos <- rowSums(ltp[components$var] >= 0.35)
total_ige <- exp(log(70) + 0.25 * n_pos + rnorm(n, 0, 1.2))
total_ige <- ifelse(total_ige < 20, 10, pmin(total_ige, 2500))
ltp$total_ige <- round(total_ige, 1)

ltp <- ltp[, c("id", "total_ige", components$var)]

dir.create("data", showWarnings = FALSE)
write.csv(ltp, "data/cannabis_ltp_simulado.csv", row.names = FALSE)
