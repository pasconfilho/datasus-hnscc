# ============================================================================
#  ANÁLISE DE SOBREVIDA — HNSCC NO SUS (2010–2023)  •  VERSÃO 4.0 FINAL
#  Fonte: IntegradorRHC (INCA/DATASUS)
#  Nível: The Lancet / Oral Oncology
#  Última atualização: março/2026
#
#  CONTEÚDO COMPLETO:
#  ✓ Importação e limpeza dos .dbf
#  ✓ Construção completa de variáveis (idade, região, tratamento, etc.)
#  ✓ Diagnóstico de viés de referenciamento regional
#  ✓ KM globais: tratamento, estadio, sítio, região, raça, coorte, tempo diag
#  ✓ KM por tabagismo e etilismo (globais)
#  ✓ KM estratificados por sítio anatômico:
#      → por tratamento, região, estadio, tempo diag-trat,
#        tabagismo, etilismo
#      → Exposição × Tratamento dentro de cada sítio (painel 2×3)
#  ✓ Painéis compostos: trat/região/tabaco/álcool por sítio
#  ✓ Disparidade regional: barras, violinos, acesso RT, heatmap
#  ✓ Mapas coropléticos (geobr): casos, sobrevida 3a, RT, tempo diag
#  ✓ Regressão de Cox multivariada + forest plot
#  ✓ Tabela 1 (gtsummary)
#  ✓ Banco analítico salvo em .rds e .csv
#
#  ARQUITETURA:
#  - tema_pub  → ggplot2 puro (usa element_markdown)
#  - tema_surv → ggsurvplot (sem element_markdown — evita conflito)
#  - salvar()    → ggplot2 puro
#  - salvar_km() → objetos ggsurvplot com risk table (png/print/dev.off)
#  - theme_set() NÃO é chamado globalmente
#  - wrap_plots sempre com old_theme guard
# ============================================================================


# ── 0. PACOTES ───────────────────────────────────────────────────────────────
# install.packages(c(
#   "foreign","tidyverse","lubridate","janitor","survival","survminer",
#   "gtsummary","ggplot2","patchwork","viridis","scales","broom",
#   "geobr","sf","ggtext"
# ))

library(foreign);   library(tidyverse);  library(lubridate); library(janitor)
library(survival);  library(survminer);  library(gtsummary)
library(ggplot2);   library(patchwork);  library(viridis)
library(scales);    library(broom);      library(geobr)
library(sf);        library(ggtext)

dir.create("graficos", showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# ── PALETAS ───────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

PALETA_TRAT <- c(
  "Cirurgia (+-adj)" = "#003C67",
  "QT+RT"            = "#E6A817",
  "RT exclusiva"     = "#C0392B",
  "QT exclusiva"     = "#7F8C8D",
  "Nenhum"           = "#BDC3C7"
)

PALETA_REGIAO <- c(
  "Norte"        = "#4E79A7",
  "Nordeste"     = "#F28E2B",
  "Sudeste"      = "#E15759",
  "Sul"          = "#499894",
  "Centro-Oeste" = "#59A14F"
)

PALETA_ESTADIO <- c(
  "0 (in situ)" = "#27AE60",
  "I"           = "#2980B9",
  "II"          = "#F39C12",
  "III"         = "#E74C3C",
  "IV"          = "#6E2222"
)

PALETA_SITIO <- c(
  "Cavidade Oral" = "#003C67",
  "Orofaringe"    = "#E6A817",
  "Laringe"       = "#C0392B",
  "Hipofaringe"   = "#1A7A8A",
  "Nasofaringe"   = "#7F8C8D"
)

PALETA_EXPO <- c(
  "Nunca"        = "#27AE60",
  "Ex-fumante"   = "#F39C12",
  "Fumante"      = "#C0392B",
  "Ex-consumidor"= "#F39C12",
  "Sim"          = "#C0392B"
)


# ══════════════════════════════════════════════════════════════════════════════
# ── TEMAS ─────────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

# Para ggplot2 puro — usa element_markdown
tema_pub <- theme_bw(base_size = 13) +
  theme(
    plot.title        = element_markdown(face = "bold", size = 14, margin = margin(b=4)),
    plot.subtitle     = element_text(color = "grey35", size = 11, margin = margin(b=8)),
    plot.caption      = element_text(color = "grey50", size = 9, hjust = 0, margin = margin(t=8)),
    axis.title        = element_text(face = "bold", size = 12),
    axis.text         = element_text(size = 11),
    axis.line         = element_line(color = "grey30"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 11),
    legend.text       = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "grey85"),
    strip.background  = element_rect(fill = "grey15", color = NA),
    strip.text        = element_text(color = "white", face = "bold", size = 11),
    plot.margin       = margin(12, 16, 12, 12)
  )

# Para ggsurvplot — SEM element_markdown (evita erro de merge de tema)
tema_surv <- theme_bw(base_size = 13) +
  theme(
    plot.title        = element_text(face = "bold", size = 14, margin = margin(b=4)),
    plot.subtitle     = element_text(color = "grey35", size = 11, margin = margin(b=8)),
    plot.caption      = element_text(color = "grey50", size = 9, hjust = 0),
    axis.title        = element_text(face = "bold", size = 12),
    axis.text         = element_text(size = 11),
    axis.line         = element_line(color = "grey30"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 11),
    legend.text       = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "grey85"),
    plot.margin       = margin(12, 16, 12, 12)
  )

# Tema para painéis pequenos dentro de wrap_plots
tema_mini <- tema_surv + theme(
  plot.title      = element_text(size = 11, face = "bold"),
  axis.title      = element_text(size = 9),
  axis.text       = element_text(size = 8),
  legend.text     = element_text(size = 7),
  legend.position = "bottom",
  legend.key.size = unit(0.35, "cm")
)

# Tema para mapas
tema_mapa <- theme_void(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle   = element_text(size = 10, color = "grey40", hjust = 0.5),
    plot.caption    = element_text(size = 9, color = "grey50", hjust = 0),
    legend.position = "right"
  )


# ══════════════════════════════════════════════════════════════════════════════
# ── FUNÇÕES UTILITÁRIAS ───────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

# Salvar ggplot2 puro
salvar <- function(p, nome, w = 14, h = 9, dpi = 300) {
  ggsave(paste0("graficos/", nome, ".png"), p,
         width = w, height = h, dpi = dpi, bg = "white")
  message("  ✓ ", nome)
}

# Salvar ggsurvplot (com ou sem risk table)
salvar_km <- function(obj, nome, w = 14, h = 10, dpi = 300) {
  png(paste0("graficos/", nome, ".png"),
      width = w * dpi, height = h * dpi, res = dpi)
  print(obj)
  dev.off()
  message("  ✓ ", nome)
}

# Gerar objeto ggsurvplot padronizado
km_pub <- function(formula, data, titulo, subtitulo = NULL,
                   paleta, xlim = 5, break_by = 1,
                   risk_table = TRUE, pval_coord = c(0.3, 0.08)) {

  fit      <- do.call(survfit, list(formula = formula, data = data))
  labs_uso <- gsub("^.*=", "", names(fit$strata))

  gg <- ggsurvplot(
    fit, data = data,
    palette           = paleta,
    conf.int          = TRUE,
    conf.int.alpha    = 0.12,
    conf.int.style    = "ribbon",
    risk.table        = risk_table,
    risk.table.col    = "strata",
    risk.table.height = 0.28,
    risk.table.y.text = FALSE,
    pval              = TRUE,
    pval.size         = 4,
    pval.coord        = pval_coord,
    log.rank.weights  = "1",
    xlim              = c(0, xlim),
    break.time.by     = break_by,
    xlab              = "Tempo desde o diagnóstico (anos)",
    ylab              = "Probabilidade de Sobrevida Global",
    surv.median.line  = "hv",
    title             = titulo,
    legend.title      = "",
    legend.labs       = labs_uso,
    ggtheme           = tema_surv,
    fontsize          = 4,
    censor            = TRUE,
    censor.shape      = "|",
    censor.size       = 2
  )

  gg$plot <- gg$plot +
    labs(subtitle = subtitulo) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1), oob = scales::squish) +
    theme(legend.position = "right")

  if (risk_table)
    gg$table <- gg$table +
      theme(axis.title.y = element_blank(),
            plot.title   = element_text(size = 10, face = "bold", color = "grey30"))
  gg
}

# Gerar mini-KM para painéis compostos (sem risk table)
km_mini <- function(formula, data, titulo, paleta, legend_labs) {
  fit <- do.call(survfit, list(formula = formula, data = data))

  ggsurvplot(
    fit, data = data,
    palette          = paleta,
    conf.int         = TRUE, conf.int.alpha = 0.10,
    pval             = TRUE, pval.size = 3,
    risk.table       = FALSE,
    xlim             = c(0, 5), break.time.by = 1,
    surv.median.line = "hv",
    xlab = "Anos", ylab = "S(t)",
    title = titulo, legend.title = "", legend.labs = legend_labs,
    ggtheme = tema_mini
  )$plot +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1), oob = scales::squish)
}

# Salvar wrap_plots com guard de tema
salvar_painel <- function(lista_plots, nome, titulo, subtitulo, caption,
                          ncol = 3, w = 18, h = 14) {
  lista_plots <- compact(lista_plots)
  if (length(lista_plots) < 2) {
    message("  ⚠ Painel '", nome, "' ignorado — plots insuficientes")
    return(invisible(NULL))
  }
  old <- theme_get(); theme_set(theme_bw(base_size = 12))
  p <- wrap_plots(lista_plots, ncol = ncol) +
    plot_annotation(
      title = titulo, subtitle = subtitulo, caption = caption,
      theme = theme(
        plot.title    = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        plot.caption  = element_text(size = 9,  color = "grey50")
      )
    )
  salvar(p, nome, w = w, h = h)
  theme_set(old)
}


# ══════════════════════════════════════════════════════════════════════════════
# ── 1. IMPORTAR E COMBINAR .DBF ───────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

PASTA_DADOS <- "~/Downloads/download_t...3105844669/"   # ← AJUSTE AQUI

arquivos_dbf <- list.files(PASTA_DADOS, pattern = "rhc\\d+\\.dbf$",
                            full.names = TRUE)
cat("Arquivos encontrados:", length(arquivos_dbf), "\n")

rhc_raw <- arquivos_dbf |>
  map(\(f) { cat("  Lendo:", basename(f), "\n"); read.dbf(f, as.is = TRUE) }) |>
  list_rbind()

cat("Total bruto:", nrow(rhc_raw), "\n")


# ══════════════════════════════════════════════════════════════════════════════
# ── 2. FILTRAR HNSCC EPIDERMÓIDE ─────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

hnscc_cids <- c(paste0("C", sprintf("%02d", 0:14)), "C30","C31","C32")
hnscc_hist <- c("8070","8071","8072","8073","8074","8075","8076","8077")

rhc_hnscc <- rhc_raw |>
  clean_names() |>
  mutate(across(where(is.character), str_trim)) |>
  filter(
    tpcaso == 1,
    loctudet %in% hnscc_cids,
    str_sub(tipohist, 1, 4) %in% hnscc_hist
  )

cat("HNSCC epidermóide analítico:", nrow(rhc_hnscc), "\n")


# ══════════════════════════════════════════════════════════════════════════════
# ── 3. CONSTRUÇÃO DE VARIÁVEIS ────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

DATA_CORTE <- as.Date("2023-12-31")

rhc_surv <- rhc_hnscc |>
  mutate(
    # Datas
    dt_consulta = dmy(datapricon),
    dt_obito    = dmy(dataobito),
    dt_trat     = dmy(datainitrt),
    dt_diag     = dmy(dtdiagno),
    dt_inicio   = coalesce(dt_consulta, dt_diag),
    ano_diag    = year(dt_inicio),

    # Desfecho
    dt_saida   = coalesce(dt_obito, DATA_CORTE),
    status     = if_else(!is.na(dt_obito), 1L, 0L),
    tempo_dias = as.numeric(dt_saida - dt_inicio),
    tempo_anos = tempo_dias / 365.25,

    # Tempo diagnóstico → tratamento
    tempo_diag_trat = as.numeric(dt_trat - dt_inicio),
    tempo_diag_trat_cat = case_when(
      tempo_diag_trat <= 30 ~ "<=30 dias",
      tempo_diag_trat <= 60 ~ "31-60 dias",
      tempo_diag_trat <= 90 ~ "61-90 dias",
      tempo_diag_trat >  90 ~ ">90 dias",
      TRUE                  ~ NA_character_
    ) |> factor(levels = c("<=30 dias","31-60 dias","61-90 dias",">90 dias")),

    # Sexo
    sexo_label = case_when(
      sexo == 1 ~ "Masculino", sexo == 2 ~ "Feminino", TRUE ~ NA_character_
    ) |> factor(levels = c("Masculino","Feminino")),

    # Idade — corrige "050" → 50
    idade_num = as.numeric(as.character(idade)),
    faixa_etaria = case_when(
      idade_num < 40  ~ "<40 anos",
      idade_num < 50  ~ "40-49 anos",
      idade_num < 60  ~ "50-59 anos",
      idade_num < 70  ~ "60-69 anos",
      idade_num >= 70 ~ ">=70 anos",
      TRUE            ~ NA_character_
    ) |> factor(levels = c("<40 anos","40-49 anos","50-59 anos","60-69 anos",">=70 anos")),

    # Raça/Cor
    raca_grupo = case_when(
      racacor == 1          ~ "Branca",
      racacor %in% c(2, 4) ~ "Preta/Parda",
      racacor %in% c(3, 5) ~ "Outras",
      TRUE                  ~ NA_character_
    ) |> factor(levels = c("Branca","Preta/Parda","Outras")),

    # Estadiamento
    estadio_label = case_when(
      estadiam %in% c("1","1A","1B","1C","01","I")                 ~ "I",
      estadiam %in% c("2","2A","2B","2C","2D","02","II")           ~ "II",
      estadiam %in% c("3","3A","3B","3C","03","III")               ~ "III",
      estadiam %in% c("4","4A","4B","4C","4D","04","IV","44","49") ~ "IV",
      estadiam == "0"                                               ~ "0 (in situ)",
      TRUE                                                          ~ "Desconhecido"
    ) |> factor(levels = c("0 (in situ)","I","II","III","IV","Desconhecido")),

    estadio_grupo = case_when(
      estadio_label %in% c("0 (in situ)","I","II") ~ "Precoce (0-II)",
      estadio_label %in% c("III","IV")             ~ "Avancado (III-IV)",
      TRUE                                         ~ "Desconhecido"
    ) |> factor(levels = c("Precoce (0-II)","Avancado (III-IV)","Desconhecido")),

    # Tratamento
    trat_raw      = as.character(pritrath),
    tem_cirurgia  = str_detect(trat_raw, "2"),
    tem_rt        = str_detect(trat_raw, "3"),
    tem_qt        = str_detect(trat_raw, "4"),
    nenhum_trat   = trat_raw == "1",
    sem_info_trat = trat_raw == "9" | is.na(pritrath),

    trat_grupo = case_when(
      sem_info_trat                                    ~ "Sem informacao",
      nenhum_trat                                      ~ "Nenhum",
      tem_cirurgia & (tem_rt | tem_qt)                 ~ "Cirurgia (+-adj)",
      tem_cirurgia & !tem_rt & !tem_qt                 ~ "Cirurgia (+-adj)",
      !tem_cirurgia & tem_rt & tem_qt                  ~ "QT+RT",
      !tem_cirurgia & tem_rt & !tem_qt                 ~ "RT exclusiva",
      !tem_cirurgia & !tem_rt & tem_qt                 ~ "QT exclusiva",
      TRUE                                             ~ "Sem informacao"
    ) |> factor(levels = c("Cirurgia (+-adj)","QT+RT","RT exclusiva",
                            "QT exclusiva","Nenhum","Sem informacao")),

    # Sítio anatômico
    sitio = case_when(
      loctudet %in% paste0("C", sprintf("%02d", 0:6)) ~ "Cavidade Oral",
      loctudet %in% c("C07","C08")                    ~ "Gl. Salivares",
      loctudet %in% c("C09","C10")                    ~ "Orofaringe",
      loctudet == "C11"                               ~ "Nasofaringe",
      loctudet %in% c("C12","C13")                    ~ "Hipofaringe",
      loctudet == "C14"                               ~ "Faringe NE",
      loctudet %in% c("C30","C31")                    ~ "Cav. Nasal/Seios",
      loctudet == "C32"                               ~ "Laringe",
      TRUE                                            ~ "Outros"
    ) |> factor(levels = c("Cavidade Oral","Orofaringe","Laringe","Hipofaringe",
                            "Nasofaringe","Gl. Salivares","Cav. Nasal/Seios",
                            "Faringe NE","Outros")),

    # Região de residência e do hospital
    uf_res_clean  = str_trim(estadres),
    uf_hosp_clean = str_trim(ufuh),

    regiao = case_when(
      uf_res_clean %in% c("AC","AM","RO","RR","PA","AP","TO")           ~ "Norte",
      uf_res_clean %in% c("MA","PI","CE","RN","PB","PE","AL","SE","BA") ~ "Nordeste",
      uf_res_clean %in% c("MG","ES","RJ","SP")                          ~ "Sudeste",
      uf_res_clean %in% c("PR","SC","RS")                               ~ "Sul",
      uf_res_clean %in% c("MT","MS","GO","DF")                          ~ "Centro-Oeste",
      TRUE                                                               ~ NA_character_
    ) |> factor(levels = c("Norte","Nordeste","Sudeste","Sul","Centro-Oeste")),

    regiao_hosp = case_when(
      uf_hosp_clean %in% c("AC","AM","RO","RR","PA","AP","TO")           ~ "Norte",
      uf_hosp_clean %in% c("MA","PI","CE","RN","PB","PE","AL","SE","BA") ~ "Nordeste",
      uf_hosp_clean %in% c("MG","ES","RJ","SP")                          ~ "Sudeste",
      uf_hosp_clean %in% c("PR","SC","RS")                               ~ "Sul",
      uf_hosp_clean %in% c("MT","MS","GO","DF")                          ~ "Centro-Oeste",
      TRUE                                                                ~ NA_character_
    ) |> factor(levels = c("Norte","Nordeste","Sudeste","Sul","Centro-Oeste")),

    tratado_fora = !is.na(regiao) & !is.na(regiao_hosp) & (regiao != regiao_hosp),

    # Tabagismo
    tabaco_label = case_when(
      tabagism == 1 ~ "Nunca",
      tabagism == 2 ~ "Ex-fumante",
      tabagism == 3 ~ "Fumante",
      TRUE          ~ NA_character_
    ) |> factor(levels = c("Nunca","Ex-fumante","Fumante")),

    # Etilismo
    alcool_label = case_when(
      alcoolis == 1 ~ "Nunca",
      alcoolis == 2 ~ "Ex-consumidor",
      alcoolis == 3 ~ "Sim",
      TRUE          ~ NA_character_
    ) |> factor(levels = c("Nunca","Ex-consumidor","Sim")),

    # Coorte temporal
    coorte = case_when(
      ano_diag %in% 2010:2012 ~ "2010-2012",
      ano_diag %in% 2013:2015 ~ "2013-2015",
      ano_diag %in% 2016:2019 ~ "2016-2019",
      ano_diag %in% 2020:2023 ~ "2020-2023",
      TRUE                    ~ NA_character_
    ) |> factor(levels = c("2010-2012","2013-2015","2016-2019","2020-2023"))

  ) |>
  filter(!is.na(dt_inicio), tempo_dias > 0, ano_diag >= 2010)

cat("\nRegistros analiticos:", nrow(rhc_surv), "\n")
cat("Obitos:", sum(rhc_surv$status),
    "(", round(mean(rhc_surv$status)*100, 1), "%)\n")


# ══════════════════════════════════════════════════════════════════════════════
# ── 4. DIAGNÓSTICO REGIONAL ───────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── % tratados fora da regiao de domicilio ──\n")
rhc_surv |>
  filter(!is.na(regiao), !is.na(regiao_hosp)) |>
  group_by(regiao) |>
  summarise(n = n(), fora = sum(tratado_fora),
            pct_fora = round(mean(tratado_fora)*100, 1)) |>
  print()

cat("\n── % Estadio III/IV por regiao ──\n")
rhc_surv |>
  filter(!is.na(regiao), estadio_label != "Desconhecido") |>
  count(regiao, estadio_label) |>
  group_by(regiao) |>
  mutate(pct = n/sum(n)*100) |>
  filter(estadio_label %in% c("III","IV")) |>
  summarise(pct_avancado = round(sum(pct), 1)) |>
  arrange(desc(pct_avancado)) |>
  print()

p_fluxo <- rhc_surv |>
  filter(!is.na(regiao), !is.na(regiao_hosp), tratado_fora) |>
  count(regiao_origem = regiao, regiao_destino = regiao_hosp) |>
  group_by(regiao_origem) |>
  mutate(pct = n/sum(n)) |>
  ggplot(aes(x = regiao_destino, y = regiao_origem,
             size = n, color = regiao_origem)) +
  geom_point(alpha = 0.85) +
  geom_text(aes(label = comma(n)), size = 3, color = "white", fontface = "bold") +
  scale_size_continuous(name = "N referenciados", range = c(6,28), labels = comma) +
  scale_color_manual(values = PALETA_REGIAO, guide = "none") +
  scale_x_discrete(position = "top") +
  labs(
    title    = "Referenciamento Inter-Regional — HNSCC/SUS",
    subtitle = "Pacientes tratados fora de sua regiao de domicilio",
    x = "Regiao do hospital", y = "Regiao de residencia",
    caption = "Fonte: IntegradorRHC/INCA • 2010-2023"
  ) + tema_pub

salvar(p_fluxo, "00_fluxo_referenciamento", w = 12, h = 7)


# ══════════════════════════════════════════════════════════════════════════════
# ── 5. TABELA 1 ───────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

tabela1 <- rhc_surv |>
  select(sexo_label, faixa_etaria, raca_grupo, regiao,
         tabaco_label, alcool_label, sitio, estadio_label,
         trat_grupo, tempo_diag_trat_cat, coorte, status, tempo_anos) |>
  tbl_summary(
    missing = "no",
    label = list(
      sexo_label          ~ "Sexo",
      faixa_etaria        ~ "Faixa etaria",
      raca_grupo          ~ "Raca/Cor",
      regiao              ~ "Regiao de residencia",
      tabaco_label        ~ "Tabagismo",
      alcool_label        ~ "Etilismo",
      sitio               ~ "Sitio anatomico",
      estadio_label       ~ "Estadiamento clinico",
      trat_grupo          ~ "Modalidade de tratamento",
      tempo_diag_trat_cat ~ "Tempo diag-tratamento",
      coorte              ~ "Coorte de diagnostico",
      status              ~ "Obito registrado",
      tempo_anos          ~ "Seguimento (anos)"
    ),
    statistic = list(
      all_continuous()  ~ "{median} ({p25}-{p75})",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) |>
  bold_labels() |>
  modify_caption("**Tabela 1.** Caracteristicas da coorte HNSCC/SUS 2010-2023")

print(tabela1)


# ══════════════════════════════════════════════════════════════════════════════
# ── 6. SOBREVIDA GLOBAL — CURVAS PRINCIPAIS ───────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Sobrevida global em 1, 3 e 5 anos ──\n")
km_geral <- survfit(Surv(tempo_anos, status) ~ 1, data = rhc_surv)
s <- summary(km_geral, times = c(1,3,5))
print(data.frame(anos = s$time, surv = round(s$surv,3),
                 ic95_low = round(s$lower,3), ic95_high = round(s$upper,3)))

# ── 6a. Por tratamento
dados_trat <- rhc_surv |>
  filter(trat_grupo %in% names(PALETA_TRAT)) |> droplevels()
km_01 <- km_pub(Surv(tempo_anos, status) ~ trat_grupo, dados_trat,
                titulo = "Sobrevida Global por Modalidade de Tratamento",
                subtitulo = "HNSCC/SUS 2010-2023 • Todos os sitios",
                paleta = unname(PALETA_TRAT[levels(dados_trat$trat_grupo)]))
salvar_km(km_01, "01_km_tratamento")

# ── 6b. Por estadiamento
dados_est <- rhc_surv |> filter(estadio_label != "Desconhecido") |> droplevels()
km_02 <- km_pub(Surv(tempo_anos, status) ~ estadio_label, dados_est,
                titulo = "Sobrevida Global por Estadiamento Clinico",
                subtitulo = "HNSCC/SUS 2010-2023",
                paleta = unname(PALETA_ESTADIO[levels(dados_est$estadio_label)]))
salvar_km(km_02, "02_km_estadiamento")

# ── 6c. Por sítio (todos juntos)
dados_sitio_p <- rhc_surv |>
  filter(sitio %in% names(PALETA_SITIO)) |> droplevels()
km_03 <- km_pub(Surv(tempo_anos, status) ~ sitio, dados_sitio_p,
                titulo = "Sobrevida Global por Sitio Anatomico",
                subtitulo = "HNSCC/SUS 2010-2023 • Principais localizacoes",
                paleta = unname(PALETA_SITIO[levels(dados_sitio_p$sitio)]))
salvar_km(km_03, "03_km_sitio_todos")

# ── 6d. Por sítio — painéis 1 vs. resto
sitios_analise <- names(PALETA_SITIO)
nomes_arq      <- c("cavidade_oral","orofaringe","laringe","hipofaringe","nasofaringe")

plots_04 <- map(sitios_analise, function(s) {
  d <- rhc_surv |>
    filter(!is.na(sitio)) |>
    mutate(grupo = factor(if_else(sitio == s, s, "Outros sitios"),
                          levels = c(s, "Outros sitios")))
  km_mini(Surv(tempo_anos, status) ~ grupo, d, s,
          paleta = c("#003C67","#C0392B"),
          legend_labs = c(s, "Outros sitios"))
})
salvar_painel(plots_04, "04_km_sitio_paineis",
              "Sobrevida por Sitio Anatomico — HNSCC/SUS 2010-2023",
              "Cada painel: sitio indicado vs. todos os demais",
              "Fonte: IntegradorRHC/INCA • p-valor: teste log-rank")

# ── 6e. Por região
km_05 <- km_pub(Surv(tempo_anos, status) ~ regiao,
                rhc_surv |> filter(!is.na(regiao)) |> droplevels(),
                titulo = "Sobrevida Global por Regiao de Residencia",
                subtitulo = "HNSCC/SUS 2010-2023 • Por UF de domicilio",
                paleta = unname(PALETA_REGIAO))
salvar_km(km_05, "05_km_regiao")

# ── 6f. Tratamento em estádio avançado
dados_av <- rhc_surv |>
  filter(estadio_grupo == "Avancado (III-IV)",
         trat_grupo %in% c("QT+RT","RT exclusiva","Cirurgia (+-adj)")) |>
  droplevels()
km_06 <- km_pub(Surv(tempo_anos, status) ~ trat_grupo, dados_av,
                titulo = "Sobrevida por Tratamento — Estadio III/IV",
                subtitulo = "HNSCC/SUS 2010-2023 • Casos avancados",
                paleta = unname(PALETA_TRAT[c("Cirurgia (+-adj)","QT+RT","RT exclusiva")]))
salvar_km(km_06, "06_km_avancado_trat")

# ── 6g. Por raça/cor
km_07 <- km_pub(Surv(tempo_anos, status) ~ raca_grupo,
                rhc_surv |> filter(!is.na(raca_grupo)) |> droplevels(),
                titulo = "Sobrevida Global por Raca/Cor",
                subtitulo = "HNSCC/SUS 2010-2023",
                paleta = c("#003C67","#E6A817","#7F8C8D"))
salvar_km(km_07, "07_km_raca")

# ── 6h. Por coorte (xlim=3 — evita artefato de seguimento curto)
km_08 <- km_pub(Surv(tempo_anos, status) ~ coorte,
                rhc_surv |> filter(!is.na(coorte)) |> droplevels(),
                titulo = "Tendencia de Sobrevida por Coorte de Diagnostico",
                subtitulo = "HNSCC/SUS 2010-2023 • Seguimento limitado a 3 anos",
                paleta = c("#001E3C","#003C67","#0073C2","#7AA6DC"),
                xlim = 3, break_by = 0.5)
salvar_km(km_08, "08_km_coorte")

# ── 6i. Por tempo diagnóstico-tratamento
km_09a <- km_pub(Surv(tempo_anos, status) ~ tempo_diag_trat_cat,
                 rhc_surv |> filter(!is.na(tempo_diag_trat_cat)) |> droplevels(),
                 titulo = "Sobrevida por Tempo Diagnostico-Tratamento",
                 subtitulo = "HNSCC/SUS 2010-2023",
                 paleta = c("#27AE60","#F39C12","#E74C3C","#6E2222"))
salvar_km(km_09a, "09a_km_tempo_diagnostico")

# ── 6j. Por tabagismo (global)
km_09b <- km_pub(Surv(tempo_anos, status) ~ tabaco_label,
                 rhc_surv |> filter(!is.na(tabaco_label)) |> droplevels(),
                 titulo = "Sobrevida Global por Historico de Tabagismo",
                 subtitulo = "HNSCC/SUS 2010-2023",
                 paleta = c("#27AE60","#F39C12","#C0392B"))
salvar_km(km_09b, "09b_km_tabagismo")

# ── 6k. Por etilismo (global)
km_09c <- km_pub(Surv(tempo_anos, status) ~ alcool_label,
                 rhc_surv |> filter(!is.na(alcool_label)) |> droplevels(),
                 titulo = "Sobrevida Global por Historico de Etilismo",
                 subtitulo = "HNSCC/SUS 2010-2023",
                 paleta = c("#27AE60","#F39C12","#C0392B"))
salvar_km(km_09c, "09c_km_etilismo")

# ── 6l. Painel exposição × tratamento (global — tabaco linha 1, álcool linha 2)
niveis_tab <- c("Nunca","Ex-fumante","Fumante")
niveis_alc <- c("Nunca","Ex-consumidor","Sim")

plots_09d <- c(
  map(niveis_tab, function(nivel) {
    d <- rhc_surv |>
      filter(tabaco_label == nivel, trat_grupo %in% names(PALETA_TRAT)) |>
      droplevels()
    if (nrow(d) < 100 || n_distinct(d$trat_grupo) < 2) return(NULL)
    km_mini(Surv(tempo_anos, status) ~ trat_grupo, d,
            paste0("Tabaco: ", nivel),
            unname(PALETA_TRAT[levels(d$trat_grupo)]),
            gsub("^.*=", "", names(survfit(Surv(tempo_anos, status) ~ trat_grupo,
                                           data = d)$strata)))
  }),
  map(niveis_alc, function(nivel) {
    d <- rhc_surv |>
      filter(alcool_label == nivel, trat_grupo %in% names(PALETA_TRAT)) |>
      droplevels()
    if (nrow(d) < 100 || n_distinct(d$trat_grupo) < 2) return(NULL)
    km_mini(Surv(tempo_anos, status) ~ trat_grupo, d,
            paste0("Alcool: ", nivel),
            unname(PALETA_TRAT[levels(d$trat_grupo)]),
            gsub("^.*=", "", names(survfit(Surv(tempo_anos, status) ~ trat_grupo,
                                           data = d)$strata)))
  })
)
salvar_painel(plots_09d, "09d_painel_exposicao_tratamento",
              "Sobrevida por Tratamento Estratificada por Exposicao",
              "Linha 1: Tabagismo | Linha 2: Etilismo — HNSCC/SUS 2010-2023",
              "Cada painel = nivel de exposicao; curvas = modalidade de tratamento\nFonte: IntegradorRHC/INCA",
              ncol = 3, w = 18, h = 14)


# ══════════════════════════════════════════════════════════════════════════════
# ── 7. ANÁLISES ESTRATIFICADAS POR SÍTIO ANATÔMICO ───────────────────────────
# Para cada sítio: KM por tratamento, região, estadio, tempo, tabaco, álcool,
# e painel composto Exposição × Tratamento
# ══════════════════════════════════════════════════════════════════════════════

for (i in seq_along(sitios_analise)) {

  s   <- sitios_analise[i]
  arq <- nomes_arq[i]
  d   <- rhc_surv |> filter(sitio == s)

  cat("\n══ Sitio:", s, "| n =", nrow(d),
      "| obitos =", sum(d$status), "══\n")

  # ── 10a. KM por tratamento
  d_trat <- d |> filter(trat_grupo %in% names(PALETA_TRAT)) |> droplevels()
  if (n_distinct(d_trat$trat_grupo) >= 2) {
    km_t <- km_pub(Surv(tempo_anos, status) ~ trat_grupo, d_trat,
                   titulo = paste0("Sobrevida por Tratamento — ", s),
                   subtitulo = "HNSCC/SUS 2010-2023",
                   paleta = unname(PALETA_TRAT[levels(d_trat$trat_grupo)]))
    salvar_km(km_t, paste0("10a_km_trat_", arq))
  }

  # ── 10b. KM por região
  d_reg <- d |> filter(!is.na(regiao)) |> droplevels()
  if (n_distinct(d_reg$regiao) >= 2) {
    km_r <- km_pub(Surv(tempo_anos, status) ~ regiao, d_reg,
                   titulo = paste0("Sobrevida por Regiao — ", s),
                   subtitulo = "HNSCC/SUS 2010-2023",
                   paleta = unname(PALETA_REGIAO[levels(d_reg$regiao)]))
    salvar_km(km_r, paste0("10b_km_regiao_", arq))
  }

  # ── 10c. KM por estadiamento
  d_est <- d |> filter(estadio_label != "Desconhecido") |> droplevels()
  if (n_distinct(d_est$estadio_label) >= 2) {
    km_e <- km_pub(Surv(tempo_anos, status) ~ estadio_label, d_est,
                   titulo = paste0("Sobrevida por Estadiamento — ", s),
                   subtitulo = "HNSCC/SUS 2010-2023",
                   paleta = unname(PALETA_ESTADIO[levels(d_est$estadio_label)]))
    salvar_km(km_e, paste0("10c_km_estadio_", arq))
  }

  # ── 10d. KM por tempo diagnóstico-tratamento
  d_td <- d |> filter(!is.na(tempo_diag_trat_cat)) |> droplevels()
  if (n_distinct(d_td$tempo_diag_trat_cat) >= 2) {
    km_td <- km_pub(Surv(tempo_anos, status) ~ tempo_diag_trat_cat, d_td,
                    titulo = paste0("Sobrevida por Tempo Diag-Trat — ", s),
                    subtitulo = "HNSCC/SUS 2010-2023",
                    paleta = c("#27AE60","#F39C12","#E74C3C","#6E2222"))
    salvar_km(km_td, paste0("10d_km_tempo_", arq))
  }

  # ── 10e. KM por tabagismo
  d_tab <- d |> filter(!is.na(tabaco_label)) |> droplevels()
  if (n_distinct(d_tab$tabaco_label) >= 2) {
    km_tab <- km_pub(Surv(tempo_anos, status) ~ tabaco_label, d_tab,
                     titulo = paste0("Sobrevida por Tabagismo — ", s),
                     subtitulo = "HNSCC/SUS 2010-2023",
                     paleta = c("#27AE60","#F39C12","#C0392B"))
    salvar_km(km_tab, paste0("10e_km_tabaco_", arq))
  }

  # ── 10f. KM por etilismo
  d_alc <- d |> filter(!is.na(alcool_label)) |> droplevels()
  if (n_distinct(d_alc$alcool_label) >= 2) {
    km_alc <- km_pub(Surv(tempo_anos, status) ~ alcool_label, d_alc,
                     titulo = paste0("Sobrevida por Etilismo — ", s),
                     subtitulo = "HNSCC/SUS 2010-2023",
                     paleta = c("#27AE60","#F39C12","#C0392B"))
    salvar_km(km_alc, paste0("10f_km_alcool_", arq))
  }

  # ── 10g. Painel Exposição × Tratamento (2×3 dentro de cada sítio)
  plots_expo_trat <- c(
    map(niveis_tab, function(nivel) {
      d_sub <- d |>
        filter(tabaco_label == nivel, trat_grupo %in% names(PALETA_TRAT)) |>
        droplevels()
      if (nrow(d_sub) < 50 || n_distinct(d_sub$trat_grupo) < 2) return(NULL)
      labs <- gsub("^.*=", "", names(
        survfit(Surv(tempo_anos, status) ~ trat_grupo, data = d_sub)$strata))
      km_mini(Surv(tempo_anos, status) ~ trat_grupo, d_sub,
              paste0("Tabaco: ", nivel),
              unname(PALETA_TRAT[levels(d_sub$trat_grupo)]),
              labs)
    }),
    map(niveis_alc, function(nivel) {
      d_sub <- d |>
        filter(alcool_label == nivel, trat_grupo %in% names(PALETA_TRAT)) |>
        droplevels()
      if (nrow(d_sub) < 50 || n_distinct(d_sub$trat_grupo) < 2) return(NULL)
      labs <- gsub("^.*=", "", names(
        survfit(Surv(tempo_anos, status) ~ trat_grupo, data = d_sub)$strata))
      km_mini(Surv(tempo_anos, status) ~ trat_grupo, d_sub,
              paste0("Alcool: ", nivel),
              unname(PALETA_TRAT[levels(d_sub$trat_grupo)]),
              labs)
    })
  )
  salvar_painel(
    plots_expo_trat,
    paste0("10g_km_exposicao_tto_", arq),
    paste0("Sobrevida por Tratamento × Exposicao — ", s),
    "Linha 1: Tabagismo | Linha 2: Etilismo — HNSCC/SUS 2010-2023",
    "Cada painel = nivel de exposicao; curvas = modalidade de tratamento\nFonte: IntegradorRHC/INCA",
    ncol = 3, w = 18, h = 13
  )
}


# ── Painéis compostos (todos os sítios lado a lado) ───────────────────────────

# Tratamento por sítio
salvar_painel(
  map(sitios_analise, function(s) {
    d <- rhc_surv |> filter(sitio == s, trat_grupo %in% names(PALETA_TRAT)) |> droplevels()
    if (n_distinct(d$trat_grupo) < 2) return(NULL)
    labs <- gsub("^.*=", "", names(survfit(Surv(tempo_anos, status) ~ trat_grupo, d)$strata))
    km_mini(Surv(tempo_anos, status) ~ trat_grupo, d, s,
            unname(PALETA_TRAT[levels(d$trat_grupo)]), labs)
  }),
  "11_painel_tratamento_por_sitio",
  "Sobrevida por Tratamento — Estratificada por Sitio Anatomico",
  "HNSCC/SUS 2010-2023", "Fonte: IntegradorRHC/INCA • p-valor: teste log-rank",
  ncol = 3, w = 18, h = 14
)

# Região por sítio
salvar_painel(
  map(sitios_analise, function(s) {
    d <- rhc_surv |> filter(sitio == s, !is.na(regiao)) |> droplevels()
    if (n_distinct(d$regiao) < 2) return(NULL)
    labs <- gsub("^.*=", "", names(survfit(Surv(tempo_anos, status) ~ regiao, d)$strata))
    km_mini(Surv(tempo_anos, status) ~ regiao, d, s,
            unname(PALETA_REGIAO[levels(d$regiao)]), labs)
  }),
  "12_painel_regiao_por_sitio",
  "Sobrevida por Regiao — Estratificada por Sitio Anatomico",
  "HNSCC/SUS 2010-2023", "Fonte: IntegradorRHC/INCA • p-valor: teste log-rank",
  ncol = 3, w = 18, h = 14
)

# Tabagismo por sítio
salvar_painel(
  map(sitios_analise, function(s) {
    d <- rhc_surv |> filter(sitio == s, !is.na(tabaco_label)) |> droplevels()
    if (n_distinct(d$tabaco_label) < 2) return(NULL)
    labs <- gsub("^.*=", "", names(survfit(Surv(tempo_anos, status) ~ tabaco_label, d)$strata))
    km_mini(Surv(tempo_anos, status) ~ tabaco_label, d, s,
            c("#27AE60","#F39C12","#C0392B")[seq_along(labs)], labs)
  }),
  "12b_painel_tabaco_por_sitio",
  "Sobrevida por Tabagismo — Estratificada por Sitio Anatomico",
  "HNSCC/SUS 2010-2023 • Nunca / Ex-fumante / Fumante",
  "Fonte: IntegradorRHC/INCA • p-valor: teste log-rank",
  ncol = 3, w = 18, h = 14
)

# Etilismo por sítio
salvar_painel(
  map(sitios_analise, function(s) {
    d <- rhc_surv |> filter(sitio == s, !is.na(alcool_label)) |> droplevels()
    if (n_distinct(d$alcool_label) < 2) return(NULL)
    labs <- gsub("^.*=", "", names(survfit(Surv(tempo_anos, status) ~ alcool_label, d)$strata))
    km_mini(Surv(tempo_anos, status) ~ alcool_label, d, s,
            c("#27AE60","#F39C12","#C0392B")[seq_along(labs)], labs)
  }),
  "12c_painel_alcool_por_sitio",
  "Sobrevida por Etilismo — Estratificada por Sitio Anatomico",
  "HNSCC/SUS 2010-2023 • Nunca / Ex-consumidor / Sim",
  "Fonte: IntegradorRHC/INCA • p-valor: teste log-rank",
  ncol = 3, w = 18, h = 14
)


# ══════════════════════════════════════════════════════════════════════════════
# ── 8. DISPARIDADE DE TRATAMENTO ─────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

# Tratamento por região
p_13 <- rhc_surv |>
  filter(!is.na(regiao), trat_grupo %in% names(PALETA_TRAT)) |>
  count(regiao, trat_grupo) |>
  group_by(regiao) |> mutate(pct = n/sum(n)) |>
  ggplot(aes(x = regiao, y = pct, fill = trat_grupo)) +
  geom_col(position = "fill", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = if_else(pct >= 0.05, paste0(round(pct*100), "%"), "")),
            position = position_fill(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = PALETA_TRAT, name = "Modalidade") +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
  labs(title = "Distribuicao de Tratamento por Regiao",
       subtitle = "HNSCC/SUS 2010-2023",
       x = NULL, y = "Proporcao", caption = "Fonte: IntegradorRHC/INCA") +
  tema_pub
salvar(p_13, "13_trat_por_regiao", w = 13, h = 7)

# Tratamento por sítio
p_14 <- rhc_surv |>
  filter(sitio %in% sitios_analise, trat_grupo %in% names(PALETA_TRAT)) |>
  count(sitio, trat_grupo) |>
  group_by(sitio) |> mutate(pct = n/sum(n)) |>
  ggplot(aes(x = sitio, y = pct, fill = trat_grupo)) +
  geom_col(position = "fill", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = if_else(pct >= 0.05, paste0(round(pct*100), "%"), "")),
            position = position_fill(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = PALETA_TRAT, name = "Modalidade") +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
  scale_x_discrete(labels = function(x) str_wrap(x, 10)) +
  labs(title = "Distribuicao de Tratamento por Sitio Anatomico",
       subtitle = "HNSCC/SUS 2010-2023",
       x = NULL, y = "Proporcao", caption = "Fonte: IntegradorRHC/INCA") +
  tema_pub
salvar(p_14, "14_trat_por_sitio", w = 14, h = 7)

# Heatmap tratamento × sítio
p_15 <- rhc_surv |>
  filter(sitio %in% sitios_analise, trat_grupo %in% names(PALETA_TRAT)) |>
  count(sitio, trat_grupo) |>
  group_by(sitio) |> mutate(pct = n/sum(n)) |>
  ggplot(aes(x = trat_grupo, y = sitio, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(round(pct*100), "%")),
            fontface = "bold", size = 4.5, color = "white") +
  scale_fill_viridis_c(name = "Proporcao", labels = percent_format(),
                       option = "mako", direction = -1) +
  scale_x_discrete(labels = function(x) str_wrap(x, 10)) +
  labs(title = "Modalidade de Tratamento por Sitio Anatomico",
       subtitle = "Proporcao de pacientes em cada combinacao sitio x tratamento",
       x = "Modalidade", y = NULL, caption = "Fonte: IntegradorRHC/INCA") +
  tema_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1))
salvar(p_15, "15_heatmap_trat_sitio", w = 13, h = 7)

# Acesso a RT por região
p_16 <- rhc_surv |>
  filter(!is.na(regiao)) |>
  group_by(regiao) |>
  summarise(n = n(), pct_rt = mean(tem_rt, na.rm = TRUE),
            ic_low  = pct_rt - 1.96*sqrt(pct_rt*(1-pct_rt)/n),
            ic_high = pct_rt + 1.96*sqrt(pct_rt*(1-pct_rt)/n)) |>
  ggplot(aes(x = fct_reorder(regiao, pct_rt), y = pct_rt, fill = regiao)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ic_low, ymax = ic_high), width = 0.25, color = "grey30") +
  geom_text(aes(label = paste0(round(pct_rt*100,1), "%")),
            hjust = -0.25, fontface = "bold", size = 4.5) +
  scale_fill_manual(values = PALETA_REGIAO) +
  scale_y_continuous(labels = percent_format(), limits = c(0,1),
                     expand = expansion(mult = c(0, 0.12))) +
  coord_flip() +
  labs(title = "Proporcao com Radioterapia por Regiao",
       subtitle = "HNSCC/SUS 2010-2023 • IC 95%",
       x = NULL, y = "% com RT", caption = "Fonte: IntegradorRHC/INCA") +
  tema_pub
salvar(p_16, "16_acesso_rt_regiao", w = 12, h = 6)

# Violino: tempo diag-trat por região
p_17 <- rhc_surv |>
  filter(!is.na(regiao), !is.na(tempo_diag_trat),
         tempo_diag_trat >= 0, tempo_diag_trat <= 365) |>
  ggplot(aes(x = regiao, y = tempo_diag_trat, fill = regiao)) +
  geom_violin(alpha = 0.7, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.14, fill = "white", color = "grey30",
               outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 21,
               size = 3, fill = "white", color = "black") +
  geom_hline(yintercept = 60, linetype = "dashed",
             color = "firebrick", linewidth = 0.8) +
  annotate("text", x = 5.5, y = 63, label = "60 dias (meta)",
           hjust = 0, color = "firebrick", size = 3.5, fontface = "italic") +
  scale_fill_manual(values = PALETA_REGIAO, guide = "none") +
  scale_y_continuous(limits = c(0,365), breaks = seq(0,365,60)) +
  labs(title = "Tempo Diagnostico-Tratamento por Regiao",
       subtitle = "HNSCC/SUS 2010-2023",
       x = NULL, y = "Dias ate inicio do tratamento",
       caption = "Excluidos valores >365 dias. Fonte: IntegradorRHC/INCA") +
  tema_pub
salvar(p_17, "17_tempo_trat_regiao", w = 12, h = 7)

# Violino: tempo diag-trat por sítio
p_18 <- rhc_surv |>
  filter(sitio %in% sitios_analise, !is.na(tempo_diag_trat),
         tempo_diag_trat >= 0, tempo_diag_trat <= 365) |>
  ggplot(aes(x = sitio, y = tempo_diag_trat, fill = sitio)) +
  geom_violin(alpha = 0.7, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.14, fill = "white", color = "grey30",
               outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 21,
               size = 3, fill = "white", color = "black") +
  geom_hline(yintercept = 60, linetype = "dashed",
             color = "firebrick", linewidth = 0.8) +
  scale_fill_manual(values = PALETA_SITIO, guide = "none") +
  scale_y_continuous(limits = c(0,365), breaks = seq(0,365,60)) +
  scale_x_discrete(labels = function(x) str_wrap(x, 10)) +
  labs(title = "Tempo Diagnostico-Tratamento por Sitio Anatomico",
       subtitle = "HNSCC/SUS 2010-2023",
       x = NULL, y = "Dias ate inicio do tratamento",
       caption = "Excluidos valores >365 dias. Fonte: IntegradorRHC/INCA") +
  tema_pub
salvar(p_18, "18_tempo_trat_sitio", w = 13, h = 7)


# ══════════════════════════════════════════════════════════════════════════════
# ── 9. MAPAS COROPLÉTICOS ─────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

estados_br <- read_state(year = 2020, simplified = TRUE) |>
  mutate(abbrev_state = str_trim(as.character(abbrev_state)))

# Helper mapa
fazer_mapa <- function(df_join, fill_var, titulo, subtitulo, caption,
                        scale_fn, nome, w = 12, h = 9) {
  p <- estados_br |>
    left_join(df_join, by = "abbrev_state") |>
    ggplot() +
    geom_sf(aes(fill = .data[[fill_var]]), color = "white", linewidth = 0.3) +
    geom_sf_text(aes(label = abbrev_state), size = 2.2, color = "white",
                 fontface = "bold", check_overlap = TRUE) +
    scale_fn +
    labs(title = titulo, subtitle = subtitulo, caption = caption) +
    tema_mapa
  salvar(p, nome, w = w, h = h)
  p
}

# Casos por UF
df_casos <- rhc_surv |>
  filter(!is.na(uf_res_clean)) |>
  count(abbrev_state = uf_res_clean, name = "n_casos")

p_19 <- fazer_mapa(df_casos, "n_casos",
  "Distribuicao de Casos HNSCC por Estado",
  "IntegradorRHC/SUS 2010-2023 • Por UF de residencia",
  "Fonte: IntegradorRHC/INCA",
  scale_fill_viridis_c(name = "N casos", labels = comma_format(),
                        option = "plasma", na.value = "grey90", direction = -1),
  "19_mapa_casos_uf")

# Sobrevida em 3 anos por UF
df_surv3 <- rhc_surv |>
  filter(!is.na(uf_res_clean)) |>
  group_by(uf_res_clean) |> filter(n() >= 50) |>
  group_modify(~ {
    fit <- survfit(Surv(tempo_anos, status) ~ 1, data = .x)
    s   <- summary(fit, times = 3, extend = TRUE)
    tibble(surv3 = s$surv, n = nrow(.x))
  }) |> ungroup() |> rename(abbrev_state = uf_res_clean)

p_20 <- fazer_mapa(df_surv3, "surv3",
  "Sobrevida Global em 3 Anos por Estado — HNSCC/SUS",
  "KM por UF de residencia • Estados com >=50 casos",
  "Cinza = dados insuficientes. Fonte: IntegradorRHC/INCA",
  scale_fill_gradientn(name = "Sobrevida\nem 3 anos",
    colors = c("#6E2222","#C0392B","#E6A817","#27AE60"),
    limits = c(0.1, 0.8), labels = percent_format(accuracy = 1),
    na.value = "grey90"),
  "20_mapa_sobrevida_3anos")

# % com RT por UF
df_rt <- rhc_surv |>
  filter(!is.na(uf_res_clean)) |>
  group_by(abbrev_state = uf_res_clean) |> filter(n() >= 50) |>
  summarise(pct_rt = mean(tem_rt, na.rm = TRUE))

p_21 <- fazer_mapa(df_rt, "pct_rt",
  "Proporcao com Radioterapia por Estado",
  "HNSCC/SUS 2010-2023",
  "Cinza = dados insuficientes. Fonte: IntegradorRHC/INCA",
  scale_fill_viridis_c(name = "% com RT", labels = percent_format(accuracy = 1),
                        option = "cividis", na.value = "grey90"),
  "21_mapa_acesso_rt")

# Tempo mediano diag-trat por UF
df_tempo <- rhc_surv |>
  filter(!is.na(uf_res_clean), !is.na(tempo_diag_trat),
         tempo_diag_trat >= 0, tempo_diag_trat <= 365) |>
  group_by(abbrev_state = uf_res_clean) |> filter(n() >= 50) |>
  summarise(mediana_dias = median(tempo_diag_trat))

p_22 <- fazer_mapa(df_tempo, "mediana_dias",
  "Tempo Mediano Diagnostico-Tratamento por Estado",
  "Verde = rapido (<30 dias); Vermelho = demorado (>90 dias)",
  "Cinza = dados insuficientes. Fonte: IntegradorRHC/INCA",
  scale_fill_gradientn(name = "Mediana\n(dias)",
    colors = c("#27AE60","#E6A817","#C0392B","#6E2222"),
    breaks = c(20, 45, 75, 110), limits = c(0, 130), na.value = "grey90"),
  "22_mapa_tempo_tratamento")

# Figura composta 2×2
old <- theme_get(); theme_set(theme_bw(base_size = 12))
fig_disp <- (p_20 | p_21) / (p_22 | p_13) +
  plot_annotation(
    title   = "Desigualdades Regionais no Manejo de HNSCC no SUS",
    caption = "A: Sobrevida 3 anos | B: % com RT | C: Tempo diag-trat | D: Modalidade por regiao\nFonte: IntegradorRHC/INCA 2010-2023",
    tag_levels = "A",
    theme = theme(plot.title   = element_text(size = 15, face = "bold"),
                  plot.caption = element_text(size = 9, color = "grey50"))
  )
salvar(fig_disp, "23_figura_disparidade_composta", w = 18, h = 14)
theme_set(old)


# ══════════════════════════════════════════════════════════════════════════════
# ── 10. REGRESSÃO DE COX MULTIVARIADA ────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

rhc_cox <- rhc_surv |>
  filter(estadio_grupo != "Desconhecido", !is.na(sexo_label),
         !is.na(regiao), trat_grupo != "Sem informacao",
         !is.na(tabaco_label), !is.na(alcool_label)) |>
  droplevels()

cat("\nRegistros no modelo de Cox:", nrow(rhc_cox), "\n")

cox_multi <- coxph(
  Surv(tempo_anos, status) ~
    sexo_label + faixa_etaria + raca_grupo + regiao +
    estadio_grupo + sitio + tabaco_label + alcool_label +
    trat_grupo + tempo_diag_trat_cat + strata(coorte),
  data = rhc_cox
)

cat("\n── Teste de proporcionalidade (Schoenfeld) ──\n")
print(cox.zph(cox_multi))

# Tabela formatada
tbl_cox <- tbl_regression(
  cox_multi, exponentiate = TRUE,
  label = list(
    sexo_label          ~ "Sexo",
    faixa_etaria        ~ "Faixa etaria",
    raca_grupo          ~ "Raca/Cor",
    regiao              ~ "Regiao de residencia",
    estadio_grupo       ~ "Estadiamento",
    sitio               ~ "Sitio anatomico",
    tabaco_label        ~ "Tabagismo",
    alcool_label        ~ "Etilismo",
    trat_grupo          ~ "Modalidade de tratamento",
    tempo_diag_trat_cat ~ "Tempo diagnostico-tratamento"
  )
) |>
  bold_p(t = 0.05) |> bold_labels() |>
  modify_caption("**Tabela 2.** Preditores independentes de mortalidade — Cox multivariado")

print(tbl_cox)

# Forest plot — ggsave direto (ggforest retorna gtable, não ggplot)
p_forest <- ggforest(cox_multi, data = rhc_cox,
                     main = "Preditores de Mortalidade — HNSCC/SUS",
                     fontsize = 0.82, refLabel = "Referencia", noDigits = 2)
ggsave("graficos/24_cox_forest_plot.png", p_forest,
       width = 16, height = 14, dpi = 300, bg = "white")
message("  ✓ 24_cox_forest_plot")


# ══════════════════════════════════════════════════════════════════════════════
# ── 11. SALVAR BANCO ANALÍTICO ────────────────────────────────════════════════
# ══════════════════════════════════════════════════════════════════════════════

saveRDS(rhc_surv, "rhc_hnscc_analitico.rds")
write_csv(rhc_surv, "rhc_hnscc_analitico.csv")

cat("\n\n=======================================================\n")
cat(" Pipeline HNSCC v4.0 concluido!\n")
cat("=======================================================\n")
cat(" Registros:  ", nrow(rhc_surv), "\n")
cat(" Obitos:     ", sum(rhc_surv$status),
    "(", round(mean(rhc_surv$status)*100,1), "%)\n")
cat(" Graficos:    ./graficos/\n\n")
cat(" GLOBAIS (00-09)\n")
cat("  00  Fluxo referenciamento inter-regional\n")
cat("  01  KM tratamento\n")
cat("  02  KM estadiamento\n")
cat("  03  KM sitio (todos juntos)\n")
cat("  04  KM sitio (paineis 1 vs. resto)\n")
cat("  05  KM regiao\n")
cat("  06  KM tratamento — estadio III/IV\n")
cat("  07  KM raca/cor\n")
cat("  08  KM coorte (xlim=3)\n")
cat("  09a KM tempo diagnostico-tratamento\n")
cat("  09b KM tabagismo\n")
cat("  09c KM etilismo\n")
cat("  09d Painel exposicao x tratamento (global)\n\n")
cat(" POR SITIO — 5 sitios x 7 analises = ate 35 graficos (10a-10g)\n")
cat("  10a KM por tratamento\n")
cat("  10b KM por regiao\n")
cat("  10c KM por estadiamento\n")
cat("  10d KM por tempo diag-trat\n")
cat("  10e KM por tabagismo\n")
cat("  10f KM por etilismo\n")
cat("  10g Painel exposicao x tratamento (2x3 por sitio)\n\n")
cat(" PAINEIS COMPOSTOS (11-12)\n")
cat("  11  Tratamento por sitio\n")
cat("  12  Regiao por sitio\n")
cat("  12b Tabagismo por sitio\n")
cat("  12c Etilismo por sitio\n\n")
cat(" DISPARIDADE (13-23)\n")
cat("  13  Trat por regiao (barras)\n")
cat("  14  Trat por sitio (barras)\n")
cat("  15  Heatmap trat x sitio\n")
cat("  16  Acesso RT por regiao\n")
cat("  17  Tempo diag-trat por regiao (violino)\n")
cat("  18  Tempo diag-trat por sitio (violino)\n")
cat("  19  Mapa: casos por UF\n")
cat("  20  Mapa: sobrevida 3 anos por UF\n")
cat("  21  Mapa: acesso RT por UF\n")
cat("  22  Mapa: tempo diag-trat por UF\n")
cat("  23  Figura composta disparidade regional\n\n")
cat(" COX (24)\n")
cat("  24  Forest plot multivariado\n")
cat("=======================================================\n")

# ============================================================================
#  FIM DO SCRIPT v4.0
# ============================================================================
