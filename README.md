# Análise de Sobrevida — HNSCC no SUS (2010–2023)

> Pipeline completo (v4.0) para análise de sobrevida em **carcinomas epidermóides de cabeça e pescoço (HNSCC)** no **SUS**, a partir dos dados do **IntegradorRHC (INCA/DATASUS)**.

**Última atualização:** março/2026  
**Script principal:** `pipeline_hnscc_v4.R`

---

## Visão geral

Este repositório contém um *pipeline* em R focado em qualidade editorial (figuras prontas para publicação), que:

- importa e consolida arquivos `.dbf` do IntegradorRHC;
- filtra a coorte de HNSCC epidermóide (2010–2023);
- constrói variáveis clínicas, demográficas e de acesso ao tratamento;
- gera curvas de Kaplan–Meier globais e estratificadas;
- quantifica disparidades regionais (fluxos, barras, violinos, heatmaps);
- produz mapas coropléticos por UF (casos, sobrevida 3 anos, acesso à RT, tempo diag→trat);
- ajusta modelo de Cox multivariado (com `strata(coorte)`) e exporta *forest plot*;
- salva banco analítico final em `.rds` e `.csv`.

---

## Conteúdo gerado (saídas)

Ao rodar o script, serão criados:

- Pasta `./graficos/` com **figuras em PNG** (alta resolução);
- `rhc_hnscc_analitico.rds` (objeto R com a coorte analítica);
- `rhc_hnscc_analitico.csv` (coorte analítica em CSV).

### Lista de gráficos (nomes de arquivo)

**Globais (00–09)**
- `00_fluxo_referenciamento`
- `01_km_tratamento`
- `02_km_estadiamento`
- `03_km_sitio_todos`
- `04_km_sitio_paineis`
- `05_km_regiao`
- `06_km_avancado_trat`
- `07_km_raca`
- `08_km_coorte`
- `09a_km_tempo_diagnostico`
- `09b_km_tabagismo`
- `09c_km_etilismo`
- `09d_painel_exposicao_tratamento`

**Por sítio anatômico (10a–10g)**
Para cada um dos 5 sítios (Cavidade Oral, Orofaringe, Laringe, Hipofaringe, Nasofaringe):
- KM por tratamento, região, estadiamento, tempo diag→trat, tabagismo, etilismo;
- painel 2×3 de **Exposição × Tratamento**.

**Painéis compostos (11–12c)**
- `11_painel_tratamento_por_sitio`
- `12_painel_regiao_por_sitio`
- `12b_painel_tabaco_por_sitio`
- `12c_painel_alcool_por_sitio`

**Disparidade / mapas (13–23)**
- `13_trat_por_regiao`
- `14_trat_por_sitio`
- `15_heatmap_trat_sitio`
- `16_acesso_rt_regiao`
- `17_tempo_trat_regiao`
- `18_tempo_trat_sitio`
- `19_mapa_casos_uf`
- `20_mapa_sobrevida_3anos`
- `21_mapa_acesso_rt`
- `22_mapa_tempo_tratamento`
- `23_figura_disparidade_composta`

**Cox (24)**
- `24_cox_forest_plot`

---

## Como rodar

### 1) Requisitos

- R (recomendado: R >= 4.2)
- Pacotes (instalação rápida):

```r
install.packages(c(
  "foreign","tidyverse","lubridate","janitor",
  "survival","survminer","gtsummary",
  "ggplot2","patchwork","viridis","scales","broom",
  "geobr","sf","ggtext"
))
```

### 2) Ajuste do caminho dos dados

No início do script, edite a variável:

```r
PASTA_DADOS <- "~/Downloads/SEU_DIRETORIO_AQUI/"
```

O script espera encontrar arquivos com padrão:

- `rhc\d+\.dbf` (ex.: `rhc2010.dbf`, `rhc2011.dbf`, ...)

### 3) Execute

No R/RStudio:

```r
source("pipeline_hnscc_v4.R")
```

---

## Principais definições (documentação rápida)

### Coorte
- **Período:** 2010–2023 (corte administrativo em `2023-12-31`)
- **HNSCC (CID):** `C00–C14`, `C30`, `C31`, `C32`
- **Histologia (epidermóide):** `8070–8077`

### Desfecho e tempo
- `status`: 1 = óbito registrado, 0 = censura no corte (`2023-12-31`)
- `tempo_anos`: tempo de seguimento em anos (dias/365.25)
- `tempo_diag_trat`: dias entre início (consulta/diagnóstico) e início de RT (`datainitrt`)

### Tratamento (grupo)
Derivado de `pritrath`:
- Cirurgia (+/- adj)
- QT+RT
- RT exclusiva
- QT exclusiva
- Nenhum

---

## Arquitetura de figuras (por que o script é robusto)

O pipeline separa intencionalmente temas para evitar conflitos entre **ggplot2** e **survminer**:

- `tema_pub`: para *ggplot2 puro* (usa `element_markdown` via **ggtext**) 
- `tema_surv`: para objetos `ggsurvplot` (*sem* `element_markdown`)
- `salvar()`: salva ggplot (`ggsave`)
- `salvar_km()`: salva `ggsurvplot` usando dispositivo `png()` (inclui *risk table*)
- `salvar_painel()`: painéis via `patchwork::wrap_plots()` com *guard* do tema global

---

## Observações importantes

- O script foi escrito para rodar localmente; ele **não baixa** dados automaticamente.
- Algumas variáveis e códigos (ex.: tratamento, tabagismo, etilismo, estadiamento) dependem do dicionário do IntegradorRHC.
- Mapas por UF aplicam filtro mínimo **n >= 50** para estabilidade de estimativas.

---

## Citação / fonte

- **Fonte dos dados:** IntegradorRHC (INCA/DATASUS)
- **Uso científico:** verifique regras de uso/citação do INCA/DATASUS para publicação.

---

## Autor

- GitHub: @pasconfilho