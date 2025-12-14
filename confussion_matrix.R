# =========================
# Comparación KIR: TIPATGES vs RESULTADO
# Confusion matrix + gráficos
# =========================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(janitor)
  library(pheatmap)
})

# ---- RUTAS ----
tip_file <- "TIPATGES CONTROLS MÈTODE.csv"
res_file <- "RESULTADO_CONSOLIDADO_KIR.tsv"

out_dir <- "salidas_comparacion_kir"
dir.create(out_dir, showWarnings = FALSE)

# ---- FUNCIONES AUXILIARES ----
normalize_sample <- function(x) {
  x %>%
    str_replace_all("-", "_") %>%      # AMAI-KIR -> AMAI_KIR
    str_replace("_S\\d+$", "") %>%     # AMAI_KIR_S4 -> AMAI_KIR
    str_trim()
}

normalize_gt <- function(x) {
  x <- ifelse(is.na(x), NA_character_, as.character(x))
  x <- str_trim(x)
  x[x == ""] <- NA_character_
  
  is_plus <- !is.na(x) & x == "+"
  
  y <- x
  idx <- !is.na(y) & !is_plus & str_detect(y, "/")
  y[idx] <- sapply(y[idx], function(z) {
    alleles <- unlist(str_split(z, "/"))
    alleles <- str_trim(alleles)
    alleles <- alleles[alleles != ""]
    paste(sort(alleles), collapse = "/")
  })
  
  y
}

present_absent <- function(x) {
  ifelse(is.na(x) | str_trim(as.character(x)) == "", 0L, 1L)
}

# ---- CARGA DE DATOS ----
tip_raw <- read_csv(tip_file, show_col_types = FALSE) %>% clean_names()
res_raw <- read_tsv(res_file, show_col_types = FALSE) %>% clean_names()

# Columna de muestra en TIPATGES (normalmente la primera)
tip_sample_col <- names(tip_raw)[1]

tip <- tip_raw %>%
  rename(sample = all_of(tip_sample_col)) %>%
  mutate(sample = normalize_sample(sample))

res <- res_raw %>%
  rename(sample = muestra) %>%
  mutate(sample = normalize_sample(sample))

# Genes comunes
genes_common <- intersect(
  setdiff(names(tip), "sample"),
  setdiff(names(res), "sample")
)
genes_common <- sort(genes_common)

if (length(genes_common) == 0) {
  stop("No hay genes comunes entre TIPATGES y RESULTADO.")
}

# ---- NORMALIZACIÓN DE GENOTIPOS ----
tip_n <- tip %>%
  select(sample, all_of(genes_common)) %>%
  mutate(across(all_of(genes_common), normalize_gt))

res_n <- res %>%
  select(sample, all_of(genes_common)) %>%
  mutate(across(all_of(genes_common), normalize_gt))

# ---- FIX CRÍTICO: forzar todo a character ----
tip_n <- tip_n %>% mutate(across(-sample, as.character))
res_n <- res_n %>% mutate(across(-sample, as.character))

# ---- FILTRAR MUESTRAS COMUNES ----
samples_common <- intersect(tip_n$sample, res_n$sample)

tip_n <- tip_n %>% filter(sample %in% samples_common) %>% arrange(sample)
res_n <- res_n %>% filter(sample %in% samples_common) %>% arrange(sample)

# ---- FORMATO LARGO ----
tip_long <- tip_n %>%
  pivot_longer(-sample, names_to = "gene", values_to = "truth_gt")

res_long <- res_n %>%
  pivot_longer(-sample, names_to = "gene", values_to = "pred_gt")

# ---- COMPARACIÓN ----
cmp <- tip_long %>%
  inner_join(res_long, by = c("sample", "gene")) %>%
  mutate(
    truth_pa = present_absent(truth_gt),
    pred_pa  = present_absent(pred_gt),
    exact_eligible = !is.na(truth_gt) & truth_gt != "+" & !is.na(pred_gt),
    exact_match = ifelse(exact_eligible, truth_gt == pred_gt, NA)
  )

# =========================
# CONFUSION MATRIX GLOBAL
# =========================
global_cm <- cmp %>%
  mutate(
    truth = factor(ifelse(truth_pa == 1, "Present", "Absent"),
                   levels = c("Absent", "Present")),
    pred  = factor(ifelse(pred_pa == 1, "Present", "Absent"),
                   levels = c("Absent", "Present"))
  ) %>%
  count(truth, pred) %>%
  complete(truth, pred, fill = list(n = 0))

write_csv(global_cm,
          file.path(out_dir, "confusion_matrix_global_presence_absence.csv"))

p_cm <- ggplot(global_cm, aes(x = pred, y = truth, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), size = 6) +
  labs(
    title = "Confusion matrix global (Presencia/Ausencia)",
    x = "Predicción (RESULTADO)",
    y = "Verdad (TIPATGES)"
  ) +
  theme_minimal()

ggsave(file.path(out_dir,
                 "confusion_matrix_global_presence_absence.png"),
       p_cm, width = 7, height = 5, dpi = 200)

# =========================
# MÉTRICAS POR GEN
# =========================
metrics_by_gene <- cmp %>%
  group_by(gene) %>%
  summarise(
    TP = sum(truth_pa == 1 & pred_pa == 1, na.rm = TRUE),
    TN = sum(truth_pa == 0 & pred_pa == 0, na.rm = TRUE),
    FP = sum(truth_pa == 0 & pred_pa == 1, na.rm = TRUE),
    FN = sum(truth_pa == 1 & pred_pa == 0, na.rm = TRUE),
    accuracy    = (TP + TN) / (TP + TN + FP + FN),
    sensitivity = ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_),
    specificity = ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_),
    ppv         = ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_),
    .groups = "drop"
  ) %>%
  arrange(desc(accuracy))

write_csv(metrics_by_gene,
          file.path(out_dir, "metrics_by_gene_presence_absence.csv"))

metrics_long <- metrics_by_gene %>%
  select(gene, accuracy, sensitivity, specificity, ppv) %>%
  pivot_longer(-gene, names_to = "metric", values_to = "value")

p_metrics <- ggplot(metrics_long,
                    aes(x = reorder(gene, value), y = value)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ metric, scales = "free_x") +
  labs(title = "Métricas por gen (Presencia/Ausencia)",
       x = "Gen", y = "Valor") +
  theme_minimal()

ggsave(file.path(out_dir,
                 "metrics_by_gene_presence_absence.png"),
       p_metrics, width = 10, height = 7, dpi = 200)

# =========================
# MÉTRICAS POR MUESTRA
# =========================
sample_acc <- cmp %>%
  group_by(sample) %>%
  summarise(
    n = n(),
    acc = mean(truth_pa == pred_pa, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(acc))

write_csv(sample_acc,
          file.path(out_dir, "accuracy_by_sample_presence_absence.csv"))

p_sample <- ggplot(sample_acc,
                   aes(x = reorder(sample, acc), y = acc)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Accuracy por muestra (Presencia/Ausencia)",
       x = "Muestra", y = "Accuracy") +
  theme_minimal()

ggsave(file.path(out_dir,
                 "accuracy_by_sample_presence_absence.png"),
       p_sample, width = 9, height = 8, dpi = 200)

# =========================
# HEATMAP DE DISCREPANCIAS
# =========================
disc_mat <- cmp %>%
  mutate(disagree = as.integer(truth_pa != pred_pa)) %>%
  select(sample, gene, disagree) %>%
  pivot_wider(names_from = gene,
              values_from = disagree,
              values_fill = 0) %>%
  arrange(sample)

disc_matrix <- disc_mat %>%
  select(-sample) %>%
  as.data.frame()

rownames(disc_matrix) <- disc_mat$sample

png(file.path(out_dir,
              "heatmap_discrepancias_presence_absence.png"),
    width = 1200, height = 900)

pheatmap(as.matrix(disc_matrix),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Heatmap discrepancias (1 = distinto)")

dev.off()

# =========================
# EXACT MATCH (opcional)
# =========================
exact_summary <- cmp %>%
  filter(exact_eligible) %>%
  summarise(
    n = n(),
    exact_acc = mean(exact_match, na.rm = TRUE)
  )

write_csv(exact_summary,
          file.path(out_dir, "exact_match_summary.csv"))

exact_by_gene <- cmp %>%
  filter(exact_eligible) %>%
  group_by(gene) %>%
  summarise(
    n = n(),
    exact_acc = mean(exact_match, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(exact_acc))

write_csv(exact_by_gene,
          file.path(out_dir, "exact_match_by_gene.csv"))

p_exact <- ggplot(exact_by_gene,
                  aes(x = reorder(gene, exact_acc), y = exact_acc)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Exact match por gen (excluye '+')",
       x = "Gen", y = "Exact match accuracy") +
  theme_minimal()

ggsave(file.path(out_dir,
                 "exact_match_by_gene.png"),
       p_exact, width = 9, height = 6, dpi = 200)

# ---- FIN ----
cat(
  "\nANÁLISIS COMPLETADO\n",
  "Directorio:", out_dir, "\n",
  "Genes comparados:", length(genes_common), "\n",
  "Muestras comparadas:", length(samples_common), "\n"
)
