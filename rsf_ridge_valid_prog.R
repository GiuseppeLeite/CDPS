# ============================================================
# Apply CDPS to GSE91061 and generate Pre vs On-treatment plot
# ============================================================

library(tidyverse)
library(glmnet)


# 1) Load model and dataset --------------------------------------


model_rsf_ridge <- readRDS("model_rsf_ridge.rds")
GSE91061_mime <- readRDS("GSE91061_mime.rds")

GSE91061_mime <- as.data.frame(GSE91061_mime, check.names = FALSE)


# 2) Apply CDPS  --------------------------------------


genes_model <- rownames(model_rsf_ridge$glmnet.fit$beta)

if (length(genes_model) != 7) {
  stop("The model is expected to contain 7 genes, but found: ", length(genes_model))
}

missing_genes <- setdiff(genes_model, colnames(GSE91061_mime))

if (length(missing_genes) > 0) {
  stop("Missing model genes in GSE91061_mime: ", paste(missing_genes, collapse = ", "))
}

expr_mat <- GSE91061_mime[, genes_model, drop = FALSE]
expr_mat[] <- lapply(expr_mat, function(x) as.numeric(as.character(x)))
expr_mat <- as.matrix(expr_mat)

if (anyNA(expr_mat)) {
  stop("Missing or non-numeric values detected in the model gene expression matrix.")
}

expr_z <- scale(expr_mat, center = TRUE, scale = TRUE)

if (any(!is.finite(expr_z))) {
  stop("At least one model gene has zero variance after scaling.")
}

GSE91061_mime$CDPS <- as.numeric(
  predict(
    model_rsf_ridge$glmnet.fit,
    newx = expr_z,
    s = model_rsf_ridge$lambda.min,
    type = "link"
  )
)

saveRDS(GSE91061_mime, "GSE91061_mime_CDPS.rds")


# 3) Prepare data for Pre vs On-treatment analysis  --------------------------------------


df <- GSE91061_mime %>%
  mutate(
    Response_bin = case_when(
      Response == "PRCR" ~ "Responder",
      TRUE ~ "Non-responder"
    ),
    Response_bin = factor(Response_bin, levels = c("Responder", "Non-responder")),
    Group = factor(Group, levels = c("Pre", "On"))
  )

paired_ids <- df %>%
  filter(Group %in% c("Pre", "On")) %>%
  count(Patient, Group) %>%
  pivot_wider(names_from = Group, values_from = n) %>%
  filter(Pre == 1, On == 1) %>%
  pull(Patient)

df_wide <- df %>%
  filter(Patient %in% paired_ids, Group %in% c("Pre", "On")) %>%
  select(Patient, Response_bin, Group, CDPS) %>%
  pivot_wider(
    id_cols = c(Patient, Response_bin),
    names_from = Group,
    values_from = CDPS
  ) %>%
  drop_na(Pre, On)

df_long <- df_wide %>%
  pivot_longer(
    cols = c(Pre, On),
    names_to = "Group",
    values_to = "CDPS"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Pre", "On")),
    Group_label = recode(
      Group,
      "Pre" = "Pre-treatment",
      "On" = "On-treatment"
    ),
    Group_label = factor(
      Group_label,
      levels = c("Pre-treatment", "On-treatment")
    )
  )


# 4) Paired Wilcoxon tests  --------------------------------------


y_offset <- 0.08 * diff(range(df_long$CDPS, na.rm = TRUE))

pvals <- df_wide %>%
  group_by(Response_bin) %>%
  summarise(
    p = wilcox.test(Pre, On, paired = TRUE, exact = FALSE)$p.value,
    y = max(c(Pre, On), na.rm = TRUE) + y_offset,
    .groups = "drop"
  ) %>%
  mutate(
    x = 1.5,
    label = paste0("p = ", signif(p, 2))
  )


# 5) Plot --------------------------------------

p_pre_on_response <- ggplot(
  df_long,
  aes(x = Group_label, y = CDPS, fill = Group_label)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.85,
    color = "black"
  ) +
  geom_line(
    aes(group = Patient),
    color = "gray45",
    alpha = 0.55,
    linewidth = 0.4
  ) +
  geom_point(
    size = 1.2,
    alpha = 0.55,
    color = "black"
  ) +
  geom_text(
    data = pvals,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    fontface = "italic",
    size = 4.5
  ) +
  facet_wrap(
    ~ Response_bin,
    nrow = 1,
    strip.position = "bottom"
  ) +
  scale_fill_manual(
    name = "Groups",
    values = c(
      "Pre-treatment" = "#F2C47E",
      "On-treatment" = "#7EAFE8"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(face = "italic", size = 13),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "GSE91061",
    x = NULL,
    y = "CDPS"
  )

p_pre_on_response

ggsave(
  filename = "GSE91061_PRE_ON_ByResponse_CDPS.pdf",
  plot = p_pre_on_response,
  width = 5.5,
  height = 4
)
