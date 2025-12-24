library(MotrpacRatTraining6moMuscleData)
library(dplyr)
library(tidyr)


DA_objects = list(
  METAB_GN = METAB_GN_DA[["trained_vs_SED"]],
  METAB_VL = METAB_VL_DA[["trained_vs_SED"]],
  PROT_GN = PROT_GN_DA[["trained_vs_SED"]],
  TRNSCRPT_GN = TRNSCRPT_GN_DA[["trained_vs_SED"]],
  TRNSCRPT_VL = TRNSCRPT_VL_DA[["trained_vs_SED"]]
)

count_df =  lapply(names(DA_objects), function(obj) {
  DA_objects[[obj]] %>%
    filter(adj.P.Val < 0.05) %>%
    mutate(
      ome = sub("_(GN|VL)$", "", obj),
      GN_VL = sub("^.*_", "", obj)
    ) %>%
    count(contrast, ome, GN_VL)
}) %>%
  bind_rows() %>%
  mutate(
    sex = ifelse(grepl("^M", contrast), "Male", "Female"),
    comparison = gsub("M_|F_", "", contrast),
    ome = factor(ome,
                 levels = c("PROT", "TRNSCRPT", "METAB"))
  ) %>%
  mutate(comparison = factor(comparison,
                             levels = c("8W - SED",
                                        "4W - SED",
                                        "2W - SED",
                                        "1W - SED")),
         GN_VL = factor(GN_VL,
                        levels = c("VL", "GN"))
  )

p = ggplot(
  count_df,
  aes(
    x = n,
    y = comparison,
    fill = GN_VL
  )
) +
  geom_col(position = position_dodge(width = 0.5), width = 0.4) +
  facet_grid(
    sex ~ ome,
    scales = "free_x"
  ) +
  scale_fill_manual(
    name = "Skeletal Muscle",
    values = c(GN = "#8B0000", VL = "grey60"),
    labels = c(GN = "GN", VL = "VL")
  ) +
  labs(
    x = "N Significant Features (FDR = 0.05)",
    y = NULL,
    fill = NULL
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90", colour = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )



pdf(file = file.path(
  here(),
  "plots",
  "S1B.pdf"),
  width = 8,
  height = 8
)
print(p)
dev.off()

