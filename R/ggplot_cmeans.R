library(ggplot2)
library(dplyr)
library(ggpubr)
library(grid)

## Iridescent color scheme from https://personal.sron.nl/~pault/
palette <- c("#fefbe9", "#fcf7d5", "#f5f3c1", "#eaf0b5", "#ddecbf",
             "#d0e7ca", "#c2e3d2", "#b5ddd8", "#a8d8dc", "#9bd2e1",
             "#8dcbe4", "#81c4e7", "#7bbce7", "#7eb2e4", "#88a5dd",
             "#9398d2", "#9b8ac4", "#9d7db2", "#9a709e", "#906388",
             "#805770", "#684957", "#46353a")


ggplot_cmeans <- function(cl,
                          zmat,
                          keep_clusters = as.integer(rownames(cl$centers)),
                          min.membership = 0,
                          common_ylim = FALSE,
                          ncol = 6L,
                          filename,
                          dpi = 200) {

  cluster_specific_sex = ifelse(grep("F", names(zmat)))

  mem <- cl$membership[, keep_clusters, drop = FALSE]
  centers <- cl$centers[keep_clusters, , drop = FALSE]
  colnames(mem) <- rownames(centers) <- paste0("Cluster ", colnames(mem))
  colnames(zmat) <- colnames(centers) <- paste0("timepoint_", colnames(zmat))

  # Include SED timepoints for both sexes
  sed <- matrix(data = 0, nrow = nrow(zmat), ncol = 2L,
                   dimnames = list(rownames(zmat),
                                   paste0("timepoint_",
                                          c("F", "M"),
                                          "_SED")))

  zmat <- cbind(sed, zmat)
  zmat <- zmat[, c(grep("F_", colnames(zmat)),
                   grep("M_", colnames(zmat)))]

  sed_center <- matrix(data = 0, nrow = nrow(centers), ncol = 2L)
  rownames(sed_center) <- rownames(centers)
  colnames(sed_center) <- colnames(sed)

  centers <- cbind(centers, sed_center)
  centers <- centers[, colnames(zmat)]

  cluster_assignment <- data.frame(feature = names(cl$cluster),
                                   cluster = paste0("Cluster ",
                                                    cl$cluster))

  z_df <- zmat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature") %>%
    tidyr::pivot_longer(cols = -feature,
                        names_to = "timepoint",
                        names_pattern = "timepoint_(.*)",
                        values_to = "z") %>%
    mutate(sex = ifelse(sub("_.*", "", timepoint) == "F", "Female", "Male"),
           timepoint = sub("[MF]_", "", timepoint),
           across(.cols = c(sex, timepoint),
                  .fns = ~ factor(.x, levels = unique(.x))))

  cluster_df <- mem %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature") %>%
    tidyr::pivot_longer(cols = -feature,
                        names_to = "cluster",
                        values_to = "membership") %>%
    filter(membership >= min.membership) %>%
    inner_join(cluster_assignment, by = c("feature", "cluster")) %>%
    mutate(cluster = factor(cluster, levels = rownames(centers)))

  df <- inner_join(z_df, cluster_df, by = "feature")  %>%
    arrange(cluster, membership) %>%
    mutate(feature = factor(feature, levels = unique(feature))) %>%
    droplevels.data.frame()

  # Cluster centroids
  center_df <- centers %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cluster") %>%
    tidyr::pivot_longer(cols = -cluster,
                        names_to = "timepoint",
                        names_pattern = "timepoint_(.*)",
                        values_to = "center") %>%
    mutate(sex = ifelse(sub("_.*", "", timepoint) == "F", "Female", "Male"),
           timepoint = sub("[MF]_", "", timepoint),
           across(.cols = c(sex, timepoint),
                  .fns = ~ factor(.x, levels = unique(.x))),
           cluster = factor(cluster, levels = rownames(centers)))

  if (common_ylim)
    ylim <- range(df$z)
  else
    ylim <- rep(NA, 2L)

  df_list <- tidyr::nest(df, .by = cluster)
  center_df_list <- tidyr::nest(center_df, .by = cluster)

  plotlist <- lapply(seq_len(nrow(df_list)), function(i) {
    df_i <- df_list$data[[i]]
    centroid_i <- center_df_list$data[[i]]

    ggplot(df_i, aes(x = timepoint, y = z)) +
      geom_line(aes(color = membership, group = feature)) +
      geom_line(aes(x = timepoint, y = center, group = 1),
                data = centroid_i, color = "black") +
      facet_grid(~ sex) +
      scale_x_discrete(name = NULL,
                       expand = expansion(add = 0.2)) +
      scale_y_continuous(limits = ylim) +
      scale_color_gradientn(name = "Membership\nProbability",
                            colors = palette,
                            values = seq(0, 1, length.out = length(palette)),
                            limits = c(0, 1),
                            breaks = seq(0, 1, 0.2)) +
      labs(y = NULL,
           title = df_list$cluster[i]) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            strip.text = element_text(face = "bold",
                                      margin = margin(t = 3)),
            plot.title = element_text(size = rel(0.9),
                                      face = "bold", hjust = 0.5,
                                      margin = margin(b = 0)),
            strip.background = element_blank(),
            axis.title.y = element_text(color = "black", face = "bold",
                                        margin = margin(r = 8)),
            axis.text.x = element_text(size = rel(0.7),
                                       color = "black", angle = 90,
                                       vjust = 0.5, hjust = 1),
            axis.text.y = element_text(color = "black"))
  })

  ncol <- min(length(plotlist), ncol)
  nrow <- ceiling(length(plotlist) / ncol)
  height <- nrow * 2
  width <- ncol * 2 + 0.5 # 0.5 = space for legend

  p <- ggarrange(plotlist = plotlist,
                 ncol = ncol,
                 nrow = nrow,
                 common.legend = TRUE,
                 legend = "right")

  # Add global y-axis title
  p <- annotate_figure(p = p,
                       left = textGrob("Z-Score", rot = 90, vjust = 1,
                                       gp = gpar(fontsize = 11,
                                                 fontface = "bold")))

  ggsave(filename = filename, plot = p,
         height = height, width = width,
         dpi = dpi, bg = "white")
}
