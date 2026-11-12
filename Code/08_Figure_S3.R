library(ggplot2)
rs<-readRDS("Data/Processed objects/rank_shift.rds")

plant_family <- c(
  Fragaria_vesca            = "Rosaceae",
  Silene_dioica             = "Caryophyllaceae",
  Campanula_sp              = "Campanulaceae",
  Ajuga_reptans             = "Lamiaceae",
  Angelica_sylvestris       = "Apiaceae",
  Scorzoneroides_autumnalis = "Asteraceae",
  Stachys_alpina            = "Lamiaceae",
  Trifolium_pratense        = "Fabaceae",
  Alliaria_petiolata        = "Brassicaceae",
  Aegopodium_podagraria     = "Apiaceae",
  Persicaria_bistorta       = "Polygonaceae",
  Daucus_carota             = "Apiaceae",
  Cirsium_palustre          = "Asteraceae",
  Lotus_corniculatus        = "Fabaceae",
  Dactylorhiza_maculata     = "Orchidaceae",
  Cirsium_arvense           = "Asteraceae",
  Chaerophyllum_hirsutum    = "Apiaceae",
  Galium_mollugo            = "Rubiaceae",
  Geranium_sylvaticum       = "Geraniaceae",
  Plantago_lanceolata       = "Plantaginaceae",
  Eupatorium_cannabinum     = "Asteraceae",
  Aquilegia_atrata          = "Ranunculaceae",
  Cirsium_oleraceum         = "Asteraceae",
  Trifolium_repens          = "Fabaceae",
  Myosotis_sp               = "Boraginaceae",
  Syringa_vulgaris          = "Oleaceae",
  Aruncus_dioicus           = "Rosaceae",
  Erigeron_annuus           = "Asteraceae",
  Cirsium_vulgare           = "Asteraceae",
  Arctium_minus             = "Asteraceae",
  Valeriana_officinalis     = "Caprifoliaceae",
  Hypochaeris_radicata      = "Asteraceae",
  Taraxacum_officinale      = "Asteraceae",
  Jacobaea_vulgaris         = "Asteraceae",
  Achillea_millefolium      = "Asteraceae",
  Cirsium_eriophorum        = "Asteraceae",
  Anthriscus_sylvestris     = "Apiaceae",
  Filipendula_ulmaria       = "Rosaceae",
  Knautia_dipsacifolia      = "Caprifoliaceae",
  Ranunculus_sp             = "Ranunculaceae",
  Hippocrepis_comosa        = "Fabaceae",
  Heracleum_sphondylium     = "Apiaceae",
  Symphoricarpos_albus      = "Caprifoliaceae",
  Epilobium_angustifolium   = "Onagraceae",
  Stachys_sylvatica         = "Lamiaceae",
  Oenothera_biennis         = "Onagraceae"
)

rs$family <- unname(plant_family[rs$plant])
stopifnot(!anyNA(rs$family)) 

plot_rank_dumbbells <- function(rs,
                                col_aft   = "#F2A541",   # afternoon
                                col_night = "#6E6AA8",   # night
                                xlab = paste0("Rank mismatch, controlled for plant abundance\n",
                                              "(negative: pollinator side; positive: herbivore side)"),
                                xlim = NULL,             # NULL = symmetric around 0
                                sort_by = c("min", "night", "shift", "afternoon")) {
  sort_by <- match.arg(sort_by)
  d <- rs[!is.na(rs$mean_shift), ]
  d$sort_value <- switch(sort_by,
                         min       = pmin(d$mean_d_aft, d$mean_d_night),  # leftmost dot
                         night     = d$mean_d_night,
                         afternoon = d$mean_d_aft,
                         shift     = d$mean_shift)
  if (is.null(xlim)) {
    m <- max(abs(c(d$mean_d_aft, d$mean_d_night)), na.rm = TRUE)
    xlim <- c(-m, m) * 1.05
  }
  d$family <- factor(d$family, levels = sort(unique(d$family)))
  
  # y order: ggplot draws the first level at the bottom, so within each family
  # sort by DEcreasing value -> most negative value ends up at the top
  d <- d[order(d$family, -d$sort_value), ]
  d$plant <- factor(d$plant, levels = unique(d$plant))
  
  # labels: "(n)  Genus species" with the name in italics; "sp" stays upright
  lab_expr <- vapply(seq_len(nrow(d)), function(i) {
    parts <- strsplit(as.character(d$plant[i]), "_")[[1]]
    name  <- if (tail(parts, 1) == "sp")
      sprintf("italic('%s')~'sp.'", parts[1]) else
        sprintf("italic('%s')", paste(parts, collapse = " "))
    sprintf("'(%d)'~~%s", d$n_comm_shared[i], name)
  }, character(1))
  names(lab_expr) <- as.character(d$plant)
  
  # long format for the two points
  pts <- rbind(
    data.frame(plant = d$plant, family = d$family, d = d$mean_d_aft,   period = "Afternoon"),
    data.frame(plant = d$plant, family = d$family, d = d$mean_d_night, period = "Night"))
  pts$period <- factor(pts$period, levels = c("Afternoon", "Night"))
  
  ggplot(d, aes(y = plant)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_segment(aes(x = mean_d_aft, xend = mean_d_night, yend = plant),
                 colour = "grey55", linewidth = 0.5) +
    geom_point(data = pts, aes(x = d, colour = period), size = 2.2) +
    scale_colour_manual(values = c(Afternoon = col_aft, Night = col_night), name = NULL) +
    scale_y_discrete(labels = function(b) parse(text = lab_expr[b])) +
    scale_x_continuous(limits = xlim) +
    facet_grid(family ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 10) +
    theme(strip.placement   = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
          strip.background  = element_blank(),
          panel.spacing.y   = unit(1.5, "mm"),
          panel.grid.minor  = element_blank(),
          panel.border      = element_rect(colour = "grey75"),
          panel.grid.major.y = element_blank(),
          legend.position   = "top",
          plot.margin       = margin(5, 25, 5, 5))
}

# USAGE
p <- plot_rank_dumbbells(rs)
ggsave("Figures/Figure S3.png", p, width = 18, height = 22, units = "cm")
