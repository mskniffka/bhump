# Figure specifications
#
# 2020-05-15
#
# Jonas Schöley

# figure specs
fig_spec <- list()

# ggplot theme ----------------------------------------------------

# ggplot theme
fig_spec$MyGGplotTheme <-
  function (
    size = 8,
    family = 'serif',
    scaler = 1,
    no_axes = FALSE,
    panel_border = FALSE,
    ar = NA
  ) {
    
    size_med = size*scaler
    size_sml = round(size*0.7)*scaler
    
    list(
      theme_classic(base_size = size_med, base_family = family),
      theme(
        # basic
        text = element_text(color = 'black'),
        line = element_line(size = 0.3*scaler, lineend = 'square'),
        # axis
        axis.title = element_text(size = size_med, face = 'bold'),
        axis.ticks = element_line(size = rel(0.5), color = 'black'),
        axis.text = element_text(size = size_med, color = 'black'),
        # strips
        strip.text = element_text(color = 'black', size = size_med),
        strip.background = element_blank(),
        # plot
        title = element_text(face = 'bold'),
        plot.subtitle = element_text(color = 'black', size = size_med, face = 'bold'),
        plot.caption = element_text(color = 'black', size = size_sml, face = 'plain'),
        plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(1, 0, 0, 0.5), units = 'mm'),
      ),
      if (isTRUE(panel_border)) {
        theme(
          panel.border =
            element_rect(fill = NA)
        )
      },
      if (isTRUE(no_axes)) {
        theme(
          axis.line = element_blank()
        )
      },
      if (!is.na(ar)) {
        theme(
          aspect.ratio = ar
        )
      }
    )
  }

# Dimensions ------------------------------------------------------

# figure width (mm)
fig_spec$width = 124.6

# figure line and point sizes
fig_spec$line_size_xs = 0.1
fig_spec$line_size_s = 0.2
fig_spec$line_size_m = 0.3
fig_spec$line_size_l = 0.5
fig_spec$line_size_xl = 1
fig_spec$point_size_xs = 0.05
fig_spec$point_size_s = 0.1
fig_spec$point_size_m = 0.2
fig_spec$stroke_size_xs = 0.3
fig_spec$text_size_s = 2
fig_spec$text_size_xs = 1.5

# Colors ----------------------------------------------------------

# color palette
fig_spec$discrete_colors <-
  c('#D23737', # red
    '#3191C9', # blue
    '#D2BC2D', # yellow
    '#4EC93B', # green
    '#881F93', # purple
    '#C5752B') # orange
fig_spec$discrete_colors_light <-
  c('#FCB3B3', # red
    '#A7DDFC', # blue
    '#FAEC8E'  # yellow
    )

# Export function -------------------------------------------------

fig_spec$ExportPDF <-
  function (figure, filename, path, ...) {
    ggsave(
      filename = paste0(filename, '.pdf'),
      plot = figure,
      path = path,
      units = 'mm',
      dpi = 300,
      useDingbats = FALSE,
      ...
    )
  }

