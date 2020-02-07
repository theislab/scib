#' Visualizing scores in an overview table
#'
#' \code{scIB_knit_table} returns a ggplot2 object with a ranking of the integration methods  
#' and their relative scores over a set of metrics.
#'
#' @param data A data frame with methods by row, metrics by columns. First column must be method's names. Score values
#' should be scaled between [0,1].
#' @param column_info A data frame describing the columns of `data`. This data frame should contain the following columns:
#'   * `id` (`character`): The corresponding column name in `data`.
#'   * `color_group` (`character`): The group name of the corresponding palette of the column.
#'   * `geom` (`character`): The geom of the column. Must be one of: `circle`, `bar`, or `text`.
#'   * `width`: Custom width for this column (default: 1).
#'   * `overlay`: Whether to overlay this column over the previous column. If so, the width of that column will be inherited.
#' @param row_info A data frame describing the rows of `data`. This data should contain the following columns:`
#'   * `id` (`character`): The corresponding row name in `data`.
#'   * `group` (`character`): The group of the row. If all are `NA`, the rows will not be split up into groups.
#' @param palettes A list of palettes to be used in the respective `color_group`.
#'
#'
#' @importFrom ggforce geom_arc_bar geom_circle geom_arc
#' @importFrom cowplot theme_nothing
#'
#' @export
scIB_knit_table <- function(
data,
column_info,
row_info,
palettes,
usability = FALSE
) {
  # no point in making these into parameters
  row_height <- 1.1
  row_space <- .1
  row_bigspace <- .5
  col_width <- 1
  col_space <- .1
  col_bigspace <- .5
  segment_data <- NULL
  
  # DETERMINE ROW POSITIONS
  if (!"group" %in% colnames(row_info) || all(is.na(row_info$group))) {
    row_info$group <- ""
    row_groups <- tibble(group = "")
    plot_row_annotation <- FALSE
  } else {
    plot_row_annotation <- TRUE
  }
  
  row_pos <-
    row_info %>%
    group_by(group) %>%
    mutate(group_i = row_number()) %>%
    ungroup() %>%
    mutate(
      row_i = row_number(),
      colour_background = group_i %% 2 == 1,
      do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
      ysep = ifelse(do_spacing, row_height + 2 * row_space, row_space),
      y = - (row_i * row_height + cumsum(ysep)),
      ymin = y - row_height / 2,
      ymax = y + row_height / 2
    )
  
  # DETERMINE COLUMN POSITIONS
  if (!"group" %in% colnames(column_info) || all(is.na(column_info$group))) {
    column_info$group <- ""
    plot_column_annotation <- FALSE
  } else {
    plot_column_annotation <- TRUE
  }
  
  column_info <-
    column_info %>%
    add_column_if_missing(width = col_width, overlay = FALSE)
  
  
  column_pos <-
    column_info %>%
    mutate(
      do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
      xsep = case_when(
        overlay ~ c(0, -head(width, -1)),
        do_spacing ~ col_bigspace,
        TRUE ~ col_space
      ),
      xwidth = case_when(
        overlay & width < 0 ~ width - xsep,
        overlay ~ -xsep,
        TRUE ~ width
      ),
      xmax = cumsum(xwidth + xsep),
      xmin = xmax - xwidth,
      x = xmin + xwidth / 2
    )
  
  
  ##########################
  #### CREATE GEOM DATA ####
  ##########################
  
  # gather circle data
  ind_circle <- which(column_info$geom == "circle")
  dat_mat <- as.matrix(data[, ind_circle])
  col_palette <- data.frame(metric = colnames(dat_mat), 
                            group = column_info[match(colnames(dat_mat), column_info$id), "group"])
  
  col_palette$palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])
  
  circle_data <- data.frame(label = unlist(lapply(colnames(dat_mat), 
                                                  function(x) rep(x, nrow(dat_mat)))), 
                            x0 = unlist(lapply(column_pos$x[ind_circle], 
                                               function(x) rep(x, nrow(dat_mat)))), 
                            y0 = rep(row_pos$y, ncol(dat_mat)),
                            r = row_height/2*as.vector(sqrt(dat_mat))
                            #r = row_height/2*as.vector(dat_mat)
  )
  colors <- NULL
  for(i in 1:ncol(dat_mat)){
    palette <- col_palette$palette[[i]]
    colors <- c(colors, palette[rank(dat_mat[,i])])
  }
  
  circle_data$colors <- colors
  
  
  
  # gather bar data
  ind_bar <- which(column_info$geom == "bar")
  dat_mat <- as.matrix(data[, ind_bar])
  
  col_palette <- data.frame(metric = colnames(dat_mat), 
                            group = column_info[match(colnames(dat_mat), column_info$id), "group"])
  
  col_palette$palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])
  
  
  rect_data <- data.frame(label = unlist(lapply(colnames(dat_mat), 
                                                function(x) rep(x, nrow(dat_mat)))),
                          method = rep(row_info$id, ncol(dat_mat)),
                          value = as.vector(dat_mat),
                          xmin = unlist(lapply(column_pos[ind_bar, "xmin"], 
                                               function(x) rep(x, nrow(dat_mat)))),
                          xmax = unlist(lapply(column_pos[ind_bar, "xmax"], 
                                               function(x) rep(x, nrow(dat_mat)))),
                          ymin = rep(row_pos$ymin, ncol(dat_mat)),
                          ymax = rep(row_pos$ymax, ncol(dat_mat)),
                          xwidth = unlist(lapply(column_pos[ind_bar, "xwidth"], 
                                                 function(x) rep(x, nrow(dat_mat))))
  )
  rect_data <- rect_data %>%
    add_column_if_missing(hjust = 0) %>%
    mutate(
      xmin = xmin + (1 - value) * xwidth * hjust,
      xmax = xmax - (1 - value) * xwidth * (1 - hjust)
    )
  
  colors <- NULL
  for(i in 1:ncol(dat_mat)){
    palette <- col_palette$palette[[i]]
    colors <- c(colors, palette[rank(dat_mat[,i])])
  }
  
  rect_data$colors <- colors
  
  
  
  # gather text data
  ind_text <- which(column_info$geom == "text")
  dat_mat <- as.matrix(data[, ind_text])
  
  text_data <- data.frame(label_value = as.vector(dat_mat), 
                          xmin = unlist(lapply(column_pos[ind_text, "xmin"], 
                                               function(x) rep(x, nrow(dat_mat)))),
                          xmax = unlist(lapply(column_pos[ind_text, "xmax"], 
                                               function(x) rep(x, nrow(dat_mat)))),
                          ymin = rep(row_pos$ymin, ncol(dat_mat)),
                          ymax = rep(row_pos$ymax, ncol(dat_mat)),
                          size = 4)
  
  text_data$colors <- "black"
  text_data[text_data$label_value == "HVG", "colors"] <- "darkgreen"
  text_data[text_data$label_value == "FULL", "colors"] <- "grey30"
  
  text_data[text_data$label_value == "scaled", "colors"] <- "darkred"
  text_data[text_data$label_value == "unscaled", "colors"] <- "grey30"
  
  text_data[(text_data$label_value == "KNN" || text_data$label_value == "features" || 
               text_data$label_value == "embeddings"), "size"] <- 2
  
  
  # ADD COLUMN NAMES
  df <- column_pos %>% filter(id != "Method")
 
  if (nrow(df) > 0) {
    segment_data <- segment_data %>% bind_rows(
      df %>% transmute(x = x, xend = x, y = -.3, yend = -.1, size = .5)
    )
    text_data <-
      bind_rows(
        text_data,
        df %>% transmute(
          xmin = x, xmax = x, ymin = 0, ymax = -0.5,
          angle = 30, vjust = 0, hjust = 0,
          label_value = id, 
          size = 3
        )
      )
  }
 
  
  # GENERATE ROW ANNOTATION
  if (plot_row_annotation) {
    row_annotation <-
      row_pos %>% 
      select(group, ymin, ymax) %>%
      group_by(group) %>%
      summarise(
        ymin = min(ymin),
        ymax = max(ymax),
        y = (ymin + ymax) / 2
      ) %>%
      ungroup() %>%
      mutate(xmin = -.5, xmax = 5) %>%
      filter(!is.na(group), group != "")
    
    text_data <- text_data %>% bind_rows(
      row_annotation %>%
        transmute(xmin, xmax, ymin = ymax + row_space, label_value = group %>% gsub("\n", " ", .), 
                  hjust = 0, vjust = .5, fontface = "bold", size = 4) %>%
        mutate(ymax = ymin + row_height)
    )
  }
  
  
  
  # determine dsize of current geoms
  suppressWarnings({
    minimum_x <- min(column_pos$xmin, segment_data$x, segment_data$xend, circle_data$x - circle_data$r, 
                     text_data$xmin, na.rm = TRUE)
    maximum_x <- max(column_pos$xmax, segment_data$x, segment_data$xend, circle_data$x + circle_data$r, 
                     text_data$xmax, na.rm = TRUE)
    minimum_y <- min(row_pos$ymin, segment_data$y, segment_data$yend, circle_data$y - circle_data$r,  
                     text_data$ymin, na.rm = TRUE)
    maximum_y <- max(row_pos$ymax, segment_data$y, segment_data$yend,circle_data$y + circle_data$r,  
                     text_data$ymax, na.rm = TRUE)
  })
  
  ####################################
  ###   CREATE HARDCODED LEGENDS   ###
  ####################################
  
  middle_pos_x <- (maximum_x + minimum_x) /2
  leg_max_y <- minimum_y - .5
  
  
  # CREATE LEGEND for ranking colors
  leg_min_x <- middle_pos_x - 5
  rank_groups <- as.character(column_info[column_info$geom == "bar", "group"])
  
  if(usability){
    rank_minimum_x <- list("Score overall" = leg_min_x, 
                           "Removal of batch effects" = leg_min_x+1, 
                           "Cell type label variance" = leg_min_x+2,
                           "Usability" = leg_min_x+3)
  }else{
  rank_minimum_x <- list("Score overall" = leg_min_x, 
                         "Removal of batch effects" = leg_min_x+1, 
                         "Cell type label variance" = leg_min_x+2)
  }
  
  rank_title_data <- data.frame(xmin = leg_min_x, 
                                xmax = leg_min_x+ 2, 
                                ymin = leg_max_y - 1, 
                                ymax = leg_max_y, 
                                label_value = "Ranking", 
                                hjust = 0, vjust = 0, 
                                fontface = "bold")
  
  for(rg in rank_groups){
    rank_palette <- palettes[[rg]][seq(1,length(palettes[[rg]]), length.out = 5)]
    
    rank_data <- data.frame(xmin = rank_minimum_x[[rg]],
                            xmax = rank_minimum_x[[rg]] + .8,
                            ymin = seq(leg_max_y-4, leg_max_y - 2, by = .5),
                            ymax = seq(leg_max_y-3.5, leg_max_y -1.5, by = .5),
                            border = FALSE,
                            colors = rank_palette
    )
    rect_data <- bind_rows(rect_data, rank_data)
    
  }
  
  # create arrow for ranking
    arrow_data <- data.frame(x = leg_min_x +3.5, 
                           xend = leg_min_x +3.5, 
                           y = leg_max_y-4, 
                           yend = leg_max_y -1.5)
  
  # add text next to the arrow
  arrow_text <- data.frame(xmin = leg_min_x +4, 
                           xmax = leg_min_x +4.5, 
                           ymin = c(leg_max_y-2, leg_max_y-4), 
                           ymax = c(leg_max_y-1.5, leg_max_y-3.5 ), 
                           label_value = c("#1", paste0("#", nrow(data))), 
                           hjust = 0, vjust = 0, size = 2.5)
  
  
  text_data <- bind_rows(text_data, rank_title_data, arrow_text)
  
  # CREATE LEGEND for circle scores
  # circle legend
  cir_minimum_x <- middle_pos_x +1
  
  cir_legend_size <- 1
  cir_legend_space <- .1

  cir_legend_dat <-
    data.frame(
      value = seq(0, 1, by = .2),
      r = row_height/2*seq(0, 1, by = .2)
    )


  x0 <- vector("integer", nrow(cir_legend_dat))
  for(i in 1:length(x0)){
    if(i == 1){
      x0[i] <- cir_minimum_x + cir_legend_space + cir_legend_dat$r[i]
    }
    else {
      x0[i] <- x0[i-1] + cir_legend_dat$r[i-1] + cir_legend_space + cir_legend_dat$r[i]
    }
  }

  cir_legend_dat$x0 <- x0
  cir_legend_min_y <- leg_max_y-4
  cir_legend_dat$y0 <- cir_legend_min_y + 1 + cir_legend_dat$r

  cir_legend_dat$colors <- colorRampPalette(rev(brewer.pal(9, "Greys")))(nrow(cir_legend_dat))

  cir_maximum_x <- max(cir_legend_dat$x0)

  cir_title_data <- data_frame(xmin = cir_minimum_x, 
                               xmax = cir_maximum_x, 
                               ymin = leg_max_y -1, 
                               ymax = leg_max_y,
                               label_value = "Score", 
                               hjust = 0, vjust = 0, fontface = "bold")
  
  cir_value_data <- data.frame(xmin = cir_legend_dat$x0 - cir_legend_dat$r,
                               xmax = cir_legend_dat$x0 + cir_legend_dat$r,
                               ymin = cir_legend_min_y,
                               ymax = cir_legend_min_y +3,
                               hjust = .5, vjust = 0, size = 2.5,
                               label_value = ifelse(cir_legend_dat$value %in% c(0, 1), 
                                                    paste0(cir_legend_dat$value*100, "%"), ""))
  
  circle_data <- bind_rows(circle_data, cir_legend_dat)
  text_data <- bind_rows(text_data, cir_title_data, cir_value_data)


  minimum_y <- min(minimum_y, min(text_data$ymin, na.rm = TRUE))
  
  
  ########################
  ##### COMPOSE PLOT #####
  ########################
  
  g <-
    ggplot() +
    coord_equal(expand = FALSE) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_size_identity() +
    scale_linetype_identity() +
    cowplot::theme_nothing()
  
  # PLOT ROW BACKGROUNDS
  df <- row_pos %>% filter(colour_background)
  if (nrow(df) > 0) {
    g <- g + geom_rect(aes(xmin = min(column_pos$xmin)-.25, xmax = max(column_pos$xmax)+.25, ymin = ymin - (row_space / 2), ymax = ymax + (row_space / 2)), df, fill = "#DDDDDD")
  } 
  
  
  # PLOT CIRCLES
  if (nrow(circle_data) > 0) {
    g <- g + ggforce::geom_circle(aes(x0 = x0, y0 = y0, fill= colors, r = r), circle_data, size=.25)
  }
  
  
  # PLOT RECTANGLES
  if (nrow(rect_data) > 0) {
    # add defaults for optional values
    rect_data <- rect_data %>%
      add_column_if_missing(alpha = 1, border = TRUE, border_colour = "black") %>%
      mutate(border_colour = ifelse(border, border_colour, NA))
    
    g <- g + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = colors, colour = border_colour, alpha = alpha), rect_data, size = .25)
  }
  
  
  # PLOT TEXT
  if (nrow(text_data) > 0) {
    # add defaults for optional values
    text_data <- text_data %>%
      add_column_if_missing(
        hjust = .5,
        vjust = .5,
        size = 3,
        fontface = "plain",
        colors = "black",
        lineheight = 1,
        angle = 0
      ) %>%
      mutate(
        angle2 = angle / 360 * 2 * pi,
        cosa = cos(angle2) %>% round(2),
        sina = sin(angle2) %>% round(2),
        alphax = ifelse(cosa < 0, 1 - hjust, hjust) * abs(cosa) + ifelse(sina > 0, 1 - vjust, vjust) * abs(sina),
        alphay = ifelse(sina < 0, 1 - hjust, hjust) * abs(sina) + ifelse(cosa < 0, 1 - vjust, vjust) * abs(cosa),
        x = (1 - alphax) * xmin + alphax * xmax,
        y = (1 - alphay) * ymin + alphay * ymax
      ) %>%
      filter(label_value != "")
    
    g <- g + geom_text(aes(x = x, y = y, label = label_value, colour = colors, hjust = hjust, vjust = vjust, size = size, fontface = fontface, angle = angle), data = text_data)
  }
  
  
  
  # PLOT SEGMENTS
  if (nrow(segment_data) > 0) {
    # add defaults for optional values
    segment_data <- segment_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")
    
    g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), segment_data)
  }
  
  # PLOT ARROW RANKING
  if (nrow(arrow_data) > 0) {
    # add defaults for optional values
    arrow_data <- arrow_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")
    
    g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), arrow_data, arrow = arrow(length = unit(0.1, "cm")), lineend = "round", linejoin = "bevel")
  }
  
  # ADD SIZE
  # reserve a bit more room for text that wants to go outside the frame
  minimum_x <- minimum_x - 2
  maximum_x <- maximum_x + 3
  minimum_y <- minimum_y - 2
  maximum_y <- maximum_y + 3
  
  g$width <- maximum_x - minimum_x
  g$height <- maximum_y - minimum_y
  
  g <- g + expand_limits(x = c(minimum_x, maximum_x), y = c(minimum_y, maximum_y))
  return(g)
  
}



add_column_if_missing <- function(df, ...) {
  column_values <- list(...)
  for (column_name in names(column_values)) {
    default_val <- rep(column_values[[column_name]], nrow(df))
    
    if (column_name %in% colnames(df)) {
      df[[column_name]] <- ifelse(is.na(df[[column_name]]), default_val, df[[column_name]])
    } else {
      df[[column_name]] <- default_val
    }
  }
  df
}