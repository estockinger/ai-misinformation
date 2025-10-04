
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(viridis)
library(circlize)
library(ggsankey)
library(ggrepel)

source("scripts/utils.R")
source("scripts/bibliometrix_wrapper.R")

colored_thematicMap <- function(df, field="KW_Merged", n=1000, minfreq=5, stemming=T, size=.1, n.labels=4, repel=T, synonyms=c(), subgraphs=T) {
  Map <- thematicMap(df, field = field, n = n, minfreq = minfreq, stemming = stemming, size = size,
                     n.labels = n.labels, repel = repel, synonyms=synonyms, subgraphs=subgraphs)
  
  # Replotting the map to match the colorscheme
  meandens=mean(Map$map$data$rdensity)
  meancentr=mean(Map$map$data$rcentrality)
  rangex=max(c(meancentr-min(Map$map$data$rcentrality),max(Map$map$data$rcentrality)-meancentr))
  rangey=max(c(meandens-min(Map$map$data$rdensity),max(Map$map$data$rdensity)-meandens))
  
  factor <- c(1.2, 1)
  xlimits=c(meancentr-(rangex*factor[1]),meancentr+rangex*factor[1])
  ylimits=c(meandens-(rangey*factor[2]),meandens+rangey*factor[2])
  
  
  annotations <- data.frame(
    xpos = c(xlimits[1], xlimits[1], xlimits[2], xlimits[2]), 
    ypos = c(ylimits[1], ylimits[2], ylimits[1], ylimits[2]),
    words = c("Emerging or Declining","Niche","Basic","Motor"),
    hjustvar = c(0,0,1,1) ,
    vjustvar = c(0,1,0,1))
  
  Map$map$data$color <- plyr::mapvalues(Map$map$data$color, unique(Map$map$data$color), viridis(length(unique(Map$map$data$color))))
  
  Map$colored_map <- ggplot(Map$map$data, aes(x=rcentrality, y=rdensity, text=c(words))) + 
    geom_hline(yintercept = meandens,linetype=2, color=adjustcolor("black",alpha.f=0.3)) +
    geom_vline(xintercept = meancentr,linetype=2, color=adjustcolor("black",alpha.f=0.3)) +
    theme(legend.position="none") +
    scale_radius(range=c(5*(1+size), 30*(1+size)))  +
    xlim(xlimits)+
    ylim(ylimits)+
    labs(x = "Relevance degree\n(Centrality)", y = "Development degree\n(Density)") +
    annotate("text", x=annotations$xpos,y= annotations$ypos,hjust=annotations$hjustvar,
             vjust=annotations$vjustvar,label=annotations$words, color=adjustcolor("gray20", alpha.f=0.5))+
    geom_point(group="NA",aes(size=log(as.numeric(freq))),shape=20, col=adjustcolor(Map$map$data$color,alpha.f=.4)) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    geom_label_repel(
      aes(group="NA",label=ifelse(freq>1,unlist(tolower(name_full)),'')),
      size=2*(1+size),
      angle=0,
      label.size = NA, 
      lineheight=1.05,
      box.padding = 0.1,
      fill = alpha(c("white"),0), 
      point.size = NA, max.overlaps=Inf, segment.colour = NA) +
    theme(axis.text.x=element_blank(),
          panel.background = element_rect(fill = '#FFFFFF'),
          axis.line.x = element_blank(), 
          axis.line.y = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), 
    )
  
  Map
}

keywords_over_time <- function(data, Tag, synonyms, show_top_n, label_top_n) {
  plotData <- topKWs(data, synonyms=synonyms, top=show_top_n, Tag=Tag) |>
    summarise(Occurrences, Prognosed = Occurrences + Extrap, .by=c(Year, Tag))
  topKwAnytime <- plotData |> summarize(total=sum(Occurrences), .by=Tag) |> slice(1:label_top_n) |> pull(Tag)
  
  p <- ggplot(data=plotData, aes(x = Year, node=Tag, fill = Tag, value = Prognosed, label=Tag)) +
    geom_sankey_bump() + 
    theme_minimal()
  
  layer_data <- layer_data(plot=p) |>
    summarize(y = mean(y), .by=c("x", "label"))
  
  g_labs_start <- layer_data |>
    filter(x == min(x) & !(label %in% topKwAnytime)) |>
    left_join(group_by(plotData[, c("Year", "Tag", "Occurrences")], Tag) |> filter(Year == min(Year)), by = c("label" = "Tag"))
    
  g_labs_end <- layer_data |>
    filter(x == max(x) & label %in% topKwAnytime) |>
    left_join(group_by(plotData[, c("Year", "Tag", "Occurrences")], Tag) |> filter(Year == max(Year)), by = c("label" = "Tag"))
  
  palette <- viridis(show_top_n)
  
  p <- ggplot(data=plotData, aes(x = Year, node=Tag, fill = Tag, value = Prognosed, next_node = Tag, next_x = Tag)) +
    geom_sankey_bump(aes(alpha = if_else(Tag %in% topKwAnytime, Tag, NA))) +
    geom_text_repel(
      data = g_labs_start, 
      aes(x, y, label = str_to_sentence(label), color = label),
      inherit.aes = F, min.segment.length = 0, size= 3, nudge_x = -.1, direction = "y", hjust = "right", force=.5, segment.size=0.2, 
    ) +
    geom_text_repel( 
      data = g_labs_end, 
      aes(x, y, label = paste(str_wrap(str_to_sentence(label), 26), "-", Occurrences), color = label),
      inherit.aes = F, min.segment.length = 0, size= 3, nudge_x = .1, direction = "y", hjust = "left", force=.5, segment.size=0.2, lineheight=.8
    ) +
    scale_x_continuous("", seq(2022,2025), seq(2022,2025), limits=c(2020.8, 2026.3)) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    scale_alpha_manual(values = rep(1, 15), na.value = .5) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
    )
  p
}

country_collaborations <- function(data, filename, tofile=T) {
  countryCountryDf <- countryCountryCollabs(data)
  countryCollabDf <- countryCollabStats(countryCountryDf, replace_smallest_n = 25)
  
  group <- countryCollabDf[, c('continent', 'chordEntry')] |> unique()
  group <- setNames(group$continent, group$chordEntry)
  grid.col <- setNames(viridis(length(group)), names(sort(group)))
  
  countryCountryDf <- countryCountryDf |> 
    left_join(countryCollabDf[, c("country", "chordEntry")], by=c("from"="country")) |> 
    left_join(countryCollabDf[, c("country", "chordEntry")], by=c("to"="country"))
  
  if (tofile) { 
    pdf(file = filename, width = 5.5, height = 5.5) 
  }
  circos.clear()
  circos.par(start.degree = 50)
  chordDiagram(
    countryCountryDf[, c("chordEntry.x", "chordEntry.y", "value")], 
    group = group,
    big.gap=5, 
    small.gap=0,
    grid.col = grid.col, 
    annotationTrack = "grid", 
    preAllocateTracks = 1)
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    index <- if_else(CELL_META$sector.index == "USA", "USA", str_wrap(str_to_title(CELL_META$sector.index), 10))
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], index, facing = "clockwise", niceFacing = TRUE, adj = c(0, .5), cex=.8)
  }, bg.border = NA) 
  if (tofile) {
    dev.off()
  }
}

topics_by_method <- function(
    data, grp, subgrp, elem, grpnames,
    factors = c(0, 1, 1.5, 2, 2.5, 3, 4) * .1,
    gridYs = c(2, 4, 6, 8, 10, 12, 14)) {
  
  plotData <- data |>
    summarise(count = n(), .by = c(!!grp, !!subgrp, !!elem)) |>
    mutate({{grp}} := factor({{grp}}, levels = grpnames)) |>
    group_by({{grp}}) |>
    group_modify(~add_row(., {{subgrp}}:="", count=0)) |>
    group_by(!!grp, !!subgrp)  |>
    mutate(id = cur_group_id()) |>
    ungroup()
  
  # Get the name and the y position of each label
  labels <- plotData |> 
    summarize(total=sum(count), .by=c(id, !!subgrp)) |> 
    mutate(total = replace_na(total, 0))
  angle <- 90 - 360 * (labels$id-0.5) / nrow(labels)
  labels$angle <- ifelse(angle < -90, angle+180, angle)
  
  # prepare a data frame for base lines
  base <- plotData |>
    summarize(start=min(id)+1, end=max(id), .by={{grp}}) |>
    rowwise() |>
    mutate(title=mean(c(start, end)),
           code = case_match({{grp}}, !!!map2(grpnames, names(grpnames), ~{.x ~.y})))
  
  # prepare a data frame for grid (scales)
  scales <- map2(gridYs, factors, ~geom_segment(
    data=data.frame(start = base$start[-1] - .8, end = base$end[1:nrow(base)-1] + .8), 
    aes(x = end -.y, y = .x, xend = start + .y, yend = .x), 
    colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE )
  )
  
  # Make the plot
  ggplot(plotData) +
    geom_bar(aes(x=as.factor(id), y=count, fill={{elem}}), stat="identity", alpha=.8) +
    scales +
    annotate("text", 
             x = rep(min(plotData$id), length(gridYs)), y = gridYs-.5, label = gridYs , 
             color="grey", size=4, angle=0, fontface="bold", hjust=.5) +
    scale_fill_viridis(discrete=TRUE, na.translate=F,  name="Study types") +
    ylim(-7, max(labels$total+2, na.rm=T)) +
    theme_minimal() + 
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    geom_text(data=labels, 
              aes(x=id, y=total + 1, label=str_wrap(theme, 23-total[1]), hjust="outward"), 
              color="black", fontface="bold",alpha=0.6, size=4, angle= labels$angle, inherit.aes = FALSE, lineheight=.8) +
    geom_segment(data=base, 
                 aes(x = start, y = -1, xend = end, yend = -1), 
                 colour = "black", alpha=0.8, size=1 , inherit.aes = FALSE ) +
    geom_text(data=base, 
              aes(x = title, y = -2, label=code, hjust="inward"), 
              colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
}

method_legend <- function(labels, filename, chunks = 1, textlendiv = 170, xjust = 0, x = 0, ymax = 1, ymin=.5) {
    slice <- ceiling(seq_along(labels)/ceiling(length(labels)/chunks))
    args <- data.frame(labels = labels, colors = viridis(length(labels)), widths = map_vec(labels, \(x) str_length(x)/textlendiv))
    
    pdf(file = filename, width = 15, height = 3) 
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    map2(split(args, slice), 
         seq(ymax, ymin, (ymin-ymax)/chunks)[1:chunks], 
         ~legend(x=x, y=.y, legend=.x$labels, fill = .x$colors, border=.x$colors, text.width=.x$widths, xjust=xjust, bty = "n", horiz=T, inset=-2))
    dev.off()
}
