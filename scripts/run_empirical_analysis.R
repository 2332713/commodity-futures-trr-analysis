#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("empirical_R_final/scripts/run_empirical_analysis.R", winslash = "/", mustWork = TRUE)
}

script_dir <- dirname(script_path)
root_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
data_path <- file.path(root_dir, "data", "analysis_dataset.csv")
output_dir <- file.path(root_dir, "output")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

COMMODITY_ORDER <- c(
  "Copper", "Aluminum", "Rebar",
  "Crude oil", "Fuel oil", "PTA", "Methanol",
  "Soybean meal", "Corn"
)

SECTOR_BY_COMMODITY <- c(
  "Copper" = "Metals",
  "Aluminum" = "Metals",
  "Rebar" = "Metals",
  "Crude oil" = "Energy/Chemicals",
  "Fuel oil" = "Energy/Chemicals",
  "PTA" = "Energy/Chemicals",
  "Methanol" = "Energy/Chemicals",
  "Soybean meal" = "Agriculture",
  "Corn" = "Agriculture"
)

SECTOR_COLORS <- c(
  "Metals" = "#1F1F1F",
  "Energy/Chemicals" = "#5A5A5A",
  "Agriculture" = "#8C8C8C",
  "Other" = "#3A3A3A"
)

SECTOR_FILLS <- c(
  "Metals" = "#D9D9D9",
  "Energy/Chemicals" = "#BDBDBD",
  "Agriculture" = "#EFEFEF",
  "Other" = "#CCCCCC"
)

BOOTSTRAP_REPLICATIONS <- 5000L

format_num <- function(x, digits = 3) {
  if (!is.finite(x)) "--" else sprintf(paste0("%.", digits, "f"), x)
}

format_pvalue <- function(x) {
  if (!is.finite(x)) {
    "--"
  } else if (x < 0.001) {
    "<0.001"
  } else {
    sprintf("%.3f", x)
  }
}

latex_escape <- function(x) {
  out <- as.character(x)
  out <- gsub("\\\\", "\\\\textbackslash{}", out)
  out <- gsub("&", "\\\\&", out, fixed = TRUE)
  out <- gsub("%", "\\\\%", out, fixed = TRUE)
  out <- gsub("_", "\\\\_", out, fixed = TRUE)
  out
}

read_analysis_dataset <- function(path) {
  if (!file.exists(path)) {
    stop("Analysis dataset not found: ", path)
  }
  data <- read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  required <- c(
    "date", "commodity", "sector", "return", "tone", "tr", "trr", "trr_ma66",
    "next_1_extreme", "next_5_extreme", "warning"
  )
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Analysis dataset is missing required columns: ", paste(missing, collapse = ", "))
  }
  data$date <- as.Date(data$date)
  data$commodity <- factor(data$commodity, levels = COMMODITY_ORDER, ordered = TRUE)
  data <- data[order(data$commodity, data$date), ]
  data$commodity <- as.character(data$commodity)
  data$next_1_extreme <- as.logical(data$next_1_extreme)
  data$next_5_extreme <- as.logical(data$next_5_extreme)
  data$warning <- as.logical(data$warning)
  rownames(data) <- NULL
  data
}

make_ranking_table <- function(data) {
  rows <- lapply(COMMODITY_ORDER, function(commodity) {
    g <- data[data$commodity == commodity, ]
    data.frame(
      commodity = commodity,
      sector = unname(SECTOR_BY_COMMODITY[[commodity]]),
      mean_tone = mean(g$tone, na.rm = TRUE),
      mean_tr = mean(g$tr, na.rm = TRUE),
      mean_trr = mean(g$trr, na.rm = TRUE),
      median_trr = stats::median(g$trr, na.rm = TRUE),
      max_trr = max(g$trr, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  table$rank <- rank(-table$mean_trr, ties.method = "first")
  table <- table[order(table$rank), ]
  rownames(table) <- NULL
  table
}

make_warning_table <- function(data) {
  bootstrap_warning_pvalue <- function(warning, next_5_extreme, repetitions = BOOTSTRAP_REPLICATIONS) {
    n <- length(warning)
    if (n == 0 || sum(warning, na.rm = TRUE) == 0) {
      return(NA_real_)
    }

    ratios <- rep(NA_real_, repetitions)
    for (b in seq_len(repetitions)) {
      idx <- sample.int(n, size = n, replace = TRUE)
      warning_b <- warning[idx]
      next_5_b <- next_5_extreme[idx]
      if (!any(warning_b, na.rm = TRUE)) next

      unconditional_b <- mean(next_5_b, na.rm = TRUE)
      if (!is.finite(unconditional_b) || unconditional_b <= 0) next

      hit_5_b <- mean(next_5_b[warning_b], na.rm = TRUE)
      ratios[[b]] <- hit_5_b / unconditional_b
    }

    ratios <- ratios[is.finite(ratios)]
    if (length(ratios) == 0) {
      return(NA_real_)
    }
    (sum(ratios <= 1) + 1) / (length(ratios) + 1)
  }

  set.seed(20260523)
  rows <- lapply(COMMODITY_ORDER, function(commodity) {
    g <- data[data$commodity == commodity, ]
    usable <- !is.na(g$warning) & !is.na(g$next_1_extreme) & !is.na(g$next_5_extreme)
    warning_count <- sum(g$warning[usable], na.rm = TRUE)
    warning_days <- if (any(usable)) 100 * mean(g$warning[usable], na.rm = TRUE) else NA_real_
    hit_1 <- if (warning_count > 0) 100 * mean(g$next_1_extreme[usable & g$warning], na.rm = TRUE) else NA_real_
    hit_5 <- if (warning_count > 0) 100 * mean(g$next_5_extreme[usable & g$warning], na.rm = TRUE) else NA_real_
    unconditional <- if (any(usable)) 100 * mean(g$next_5_extreme[usable], na.rm = TRUE) else NA_real_
    ratio <- if (is.finite(unconditional) && unconditional > 0) hit_5 / unconditional else NA_real_
    bootstrap_p <- bootstrap_warning_pvalue(g$warning[usable], g$next_5_extreme[usable])

    data.frame(
      commodity = commodity,
      warning_days = warning_days,
      hit_1 = hit_1,
      hit_5 = hit_5,
      unconditional_hit = unconditional,
      warning_ratio = ratio,
      bootstrap_p = bootstrap_p,
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  table
}

write_ranking_latex <- function(table, path) {
  lines <- c(
    "\\begin{table}[!t]",
    "\t\\centering",
    "\t\\caption{Tail-risk resilience ranking across Chinese commodity futures}",
    "\t\\label{tab:trr_ranking}",
    "\t\\resizebox{\\textwidth}{!}{%",
    "\t\t\\begin{tabular}{llrrrrrr}",
    "\t\t\t\\toprule",
    "\t\t\tCommodity & Sector & Mean tone & Mean $TR$ & Mean $TRR$ & Median $TRR$ & Max $TRR$ & Rank \\\\",
    "\t\t\t\\midrule"
  )
  for (i in seq_len(nrow(table))) {
    row <- table[i, ]
    lines <- c(lines, paste0(
      "\t\t\t",
      paste(c(
        latex_escape(row$commodity),
        latex_escape(row$sector),
        format_num(row$mean_tone, 3),
        format_num(row$mean_tr, 3),
        format_num(row$mean_trr, 3),
        format_num(row$median_trr, 3),
        format_num(row$max_trr, 3),
        as.integer(row$rank)
      ), collapse = " & "),
      " \\\\"
    ))
  }
  lines <- c(
    lines,
    "\t\t\t\\bottomrule",
    "\t\t\\end{tabular}",
    "\t}",
    "\t\\vspace{0.3em}",
    "\t\\begin{minipage}{0.95\\textwidth}",
    "\t\t\\footnotesize",
    "\t\tNotes: $TR$ denotes downside tail risk, measured as the negative of the 5\\% rolling return quantile. $TRR$ denotes the tail-risk resilience index. A larger $TRR$ indicates weaker resilience. Rank is sorted from weakest to strongest resilience according to the mean $TRR$.",
    "\t\\end{minipage}",
    "\\end{table}",
    ""
  )
  writeLines(lines, con = path, useBytes = TRUE)
}

write_warning_latex <- function(table, path) {
  lines <- c(
    "\\begin{table}[!t]",
    "\t\\centering",
    "\t\\caption{Downside-risk warning performance}",
    "\t\\label{tab:warning}",
    "\t\\resizebox{\\textwidth}{!}{%",
    "\t\t\\begin{tabular}{lrrrrrr}",
    "\t\t\t\\toprule",
    "\t\t\tCommodity & Warning days (\\%) & 1-day hit rate (\\%) & 5-day hit rate (\\%) & Unconditional hit rate (\\%) & Warning ratio & Bootstrap $p$-value \\\\",
    "\t\t\t\\midrule"
  )
  for (i in seq_len(nrow(table))) {
    row <- table[i, ]
    lines <- c(lines, paste0(
      "\t\t\t",
      paste(c(
        latex_escape(row$commodity),
        format_num(row$warning_days, 1),
        format_num(row$hit_1, 1),
        format_num(row$hit_5, 1),
        format_num(row$unconditional_hit, 1),
        format_num(row$warning_ratio, 2),
        format_pvalue(row$bootstrap_p)
      ), collapse = " & "),
      " \\\\"
    ))
  }
  lines <- c(
    lines,
    "\t\t\t\\bottomrule",
    "\t\t\\end{tabular}",
    "\t}",
    "\t\\vspace{0.3em}",
    "\t\\begin{minipage}{0.95\\textwidth}",
    "\t\t\\footnotesize",
    "\t\tNotes: Warning days are days with high pessimistic report tone and weak tail-risk resilience. The 1-day and 5-day hit rates report the frequency of subsequent extreme downside returns. The unconditional hit rate and warning ratio use the same 5-day horizon. Bootstrap $p$-values are one-sided probabilities for the 5-day warning ratio being less than or equal to one, computed from 5,000 resamples.",
    "\t\\end{minipage}",
    "\\end{table}",
    ""
  )
  writeLines(lines, con = path, useBytes = TRUE)
}

plot_trr_time <- function(data, path_pdf, path_png, path_eps) {
  draw <- function() {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(
      bg = "white",
      mfrow = c(3, 3),
      mar = c(2.7, 3.3, 2.0, 0.7),
      oma = c(3.1, 3.3, 0.4, 0.4),
      mgp = c(2.0, 0.7, 0),
      tcl = -0.25,
      las = 1,
      cex.axis = 0.78,
      cex.main = 0.88,
      font.main = 2,
      fg = "#1F1F1F",
      col.axis = "#1F1F1F",
      col.lab = "#1F1F1F",
      col.main = "#1F1F1F",
      lend = "round"
    )
    year_ticks <- as.Date(c("2020-01-01", "2022-01-01", "2024-01-01", "2026-01-01"))
    for (commodity in COMMODITY_ORDER) {
      g <- data[data$commodity == commodity, ]
      sector <- unname(SECTOR_BY_COMMODITY[[commodity]])
      color <- unname(SECTOR_COLORS[[sector]])
      ylim <- range(g$trr_ma66, na.rm = TRUE)
      pad <- diff(ylim) * 0.08
      if (!is.finite(pad) || pad == 0) pad <- 0.005
      plot(
        g$date, g$trr_ma66,
        type = "l",
        col = color,
        lwd = 1.15,
        main = commodity,
        xlab = "",
        ylab = "",
        xaxt = "n",
        ylim = c(ylim[[1]] - pad, ylim[[2]] + pad),
        bty = "l"
      )
      axis.Date(1, at = year_ticks, format = "%Y", lwd = 0.55, lwd.ticks = 0.55)
      box(bty = "l", lwd = 0.65, col = "#1F1F1F")
      lines(g$date, g$trr_ma66, col = color, lwd = 1.15)
    }
    mtext("Trading date", side = 1, outer = TRUE, line = 1.5, cex = 0.98)
    mtext("TRR", side = 2, outer = TRUE, line = 0.2, cex = 0.95)
  }

  grDevices::pdf(path_pdf, width = 7.2, height = 6.2, family = "serif", useDingbats = FALSE)
  draw()
  grDevices::dev.off()
  grDevices::postscript(
    path_eps, width = 7.2, height = 6.2, family = "serif",
    onefile = FALSE, horizontal = FALSE, paper = "special"
  )
  draw()
  grDevices::dev.off()
  grDevices::png(path_png, width = 4320, height = 3720, res = 600, type = "cairo")
  draw()
  grDevices::dev.off()
}

plot_trr_distribution <- function(data, path_pdf, path_png, path_eps) {
  draw <- function() {
    ranking <- make_ranking_table(data)
    order_desc <- ranking$commodity
    plot_order <- rev(order_desc)
    xlim_all <- range(data$trr, na.rm = TRUE)
    xpad <- diff(xlim_all) * 0.04
    xlim_all <- c(xlim_all[[1]] - xpad, xlim_all[[2]] + xpad)
    data_list <- lapply(plot_order, function(commodity) {
      data$trr[data$commodity == commodity & is.finite(data$trr)]
    })
    box_fills <- vapply(plot_order, function(commodity) {
      sector <- unname(SECTOR_BY_COMMODITY[[commodity]])
      unname(SECTOR_FILLS[[sector]])
    }, character(1))

    add_panel_title <- function(label, x_ndc) {
      old_xpd <- par("xpd")
      par(xpd = NA)
      text(
        graphics::grconvertX(x_ndc, from = "ndc", to = "user"),
        graphics::grconvertY(par("fig")[[4]] - 0.018, from = "ndc", to = "user"),
        label,
        adj = c(0, 0.5),
        font = 2,
        cex = 0.75,
        col = "#1F1F1F"
      )
      par(xpd = old_xpd)
    }

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(1.18, 1.00), heights = c(1, 1))
    par(
      bg = "white",
      mgp = c(2.2, 0.7, 0),
      tcl = -0.25,
      las = 1,
      cex.axis = 0.74,
      fg = "#1F1F1F",
      col.axis = "#1F1F1F",
      col.lab = "#1F1F1F"
    )

    par(mar = c(3.0, 5.4, 2.0, 0.8))
    plot(
      NA,
      xlim = xlim_all,
      ylim = c(0.5, length(data_list) + 0.5),
      xlab = "TRR",
      ylab = "",
      axes = FALSE
    )
    for (i in seq_along(data_list)) {
      stats <- boxplot.stats(data_list[[i]])$stats
      segments(stats[[1]], i, stats[[2]], i, col = "#2B2B2B", lwd = 0.75)
      segments(stats[[4]], i, stats[[5]], i, col = "#2B2B2B", lwd = 0.75)
      segments(stats[[1]], i - 0.17, stats[[1]], i + 0.17, col = "#2B2B2B", lwd = 0.75)
      segments(stats[[5]], i - 0.17, stats[[5]], i + 0.17, col = "#2B2B2B", lwd = 0.75)
      rect(stats[[2]], i - 0.27, stats[[4]], i + 0.27, col = box_fills[[i]], border = "#2B2B2B", lwd = 0.75)
      segments(stats[[3]], i - 0.27, stats[[3]], i + 0.27, col = "#0F0F0F", lwd = 1.05)
    }
    axis(1, lwd = 0.55, lwd.ticks = 0.55)
    axis(2, at = seq_along(plot_order), labels = plot_order, las = 1, cex.axis = 0.62, lwd = 0.55, lwd.ticks = 0.55)
    box(bty = "l", lwd = 0.65, col = "#1F1F1F")
    add_panel_title("(a) Commodity distributions", 0.070)

    par(mar = c(3.0, 4.8, 2.0, 0.8))
    ranking_asc <- ranking[order(ranking$mean_trr), ]
    y <- seq_len(nrow(ranking_asc))
    plot(
      ranking_asc$mean_trr, y,
      type = "n",
      xlim = xlim_all,
      ylim = c(0.5, length(y) + 0.5),
      xlab = "Mean TRR",
      ylab = "",
      axes = FALSE
    )
    segments(xlim_all[[1]], y, ranking_asc$mean_trr, y, col = "#B8B8B8", lwd = 0.85)
    points(ranking_asc$mean_trr, y, pch = 21, bg = "#4A4A4A", col = "#1F1F1F", cex = 0.9, lwd = 0.7)
    axis(1, lwd = 0.55, lwd.ticks = 0.55)
    axis(2, at = y, labels = ranking_asc$commodity, las = 1, cex.axis = 0.58, lwd = 0.55, lwd.ticks = 0.55)
    box(bty = "l", lwd = 0.65, col = "#1F1F1F")
    add_panel_title("(b) Mean ranking", 0.565)

    par(mar = c(3.0, 5.4, 2.0, 0.8))
    sector_order <- c("Energy/Chemicals", "Metals", "Agriculture")
    density_list <- lapply(sector_order, function(sector) {
      stats::density(data$trr[data$sector == sector & is.finite(data$trr)], na.rm = TRUE)
    })
    ymax <- max(vapply(density_list, function(d) max(d$y), numeric(1)))
    plot(
      NA,
      xlim = xlim_all,
      ylim = c(0, ymax * 1.08),
      xlab = "TRR",
      ylab = "Density",
      axes = FALSE
    )
    density_cols <- c("#2B2B2B", "#5A5A5A", "#8C8C8C")
    density_lty <- c(1, 2, 3)
    for (i in seq_along(density_list)) {
      lines(density_list[[i]]$x, density_list[[i]]$y, col = density_cols[[i]], lwd = 1.05, lty = density_lty[[i]])
    }
    axis(1, lwd = 0.55, lwd.ticks = 0.55)
    axis(2, lwd = 0.55, lwd.ticks = 0.55)
    box(bty = "l", lwd = 0.65, col = "#1F1F1F")
    legend(
      "topright",
      legend = c("Energy/Chemicals", "Metals", "Agriculture"),
      col = density_cols,
      lty = density_lty,
      lwd = 1.05,
      bty = "n",
      cex = 0.58,
      y.intersp = 0.9,
      seg.len = 1.5
    )
    add_panel_title("(c) Sector densities", 0.070)

    par(mar = c(3.0, 4.8, 2.0, 0.8))
    sector_labels <- c("Energy/\nChem.", "Metals", "Agric.")
    sector_data <- lapply(sector_order, function(sector) {
      data$trr[data$sector == sector & is.finite(data$trr)]
    })
    sector_fills <- unname(SECTOR_FILLS[sector_order])
    boxplot(
      sector_data,
      names = sector_labels,
      col = sector_fills,
      border = "#2B2B2B",
      outline = FALSE,
      xlab = "",
      ylab = "TRR",
      frame.plot = FALSE,
      whisklty = 1,
      whisklwd = 0.75,
      staplewex = 0.6,
      boxwex = 0.52,
      medcol = "#0F0F0F",
      medlwd = 1.05,
      boxlwd = 0.75,
      staplelwd = 0.75,
      outlwd = 0.65
    )
    box(bty = "l", lwd = 0.65, col = "#1F1F1F")
    add_panel_title("(d) Sector distributions", 0.565)
  }

  grDevices::pdf(path_pdf, width = 7.2, height = 6.2, family = "serif", useDingbats = FALSE)
  draw()
  grDevices::dev.off()
  grDevices::postscript(
    path_eps, width = 7.2, height = 6.2, family = "serif",
    onefile = FALSE, horizontal = FALSE, paper = "special"
  )
  draw()
  grDevices::dev.off()
  grDevices::png(path_png, width = 4320, height = 3720, res = 600, type = "cairo")
  draw()
  grDevices::dev.off()
}

analysis_data <- read_analysis_dataset(data_path)
ranking_table <- make_ranking_table(analysis_data)
warning_table <- make_warning_table(analysis_data)

utils::write.csv(ranking_table, file.path(table_dir, "tab_trr_ranking.csv"), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(warning_table, file.path(table_dir, "tab_warning.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write_ranking_latex(ranking_table, file.path(table_dir, "tab_trr_ranking.tex"))
write_warning_latex(warning_table, file.path(table_dir, "tab_warning.tex"))
plot_trr_time(
  analysis_data,
  file.path(figure_dir, "fig_trr_time.pdf"),
  file.path(figure_dir, "fig_trr_time.png"),
  file.path(figure_dir, "fig_trr_time.eps")
)
plot_trr_distribution(
  analysis_data,
  file.path(figure_dir, "fig_trr_dist.pdf"),
  file.path(figure_dir, "fig_trr_dist.png"),
  file.path(figure_dir, "fig_trr_dist.eps")
)

cat("Done.\n")
cat("Analysis dataset:       ", normalizePath(data_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Analysis sample:        ", format(min(analysis_data$date)), " to ", format(max(analysis_data$date)), "\n", sep = "")
cat("Analysis rows:          ", nrow(analysis_data), "\n", sep = "")
cat("Weakest resilience:     ", ranking_table$commodity[[1]], "\n", sep = "")
cat("Strongest resilience:   ", ranking_table$commodity[[nrow(ranking_table)]], "\n", sep = "")
