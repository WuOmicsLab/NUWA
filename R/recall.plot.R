
#' Plotting the recall values at both marker level and sample level.
#'
#' This function outputs a density plot and a scatter points plot to
#' visualize NUWA-ms accuracy by evaluating similarity between observed
#' and inferred abundances for markers with quantification in the inferred dataset.
#'
#' @param recallRes a list, the output of recall function.
#' @param level a character, one of "marker" (at marker level) and "sample" (at sample level). Default is "marker".
#' At marker level, scatter plot showing associations between marker level recall and correlation coefficients for all samples.
#' Dotted lines indicate a recall of 0.8, i.e. 80th percentile.
#' At sample level, density plot showing distributions of correlation coefficients within the same samples or between different samples.
#' Accuracy rate (AR), representing the overall accuracy in the dataset, and the number of comparisons is indicated.
#' 
#' @param ... additional arguments passed to the ggplot2::theme() function.
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' recall.plot(recallRes, "marker")
#' recall.plot(recallRes, "sample")
recall.plot <- function(recallRes, level = c("marker", "sample")[1], ...) {

    library(ggplot2)
    res <- recallRes

    relapos <- function(p, lev, pos) {
        if (lev == "y") {
            range <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
        } else {
            range <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
        }
        return(range[1] + pos * (range[2] - range [1]))
    }
    themeargs <- list(
        panel.background = element_blank(),
        # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5) * 0, "cm"),
        panel.grid = element_blank(),
        axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        axis.line.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 10, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 10)
        # axis.ticks.y.right = element_blank(),
    )

    themeargs <- modifyList(themeargs, list(...),keep.null = T)



    if (!level %in% c("marker", "sample")) {
        stop("Argument level should be one of 'marker' and 'sample'")
    }

    corMethod <- res$corMethod
    if (corMethod == 'pearson') {

        if (level == 'marker') {
            xlab <- 'Self Pearson r'
        } else {
            xlab <- 'Pearson r'
        }

    } else {
        if (level == 'marker') {
            xlab <- expression(paste('Self Spearman ', rho))
        } else {
            xlab <- expression(paste('Spearman ', rho))
        }

    }

    df0 <- res$recallTable
    # stop("sc")
    # Gene scatter plot
    if (level == "marker") {
        # pdf(file, width = inch(4), height = inch(4))
        df <- df0[df0$level == "Gene", ]
        df$den <- densCols(df[, "cor"], df[, "recall"])
        # abline(0.8,0,lty=2)
        bing <- length(df[, 1][df[, "recall"] >= 0.8])
        # t1 <- paste0("Well inferred\n", bing, " ", level, "s")
        badg <- length(df[, 1][df[, "recall"] < 0.8])
        # t2 <- paste0("Not well inferred\n", badg, " ", level, "s")
        rat <- round(bing / (bing + badg) * 100, digits = 1)
        t1 <- paste0("Accurate rate\n= ", rat, "%")
        p <- ggplot(df) +
            geom_point(aes(cor, recall, col = den), size = 2) +
            scale_color_identity() +
            annotate(geom = "text", x = -1, y = 0.6, label = t1, hjust = 0, size = 10 * 5/14) +
            # annotate(geom = "text", x = -.95, y = .65, label = t2, hjust = 0) +
            geom_hline(yintercept = 0.8, linetype = "dotted") +
            theme_bw() +
            do.call(theme, themeargs) +
            labs(x = xlab, y = 'Recall') +
            scale_y_continuous(breaks = c(.2, .5, .8), limits = c(0, 1)) +
            scale_x_continuous(breaks = c(-.8, 0, .8), limits = c(-1, 1))
        print(p)
        # dev.off()
    } else {
        # Sample level -------------------------------------
        df <- df0[df0$level == "Sample", ]
        simNull <- res$simNull$Sample
        simNull <- as.vector(simNull)
        simNull <- na.omit(simNull)
        # pdf(indir(d, "_sample.pdf"), width = inch(4), height = inch(4))
        t1 <- paste0("Non-self\nsamples\n(", length(simNull), ")")
        t2 <- paste0("Self samples\n(", nrow(df), ")")
        ###########################
        # cat("===============", d, "===============\n")
        # cat("median:", median(df$cor), "\n")
        # p <- wilcox.test(simNull, df$cor)$p.value
        # cat("pval:", p, "\n")
        # cat("==============================\n")
        self_median=median(df$cor, na.rm = T)
        nonself_median=median(simNull, na.rm = T)

        #############################
        p <- ggplot() +
            geom_density(aes(x = x), data = data.frame(x = simNull), fill = "gray", alpha = .5) +
            geom_density(aes(x = cor), data = df, fill = rgb(192, 192, 0, max = 255), alpha = .5) +
            do.call(theme, themeargs) +
            labs(x = xlab, y = 'Density') +
            scale_x_continuous(breaks = c(-.8, 0, .8), limits = c(-1, 1))+
            geom_vline(xintercept = self_median, lty=2)+
            geom_vline(xintercept = nonself_median, lty=2)
        y1 <- relapos(p, "y", 0.6)
        y2 <- relapos(p, "y", 0.95)
        p <- p +
            annotate(geom = "text", x = -1, y = y1, hjust = 0, vjust = 1,label = t1, size = 10 * 5 / 14,lineheight = .8) +
            annotate(geom = "text", x = -1, y = y2, hjust = 0,vjust = 1, label = t2, size = 10 * 5 / 14,lineheight = .8) +
            annotate(geom = "text", x = nonself_median, y = y2, hjust = 0,vjust = 1, label = round(nonself_median, 2), size = 10 * 5 / 14,lineheight = .8)+
            annotate(geom = "text", x = self_median, y = y2, hjust = 0,vjust = 1, label = round(self_median, 2), size = 10 * 5 / 14,lineheight = .8)
        print(p)
        # dev.off()
    }

}
