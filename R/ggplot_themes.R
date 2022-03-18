#' ggplot2 themes
#'
#' A ggplot theme object for white background figures +/- a legend
#' @examples
#' library(ggplot2)
#' df <- data.frame(
#'   gp = factor(rep(letters[1:3], each = 10)),
#'   y = rnorm(30))
#' ggplot(df, aes(gp, y)) +
#'   geom_point() +
#'   theme_white
#' @name theme_white
#' @seealso [Plot_Coverage]
NULL

#' @describeIn theme_white White theme without figure legend
#' @export
theme_white <- theme(
    axis.line.x = element_line(colour = "black"),
    panel.grid.major = element_line(size = rel(0.5), colour = "grey"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.title.x.top = element_blank(),
    # axis.text.x.bottom = element_blank(),
    # axis.text.y = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.text = element_text(size = rel(1.0)),
    plot.title = element_text(hjust = 0.5),
    # axis.title.x=element_blank(),
    # axis.title.y=element_blank()
)

#' @describeIn theme_white White theme but with a figure legend (if applicable)
#' @export
theme_white_legend <- theme(
    axis.line.x = element_line(colour = "black"),
    panel.grid.major = element_line(size = rel(0.5), colour = "grey"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # legend.position = "none",
    axis.title.x.top = element_blank(),
    # axis.text.x.bottom = element_blank(),
    # axis.text.y = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.text = element_text(size = rel(1.0)),
    plot.title = element_text(hjust = 0.5),
    # axis.title.x=element_blank(),
    # axis.title.y=element_blank()
)

#' @describeIn theme_white White theme with figure legend but without horizontal
#' grid lines. Used internally in PlotGenome
#' @export
theme_white_legend_plot_track <- theme(
    axis.line.x = element_line(colour = "black"),
    panel.grid.major.x = element_line(size = rel(0.5), colour = "grey"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # legend.position = "none",
    axis.title.x.top = element_blank(),
    # axis.text.x.bottom = element_blank(),
    # axis.text.y = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.text = element_text(size = rel(1.0)),
    plot.title = element_text(hjust = 0.5),
    # axis.title.x=element_blank(),
    # axis.title.y=element_blank()
)
