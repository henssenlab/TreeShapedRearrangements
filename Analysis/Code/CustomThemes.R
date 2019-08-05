library(ggplot2)
theme_kons1 <- function(base_size = 12,
                        base_family = "Helvetica",
                        base_line_size = base_size / 50,
                        base_rect_size = base_size / 50) { 
  theme_classic(base_size = base_size,
                base_family = base_family, 
                base_line_size = base_line_size, 
                base_rect_size = base_rect_size) %+replace% 
    theme(
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color = "black", angle=45, vjust = 1, hjust = 1),
      axis.text = element_text(size = 10), 
      axis.title =  element_text(size = 10, face="bold"), 
      plot.title = element_text(face="plain", size = 12, hjust=0.5)
    )
}

theme_and_color_kons1 = list(
  theme_kons1,
  scale_color_brewer(palette="Set1"),
  scale_fill_brewer(palette="Set1")
)

winsorize = function (x, fraction=.05){
  lim <- quantile(x, probs=c(fraction, 1-fraction))
  x[ x < lim[1] ] <- lim[1]
  x[ x > lim[2] ] <- lim[2]
  x
}

Set1Palette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
