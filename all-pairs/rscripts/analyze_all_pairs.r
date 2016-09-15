# Generates a png image of the median RTD of each pair of ranks

require('ggplot2')
require('reshape')

analyze_all_pairs <- function(raw_data, kernel, kind, outputdir){
  time_cols        = c(3:length(raw_data))
  cores_per_host   = 8 # Adjust to your needs
  cores_per_socket = 4 # Adjust to your needs
  upper_bound      = 600
  top_limit        = 12
  title            = paste("all pairs RTD in microsec (", kernel , ")", sep="")
  kernel_fname     = gsub(" ", "_", fixed = TRUE, kernel) # replace spaces
	plotfile				 = paste(kernel_fname, kind, "png", sep=".")
  plotfile         = paste(outputdir, plotfile, sep="/")
 
  measure_data <- raw_data 
  
  # combine all aggregated data
  clear_data <- melt(measure_data) 
  # name columns for plotting
  names(clear_data) <- c("sender", "reciever", "data") 

  # set outliers to top_limit
  clear_data$data[clear_data$data > top_limit] <- top_limit

  # set correct name
#  colnames(clear_data)[3] <- kind
  maj_grid_lines = seq(0, dim(raw_data)[1], cores_per_host)
  min_grid_lines = seq(cores_per_socket, dim(raw_data)[1], cores_per_host)
  pal <- colorRampPalette(c("blue", "cyan", "yellow", "red", "deeppink"))
  #plot
  p <- ggplot(clear_data, aes(x=sender-0.5, y=reciever-0.5))
  p <- p + ggtitle(title) + ylab("sender rank") + xlab("reciever rank")
  p <- p + geom_tile(aes(fill=data))
	p <- p + scale_fill_gradientn( colours = pal(100), limits = c(0,top_limit), guide = guide_legend(title = kind))
  p <- p + scale_x_continuous(minor_breaks = maj_grid_lines)
  p <- p + scale_y_continuous(minor_breaks = maj_grid_lines)

# mark host and socket boundaries
#  p <- p + geom_vline(xintercept=maj_grid_lines, alpha=1/8, size=.5) # major v lines
#  p <- p + geom_hline(yintercept=maj_grid_lines, alpha=1/8, size=.5) # major h lines

  p <- p + theme_bw()
  p <- p + coord_fixed(ratio = 1)
  p <- p + theme(plot.title=element_text(face="bold", size=20))
  ggsave(filename = plotfile, plot=p, units = "mm", width = 200, dpi = 400)
  
  return (clear_data)
}
