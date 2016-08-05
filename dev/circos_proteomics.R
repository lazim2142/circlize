require(circlize)

circos.trackBPlot = function(factors, y, track.height = circos.par("track.height"), track.index = NULL, force.ylim = TRUE, spacing = 0,
							 y.labels.show = FALSE, y.labels.font = NULL, y.labels.cex = NULL,
							 x.labels.show = FALSE, x.labels = NULL, x.labels.height = 0, x.labels.col = NULL, x.labels.font = NULL, x.labels.cex = NULL,
                             pos.col = "blue", pos.border = "black", pos.lty = par("lty"), pos.lwd = par("lwd"),
							 neg.col = "red", neg.border = "black", neg.lty = par("lty"), neg.lwd = par("lwd"),
                             bg.col = "white", bg.border = "black", bg.lty = par("lty"), bg.lwd = par("lwd")) 
{
	# basic check here
	if (length(y) != length(factors)) {
		stop("Length of data and length of factors differ.\n")
	}
	
	if (!is.factor(factors)) {
		factors = factor(factors)
	}
	
	# check whether there are some categories that are not in the circle
	setdiff.factors = setdiff(levels(factors), get.all.sector.index())
	if (length(setdiff.factors)) {
		stop(paste("Cannot find these categories in existed sectors:", paste(setdiff.factors, collapse = ", "), ".\n", sep = ""))
	}
	
	ylim = c(min(0, y), max(0, y))
	
	circos.trackPlotRegion(factors = factors, y = y, track.height = track.height, track.index = track.index, 
						   force.ylim = force.ylim, ylim = ylim,
						   bg.col = bg.col, bg.border = bg.border, bg.lty = bg.lty, bg.lwd = bg.lwd,
						   panel.fun = function(x, y) {
								if(x.labels.show && is.null(x.labels)){
									circos.axis(major.at = c(1:length(factors)))
								}
						   		if(y.labels.show){ 
						   			circos.yaxis(labels.cex = y.labels.cex,
						   						 labels.font = y.labels.font
						   						 )
						   			 
					   			}
						   })
	
	sectors = levels(factors)
	track.index = get.current.track.index()
	frequencies = table(factors)
	
	df = data.frame(factor = factors, y = y, x.labels = x.labels)
	
	circos.par(points.overflow.warning = FALSE)
	for (sector in sectors)
	{
		num_bars = as.vector(frequencies[sector])
		xrange = get.cell.meta.data("xrange", sector.index = sector, track.index = track.index)
		bar_width = xrange/num_bars
		xlim = get.cell.meta.data("xlim", sector.index = sector, track.index = track.index)
		
		for(bar in 1:num_bars)
		{
			y_val = df[df$factor==sector,][bar, 2]
			label = df[df$factor==sector,][bar, 3]
			circos.rect(xlim[1] + (bar-1)*bar_width + spacing, 0,
						xlim[1] + bar*bar_width - spacing, y_val,
						sector.index = sector, track.index = track.index,
						col = ifelse(y_val > 0, pos.col, neg.col), 
						border = ifelse(y_val > 0, pos.border, neg.border),
						lty = ifelse(y_val > 0, pos.lty, neg.lty), 
						lwd = ifelse(y_val > 0, pos.lwd, neg.lwd))
			
			if(x.labels.show && !is.null(x.labels))
			{
				circos.text(x = xlim[1] + (bar-0.5)*bar_width, y = x.labels.height, facing = "clockwise",
							labels = label, col = x.labels.col, niceFacing = TRUE, font = x.labels.font, cex = x.labels.cex,
							sector.index = sector, track.index = track.index)	
			}
		}
	}
	circos.par(points.overflow.warning = TRUE)
	
	return(invisible(NULL))
}

set.seed(1)

n = 70

# Data Frame Construction
a = data.frame(factor = rep(c("a","b","c","d","e"), length.out = n))
a$factor = sort(a$factor)

a$x = rep(0, length.out = nrow((a)))
for(level in levels(a$factor))
{
	a[a$factor == level,2] = c(1:nrow(a[a$factor == level,]))
}

a$y = rnorm(n, mean = 0, sd = 1)
a$x.labels = rep(c("P1", "P2", "P3", "P4"), length.out = n)
# End Data Frame Construction

par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("gap.degree" = 5, track.margin = c(0,0))
circos.initialize(factors = a$factor, x = a$x)

circos.trackBPlot(a$factor, a$y, a$xlabels, track.height = 0.15, spacing = 0.1,
				  x.labels.show = TRUE, x.labels = a$x.labels, x.labels.height = 4, x.labels.cex = 0.7,
				  y.labels.show = TRUE, y.labels.cex = 0.5,
				  pos.col = "green4", pos.border = "green4", pos.lwd = 0.1, pos.lty = 1,
				  neg.col = "green3", neg.border = "green4", neg.lwd = 0.1, neg.lty = 1,
                  bg.col = "white", bg.border = "Black", bg.lwd = 0.1, bg.lty = 1)

circos.trackBPlot(a$factor, a$y, a$xlabels, track.height = 0.15, spacing = 0.1,
				  x.labels.show = FALSE, x.labels = a$x.labels, x.labels.height = 3, x.labels.cex = 0.9,
				  y.labels.show = TRUE, y.labels.cex = 0.5,
				  pos.col = "magenta4", pos.border = "green4", pos.lwd = 0.1, pos.lty = 1,
				  neg.col = "purple4", neg.border = "green4", neg.lwd = 0.1, neg.lty = 1,
				  bg.col = "white", bg.border = "Black", bg.lwd = 0.1, bg.lty = 1)

circos.trackBPlot(a$factor, a$y, a$xlabels, track.height = 0.15, spacing = 0.1,
				  x.labels.show = FALSE, x.labels = a$x.labels, x.labels.height = 3, x.labels.cex = 0.9,
				  y.labels.show = TRUE, y.labels.cex = 0.5,
				  pos.col = "red4", pos.border = "green4", pos.lwd = 0.1, pos.lty = 1,
				  neg.col = "red3", neg.border = "green4", neg.lwd = 0.1, neg.lty = 1,
				  bg.col = "white", bg.border = "Black", bg.lwd = 0.1, bg.lty = 1)

circos.clear()





