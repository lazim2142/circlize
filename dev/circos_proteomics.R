require(circlize)

circos.trackBPlot = function(factors, y, track.height = circos.par("track.height"), track.index = NULL, force.ylim = TRUE, spacing = 0,
							 y.labels.show = FALSE, y.labels.text = NULL, y.labels.height = 0, y.labels.font = NULL, y.labels.cex = NULL, 
							 x.labels.show = FALSE, x.labels = NULL, x.labels.height = 0, x.labels.font = NULL, x.labels.cex = NULL, x.labels.col = NULL,  
                             col_pal = NULL, border = "black", lty = par("lty"), lwd = par("lwd"),
                             bg.col = "white", bg.border = "black", bg.lty = par("lty"), bg.lwd = par("lwd")) 
{
	# Length Check
	y = as.data.frame(y)
	num_samples = nrow(y)
	if (num_samples != length(factors)) {
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
	
	# Get the number of variables to iterate through
	num_vars = ncol(y)
	
	# Set a default color palette
	if(is.null(col_pal)){col_pal = rainbow(num_vars)}
		
	# Get the minimum and maximum y values across all variables
	ylim = c(min(0, apply(y, 2, min)), 
			 max(0, apply(y, 2, max)))
	
	circos.par(points.overflow.warning = FALSE)
	
	# Plot the track outline
	circos.trackPlotRegion(factors = factors, y = y[,1], track.height = track.height, track.index = track.index, 
						   force.ylim = force.ylim, ylim = ylim,
						   bg.col = bg.col, bg.border = bg.border, bg.lty = bg.lty, bg.lwd = bg.lwd,
						   panel.fun = function(x, y) {
								# X-Axis
						   		if(x.labels.show && is.null(x.labels)){
									circos.axis(major.at = c(1:length(factors)))
						   		}
						   		
						   		#Y-Axis
						   		if(y.labels.show){ 
						   			circos.yaxis(labels.cex = y.labels.cex,
						   						 labels.font = y.labels.font
						   						 )
						   		}
						   		
						   		#Y-Labels
						   		cell.xlim = get.cell.meta.data("cell.xlim")
						   		ycenter = get.cell.meta.data("ycenter")
						   		circos.text(cell.xlim[1]-y.labels.height, ycenter, paste0(y.labels.text), adj = c(0.5, 1),
						   					cex = 0.5, facing = "clockwise") 
						   })
	
	sectors = levels(factors)
	track.index = get.current.track.index()
	frequencies = table(factors)
	
	df = data.frame(factor = factors, y, x.label = x.labels)
	
	for (sector in sectors)
	{
		samples_in_sector = as.vector(frequencies[sector])
		
		xrange = get.cell.meta.data("xrange", sector.index = sector, track.index = track.index)
		sample_width = xrange/samples_in_sector
		bar_width = (sample_width-2*spacing)/num_vars
		xlim = get.cell.meta.data("xlim", sector.index = sector, track.index = track.index)
		
		for(sample in 1:samples_in_sector)
		{
			# Bars
			for(var in 1:num_vars){
				y_val = df[df$factor==sector,][sample, 1 + var]
				label = df[df$factor==sector,][sample, "x.label"]
				
				circos.rect(xlim[1] + (sample-1)*sample_width + (var-1)*bar_width + spacing, 0,
							xlim[1] + (sample-1)*sample_width + var*bar_width + spacing, y_val,
							sector.index = sector, track.index = track.index,
							col = col_pal[var], 
							border = border,
							lty = lty, 
							lwd = lwd)
			}
			
			# Lines
			if(sample > 1){
				circos.lines(c(xlim[1] + (sample-1)*sample_width, xlim[1] + (sample-1)*sample_width),
							 ylim,
							 sector.index = sector, track.index = track.index,
							 col = "Black")
			}
			
			# X-Labels
			if(x.labels.show && !is.null(x.labels))	{
				circos.text(x = xlim[1] + (sample-0.5)*sample_width, y = x.labels.height, facing = "clockwise", adj = c(0, NA),
							labels = label, col = x.labels.col, niceFacing = TRUE, font = x.labels.font, cex = x.labels.cex,
							sector.index = sector, track.index = track.index)	
			}
		}
	}
	
	# Y-Label Text
	
	
	circos.par(points.overflow.warning = TRUE)
	return(invisible(NULL))
}

													################################################
						##################################################################################################
############################################################### Data Frame Construction #################################################################
set.seed(1)
n = 50

a = data.frame(factor = rep(c("a","b","c","d","e"), length.out = n))
a$factor = sort(a$factor)

a$x = rep(0, length.out = nrow((a)))
for(level in levels(a$factor)){
	a[a$factor == level,"x"] = c(1:nrow(a[a$factor == level,]))
}

a$y1 = rnorm(n, mean = 0, sd = 1)
a$y2 = rnorm(n, mean = 1, sd = 1)
a$y3 = rnorm(n, mean = -1, sd = 1)

a$x.labels = rep(c("P1", "P2", "P3", "P4"), length.out = n)
############################################################### End Data Frame Construction ###############################################################






													################################################
						##################################################################################################
####################################################################### Load Real Data ####################################################################
if(FALSE){
nsclc_sclc = read.table(file = "~/circlize/dev/abpp_cells_qc_imputed_log2_pas_processing_12-31-15.txt", sep = '\t', header = TRUE)
adeno_sclc_m = read.table(file = "~/circlize/dev/abpp_pdx_qc_imputed_log2_pas_processing_03-01-16.txt", sep = '\t', header = TRUE)
adeno_normal_sclc_h = read.table(file = "~/circlize/dev/abpp_tissue_qc_imputed_log2_pas_processing_09-25-15.txt", sep = '\t', header = TRUE)

# Convert Rows with Scientific Notation to numeric
nsclc_sclc[,grep("e..", nsclc_sclc)] = sapply(nsclc_sclc[,grep("e..", nsclc_sclc)], as.numeric)
adeno_sclc_m[,grep("e..", adeno_sclc_m)] = sapply(adeno_sclc_m[,grep("e..", adeno_sclc_m)], as.numeric)
adeno_normal_sclc_h[,grep("e..", adeno_normal_sclc_h)] = sapply(adeno_normal_sclc_h[,grep("e..", adeno_normal_sclc_h)], as.numeric)

# Subset the data
nsclc_sclc = nsclc_sclc[abs(nsclc_sclc$log2fc_SCLC_minus_NSCLC) > 1 & abs(nsclc_sclc$wilcox.ranksum.fdr) < 0.05,]
adeno_sclc_m = adeno_sclc_m[abs(adeno_sclc_m$M_SCLC_minus_NSCLC_log2fc) > 1,]
adeno_normal_sclc_h = adeno_normal_sclc_h[abs(adeno_normal_sclc_h$H_SCLC_minus_NSCLC_log2fc) > 1,]

# Retrieve proteins common to all subsets
common <- merge(adeno_sclc_m,adeno_normal_sclc_h, row.names = c('Leading.Proteins', 'Positions'))
common <- merge(common,nsclc_sclc, row.names = c('Leading.Proteins', 'Positions'))

# Add a factors column
common$factor = factor(sample(LETTERS[1:5], nrow(common), replace = TRUE))

# Add a count of observations for each factor
for(level in levels(common$factor)){
	common[common$factor == level, "x"] = c(1:nrow(common[common$factor == level,]))
}

# Average the intensities for data set 1
sclc = common[,grep("^SCLC", colnames(common))]
common$SCLC_avg = rowMeans(sclc)

nsclc = common[,grep("^NSCLC", colnames(common))]
common$NSCLC_avg = rowMeans(nsclc)

# Average the intensities for data set 2
m_sclc = common[,grep("^M_SCLC", colnames(common))]
common$M_SCLC_avg = rowMeans(m_sclc)

m_adeno = common[,grep("^M_ADENO", colnames(common))]
common$M_ADENO_avg = rowMeans(m_adeno)

# Average the intensities for data set 3
h_sclc = common[,grep("^H_SCLC", colnames(common))]
common$H_SCLC_avg = rowMeans(h_sclc)

h_adeno = common[,grep("^H_ADENO", colnames(common))]
common$H_ADENO_avg = rowMeans(h_adeno)

h_norm = common[,grep("^H_NORM", colnames(common))]
common$H_NORM_avg = rowMeans(h_norm)
}
################################################################# Finished Loading Real Data ##############################################################




													################################################
							##################################################################################################
###################################################################### Circos Plots ########################################################################

par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("gap.degree" = 8, track.margin = c(0,0))
circos.initialize(factors = common$factor, x = common$x)

circos.trackPlotRegion(factors = common$factor, y = common$M_ADENO_1, track.height = 0.05,
					   force.ylim = TRUE, ylim = c(0,1), track.margin = c(0.1,0),
					   bg.col = c("purple","blue", "green", "red", "orange"), bg.border = "black", bg.lwd = 0.5,
					   panel.fun = function(x, y) {
					   	xcenter = get.cell.meta.data("xcenter")
					   	ycenter = get.cell.meta.data("ycenter")
					   	sector.index = get.cell.meta.data("sector.index")
					   	circos.text(xcenter, ycenter, paste0(sector.index), adj = c(0.5, 0.5),
					   				col = "white", cex = 1, niceFacing = TRUE) 
					   })

circos.trackBPlot(common$factor, common[,c("SCLC_avg","NSCLC_avg")], track.height = 0.15, spacing = 0.05,
				  x.labels.show = TRUE, x.labels = common$Leading.Proteins, x.labels.height = 30, x.labels.cex = 0.3,
				  y.labels.show = TRUE, y.labels.cex = 0.5, y.labels.text = "SCLC_avg\nNSCLC_avg", y.labels.height = 4.3,
				  col_pal = c("#084c3b", "#14c899", "#39d2f8", "#6addf9", "#045c71"), border = "white", lwd = 0.1, lty = 1,
				  bg.col = "#eefdf9", bg.border = "white", bg.lwd = 3, bg.lty = 1)

circos.trackBPlot(common$factor, common[,c("H_SCLC_avg","H_ADENO_avg","H_NORM_avg")], track.height = 0.15, spacing = 0.1,
				  x.labels.show = FALSE, x.labels = common$Leading.Proteins, x.labels.height = 30, x.labels.cex = 0.3,
				  y.labels.show = TRUE, y.labels.cex = 0.5, y.labels.text = "H_SCLC_avg\nH_ADENO_avg\nH_NORM_avg", y.labels.height = 6,
				  col_pal = c("#107382","#169aaf","#35cee5","#80e0ef"), border = "white", lwd = 0.1, lty = 1,
				  bg.col = "#eafafc", bg.border = "white", bg.lwd = 3, bg.lty = 1)

circos.trackBPlot(common$factor, common[,c("M_SCLC_avg","M_ADENO_avg")], track.height = 0.15, spacing = 0.1,
				  x.labels.show = FALSE, x.labels = common$Leading.Proteins, x.labels.height = 30, x.labels.cex = 0.3,
				  y.labels.show = TRUE, y.labels.cex = 0.5, y.labels.text = "M_SCLC_avg\nM_ADENO_avg", y.labels.height = 5.5,
				  col_pal = c("#c66e04","#fb9924","#fdcf99","fee7cb"), border = "white", lwd = 0.1, lty = 1,
				  bg.col = "#fff6ec", bg.border = "white", bg.lwd = 3, bg.lty = 1)

circos.trackBPlot(common$factor, common[,c("wilcox.ranksum.fdr")], track.height = 0.15, spacing = 0.1,
				  x.labels.show = FALSE, x.labels = common$Leading.Proteins, x.labels.height = 30, x.labels.cex = 0.3,
				  y.labels.show = TRUE, y.labels.cex = 0.5, y.labels.text = "wilcox.ranksum.fdr", y.labels.height = 7,
				  col_pal = c("magenta4","magenta3","magenta2","magenta1"), border = "white", lwd = 0.1, lty = 1,
				  bg.col = "#fef0f7", bg.border = "white", bg.lwd = 3, bg.lty = 1)




# Highlights
highlight.sector(sector.index = 'C', col = "#FF000040")

circos.clear()


if(FALSE){
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("gap.degree" = 5, track.margin = c(0,0))
circos.initialize(factors = a$factor, x = a$x)

circos.trackBPlot(a$factor, a[,c("y1","y2","y3","y2","y1")], track.height = 0.15, spacing = 0.05,
				  x.labels.show = TRUE, x.labels = a$x.labels, x.labels.height = 4, x.labels.cex = 0.7,
				  y.labels.show = TRUE, y.labels.cex = 0.5,
				  col_pal = c("#6addf9", "#045c71", "#079ec3", "#09c6f4", "#39d2f8"), border = "white", lwd = 0.1, lty = 1,
				  bg.col = "#edfbfe", bg.border = "white", bg.lwd = 3, bg.lty = 1)

circos.trackBPlot(a$factor, a[,c("y2", "y3")], track.height = 0.15, spacing = 0.1,
				  x.labels.show = FALSE, x.labels = a$x.labels, x.labels.height = 3, x.labels.cex = 0.9,
				  y.labels.show = TRUE, y.labels.cex = 0.5,
				  col_pal = c("magenta4","magenta3","magenta2","magenta1"), border = "pink", lwd = 0.1, lty = 1,
				  bg.col = "#f8e5f8", bg.border = "white", bg.lwd = 3, bg.lty = 1)

circos.trackBPlot(a$factor, a[,c("y3", "y1", "y2")], track.height = 0.15, spacing = 0.1,
				  x.labels.show = FALSE, x.labels = a$x.labels, x.labels.height = 3, x.labels.cex = 0.9,
				  y.labels.show = TRUE, y.labels.cex = 0.5,
				  col = c("#81fe5d", "green3", "green4"), border = "green4", lwd = 0.1, lty = 1,
				  bg.col = "#f5fde2", bg.border = "white", bg.lwd = 3, bg.lty = 1)
}
##################################################################### End Circos Plots #######################################################################




