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
####################################################################### Load Protein Data ####################################################################
if(FALSE){ # Start file reads and subsets for mass-spec data
	nsclc_sclc = read.table(file = "~/circlize/dev/abpp_cells_qc_imputed_log2_pas_processing_12-31-15.txt", sep = '\t', header = TRUE)
	adeno_sclc_m = read.table(file = "~/circlize/dev/abpp_pdx_qc_imputed_log2_pas_processing_03-01-16.txt", sep = '\t', header = TRUE)
	adeno_normal_sclc_h = read.table(file = "~/circlize/dev/abpp_tissue_qc_imputed_log2_pas_processing_09-25-15.txt", sep = '\t', header = TRUE)

	# Convert Rows with Scientific Notation to numeric
	nsclc_sclc[,grep("e..", nsclc_sclc)] = sapply(nsclc_sclc[,grep("e..", nsclc_sclc)], as.numeric)
	adeno_sclc_m[,grep("e..", adeno_sclc_m)] = sapply(adeno_sclc_m[,grep("e..", adeno_sclc_m)], as.numeric)
	adeno_normal_sclc_h[,grep("e..", adeno_normal_sclc_h)] = sapply(adeno_normal_sclc_h[,grep("e..", adeno_normal_sclc_h)], as.numeric)

	# Subset the mass-spec data by lg2fc and q_values
	nsclc_sclc = nsclc_sclc[abs(nsclc_sclc$log2fc_SCLC_minus_NSCLC) > 1 & abs(nsclc_sclc$wilcox.ranksum.fdr) < 0.05,]
	adeno_sclc_m = adeno_sclc_m[abs(adeno_sclc_m$M_SCLC_minus_NSCLC_log2fc) > 1 & (T || abs(adeno_sclc_m$M_SCLC_minus_NSCLC_log2fc) < 0.05),]
	adeno_normal_sclc_h = adeno_normal_sclc_h[abs(adeno_normal_sclc_h$H_SCLC_minus_NSCLC_log2fc) > 1 & (T || abs(adeno_normal_sclc_h$H_NSCLC_minus_normal_log2fc) < 0.05),]
	
	# Retrieve proteins common to all subsets
	prot_common = merge(adeno_sclc_m,adeno_normal_sclc_h, row.names = c('Leading.Proteins', 'Positions'))
	prot_common = merge(prot_common,nsclc_sclc, row.names = c('Leading.Proteins', 'Positions'))
	
	# Clean Up Environment
	rm(list = c("nsclc_sclc", "adeno_sclc_m", "adeno_normal_sclc_h"))
	
	# Average the intensities for data set 1 
	prot_common$mean_SCLC = rowMeans(prot_common[,grep("^SCLC", colnames(prot_common))])
	prot_common$mean_NSCLC = rowMeans(prot_common[,grep("^NSCLC", colnames(prot_common))])
	
	# Average the intensities for data set 2
	prot_common$mean_M_SCLC = rowMeans(prot_common[,grep("^M_SCLC_[0-9]", colnames(prot_common))])
	prot_common$mean_M_ADENO = rowMeans(prot_common[,grep("^M_ADENO_[0-9]", colnames(prot_common))])
	
	# Average the intensities for data set 3
	prot_common$mean_H_SCLC = rowMeans(prot_common[,grep("^H_SCLC_[0-9]", colnames(prot_common))])
	prot_common$mean_H_ADENO = rowMeans(prot_common[,grep("^H_ADENO_[0-9]", colnames(prot_common))])
	prot_common$mean_H_NORM = rowMeans(prot_common[,grep("^H_NORM_[0-9]", colnames(prot_common))])
	
	# Remove columns after finding the means
	prot_common[,grep("^SCLC|^NSCLC|^M_SCLC_[0-9]|^M_ADENO_[0-9]|^H_SCLC_[0-9]|^H_ADENO_[0-9]|^H_NORM_[0-9]", colnames(prot_common))] = NULL
} # End file reads and subsets for mass-spec data
##################################################################### Done Loading Protein Data ####################################################################


															################################################
									##################################################################################################
####################################################################### Load Metabolite Data ####################################################################
if(FALSE){ # Start file reads and subsets for metabolite data
	metab_pos = read.table(file = "~/circlize/dev/metab_pos_impute_avg_03-14-16.txt", sep = '\t', header = TRUE)
	metab_neg = read.table(file = "~/circlize/dev/metab_neg_impute_avg_03-14-16.txt", sep = '\t', header = TRUE)
	
	# Convert Rows with Scientific Notation to numeric
	metab_pos[,grep("e..", metab_pos)] = sapply(metab_pos[,grep("e..", metab_pos)], as.numeric)
	metab_neg[,grep("e..", metab_neg)] = sapply(metab_neg[,grep("e..", metab_neg)], as.numeric)
	
	# Subset metabolites by annotation availability
	metab_pos = droplevels(metab_pos[as.character(metab_pos$Compound) != as.character(metab_pos$Accepted.Description),])
	metab_neg = droplevels(metab_neg[as.character(metab_neg$Compound) != as.character(metab_neg$Accepted.Description),])
	
	# Drop metabolites that are listed as a class/family
	for (level in levels(metab_pos$Accepted.Description)){
		if(nrow(metab_pos[metab_pos$Accepted.Description == level,]) > 1){
			metab_pos = metab_pos[metab_pos$Accepted.Description != level,]
		}
	}
	for (level in levels(metab_neg$Accepted.Description)){
		if(nrow(metab_neg[metab_neg$Accepted.Description == level,]) > 1){
			metab_neg = metab_neg[metab_neg$Accepted.Description != level,]
		}
	}
	
	# Create common data frame of metabolites, preferring those with lower pval
	metab_common = rbind(metab_pos, metab_neg)
	metab_common = merge(metab_common, aggregate(FDR ~ Accepted.Description, metab_common, min))
	
	# Clean Up Environment
	rm(list = c("metab_pos", "metab_neg"))
}
##################################################################### Done Loading Metabolite Data ####################################################################

															################################################
								##################################################################################################
############################################################# Getting mapping from UniProtKB to Gene Names ###################################################################
if(FALSE){
	# Generate keys for look up online
	uniprotkb_names = {}
	for(proteins in prot_common$Leading.Proteins){
		proteins = gsub("CON__", '', proteins)
		split_names = strsplit(proteins, ";")
		for(name in split_names)	{
			uniprotkb_names = append(uniprotkb_names, name)
		}
	}
	uniprotkb_names = unique(uniprotkb_names)
	
	# Download UniProt Data
	require(UniProt.ws)
	uniprotkb_gene_map = select(UniProt.ws(), 
						  keys = uniprotkb_names, 
						  columns = "GENES",
						  keytype = "UNIPROTKB")
	
	# Add gene names to prot_common data frame
	prot_common$Gene.Symbol = unlist(lapply(prot_common$Gene.Symbol, as.character))
	for(proteins in prot_common$Leading.Proteins){
		# Get the names of UniProtKB items to look up
		split_names = strsplit(proteins, ";")
		
		# Get the gene names for each UniProtKB item
		mappings = uniprotkb_gene_map[uniprotkb_gene_map$UNIPROTKB == split_names,"GENES"]
		
		# Keep at most the first two mappings
		mappings = unlist(strsplit(mappings, "; "))
		if(!is.null(mappings))	{
			for(i in 1:length(mappings)){
				mapping = unlist(strsplit(mappings[i], ' '))
				if (length(mapping) > 2) {
					mappings[i] = c(mapping[1:2]) 
				}
			}
		}
		
		gene_names = ""
		gene_names = paste(gene_names, mappings, sep = " ", collapse = ";")
		
		if(gene_names != "" & gene_names != " ")	{
			prot_common[prot_common$Leading.Proteins==proteins, "Gene.Symbol"] = gene_names	
		}
	}
}
########################################################## Done getting mapping from UniProtKB to Gene Names ###################################################################


															################################################
									##################################################################################################
######################################################### Get mapping from UniProtKB to associated metabolites #######################################################
# HMDB XML Parsing to get mapping from UniProtKB to associated metabolites
if(FALSE){
	require("XML")
	# Setup to parse each XML section as a different file because the file is not properly formatted
	# Start contains a vector of line numbers where an xml file begins
	# Start contains a vector of line numbers where xml file ends
	hmdb_proteins = readLines("~/circlize/dev/hmdb_proteins.xml")
	start = grep('<?xml version="1.0" encoding="UTF-8"?>', hmdb_proteins, fixed=T)
	end = c(start[-1]-1,length(hmdb_proteins))

	# Iterate through all the xml "files" to collect information
	uniprotkb_metab_map = data.frame(Leading.Proteins = character(), Metabolites = character(), stringsAsFactors = F)
	
	for(i in 1:length(start)){
		xml = xmlParse(paste(hmdb_proteins[start[i]:end[i]],collapse="\n"))
		metab_refs = paste(unique(xpathApply(xml, "//metabolite_references/metabolite_reference/metabolite/name", xmlValue)), collapse = ", ")
		if(metab_refs != ""){
			uniprotkb_metab_map[nrow(uniprotkb_metab_map) + 1, 1] = unlist(xpathApply(xml, "//uniprot_id", xmlValue))
			uniprotkb_metab_map[nrow(uniprotkb_metab_map), 2] = metab_refs
		}
	}
	
	rm(hmdb_proteins)
	
} # End HMDB XML Parsing

################################################################# Finished Parsing HMDB Data ##############################################################
if(TRUE)
{
	common = merge(prot_common, uniprotkb_metab_map, by = "Leading.Proteins", all = T)
}


													################################################
							##################################################################################################
###################################################################### Circos Plots ########################################################################
if(FALSE){
	# Add a factors column
	prot_common$factor = factor(sample(LETTERS[1:5], nrow(prot_common), replace = TRUE))
	
	# Add a count of observations for each factor because the barplot needs an "x-axis"
	for(level in levels(prot_common$factor)){
		prot_common[prot_common$factor == level, "x"] = c(1:nrow(prot_common[prot_common$factor == level,]))
	}
	
	par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
	circos.par("gap.degree" = 8, track.margin = c(0,0))
	circos.initialize(factors = prot_common$factor, x = prot_common$x)
	
	circos.trackPlotRegion(factors = prot_common$factor, track.height = 0.05,
						   force.ylim = TRUE, ylim = c(0,1), track.margin = c(0.1,0),
						   bg.col = c("purple","blue", "green", "red", "orange"), bg.border = "black", bg.lwd = 0.5,
						   panel.fun = function(x, y) {
							   	xcenter = get.cell.meta.data("xcenter")
							   	ycenter = get.cell.meta.data("ycenter")
							   	sector.index = get.cell.meta.data("sector.index")
							   	circos.text(xcenter, ycenter, paste0(sector.index), adj = c(0.5, 0.5),
							   				col = "white", cex = 1, niceFacing = TRUE) 
						   })
	
	circos.trackBPlot(prot_common$factor, prot_common[,c("mean_SCLC","mean_NSCLC")], track.height = 0.15, spacing = 0.05,
					  x.labels.show = TRUE, x.labels = prot_common$Gene.Symbol, x.labels.height = 30, x.labels.cex = 0.3,
					  y.labels.show = TRUE, y.labels.cex = 0.5, y.labels.text = "SCLC_mean\nNSCLC_mean", y.labels.height = 4.3,
					  col_pal = c("#084c3b", "#14c899", "#39d2f8", "#6addf9", "#045c71"), border = "white", lwd = 0.1, lty = 1,
					  bg.col = "#eefdf9", bg.border = "white", bg.lwd = 3, bg.lty = 1)
	
	circos.trackBPlot(prot_common$factor, prot_common[,c("mean_H_SCLC","mean_H_ADENO","mean_H_NORM")], track.height = 0.15, spacing = 0.1,
					  x.labels.show = FALSE, x.labels = prot_common$Leading.Proteins, x.labels.height = 30, x.labels.cex = 0.3,
					  y.labels.show = TRUE, y.labels.cex = 0.5, y.labels.text = "H_SCLC_mean\nH_ADENO_mean\nH_NORM_mean", y.labels.height = 6,
					  col_pal = c("#107382","#169aaf","#35cee5","#80e0ef"), border = "white", lwd = 0.1, lty = 1,
					  bg.col = "#eafafc", bg.border = "white", bg.lwd = 3, bg.lty = 1)
	
	circos.trackBPlot(prot_common$factor, prot_common[,c("mean_M_SCLC","mean_M_ADENO")], track.height = 0.15, spacing = 0.1,
					  x.labels.show = FALSE, x.labels = prot_common$Leading.Proteins, x.labels.height = 30, x.labels.cex = 0.3,
					  y.labels.show = TRUE, y.labels.cex = 0.5, y.labels.text = "M_SCLC_mean\nM_ADENO_mean", y.labels.height = 5.5,
					  col_pal = c("#c66e04","#fb9924","#fdcf99","fee7cb"), border = "white", lwd = 0.1, lty = 1,
					  bg.col = "#fff6ec", bg.border = "white", bg.lwd = 3, bg.lty = 1)
	
	circos.trackBPlot(prot_common$factor, prot_common[,c("M_SCLC_avg","M_ADENO_avg")], track.height = 0.15, spacing = 0.1,
					  x.labels.show = FALSE, x.labels = prot_common$Leading.Proteins, x.labels.height = 30, x.labels.cex = 0.3,
					  y.labels.show = TRUE, y.labels.cex = 0.5, y.labels.text = "M_SCLC_avg\nM_ADENO_avg", y.labels.height = 5.5,
					  col_pal = c("#c66e04","#fb9924","#fdcf99","fee7cb"), border = "white", lwd = 0.1, lty = 1,
					  bg.col = "#fff6ec", bg.border = "white", bg.lwd = 3, bg.lty = 1)
	
	#circos.trackBPlot(prot_common$factor, prot_common[,c("wilcox.ranksum.fdr")], track.height = 0.15, spacing = 0.1,
	#				  x.labels.show = FALSE, x.labels = prot_common$Leading.Proteins, x.labels.height = 30, x.labels.cex = 0.3,
	#				  y.labels.show = TRUE, y.labels.cex = 0.5, y.labels.text = "wilcox.ranksum.fdr", y.labels.height = 7,
	#				  col_pal = c("magenta4","magenta3","magenta2","magenta1"), border = "white", lwd = 0.1, lty = 1,
	#				  bg.col = "#fef0f7", bg.border = "white", bg.lwd = 3, bg.lty = 1)
	
	# Highlights
	highlight.sector(sector.index = 'C', col = "#FF000040")
	
	circos.clear()
} ##################################################################### End Circos Plots #######################################################################

