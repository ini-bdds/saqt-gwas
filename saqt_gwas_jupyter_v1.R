
# Library of functions for exploring and visualizating results of SAQT-GWAS


initialize_code <- function() {
	library(devtools)
	library(plotly)
	load_all('/Users/willmatloff/Documents/LONI/R/Jupyter_GWAS/Manhattanly/manhattanly-master') #add path to Manhattanly directory here
	library(tools)
	library(plyr)
	library(dplyr)
}

initialize <- function() {
	suppressMessages(initialize_code())
}

view_phenotypes <- function(location) {
	files <- list.files(path=paste(location,"/PLINK_results_linear_annotated", sep=""), full.names=T, recursive=F)
	files <- files[which(file.info(files)$size>0)] #get files that don't have size = 0
	names <- lapply(files, function(filename) {return(file_path_sans_ext(file_path_sans_ext(basename(filename))))})
	names
}

linear_assoc_manhattan_plot <- function(location, phenotype) {
	pvaldata <- read.table(paste(location,"/PLINK_results_linear_annotated/",phenotype,".assoc.linear", sep=""), sep="\t", header=TRUE, colClasses=c("numeric", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric", "character", "character", "character"), stringsAsFactors=FALSE)
	colnames(pvaldata)[12] <- "QT"
	p <- manhattanly(pvaldata, point_size = 3, suggestiveline = -log10(1e-06), genomewideline = -log10(5e-08), snp = "SNP", gene = "GeneSymbol", annotation1 = "GeneLocation", annotation2 = "QT")
	embed_notebook(p)
}

read_table_filename <- function(filename){
	ret <- read.table(filename, sep="\t", header=TRUE, colClasses=c("numeric", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric", "character", "character", "character"), stringsAsFactors=FALSE)
	if (nrow(ret) > 0) {
		colnames(ret)[12] <- "QT"
		ret
	}
	else {
		return(NULL)
	}
}

combined_linear_assoc_manhattan_plot <- function(location) {
	files <- list.files(path=paste(location,"/PLINK_results_linear_annotated", sep=""), full.names=T, recursive=F)
	files <- files[which(file.info(files)$size>0)] #get files that don't have size = 0
	dataframe_list <- lapply(files, read_table_filename)
	combined <- do.call(rbind, dataframe_list)
	p <- manhattanly(combined, point_size = 3, suggestiveline = -log10(1e-06), genomewideline = -log10(5e-08), snp = "SNP", gene = "GeneSymbol", annotation1 = "GeneLocation", annotation2 = "QT")
	embed_notebook(p)
}

view_combined_top_snps_linear <- function(location, number) {
	files <- list.files(path=paste(location,"/PLINK_results_linear_annotated", sep=""), full.names=T, recursive=F)
	files <- files[which(file.info(files)$size>0)] #get files that don't have size = 0
	dataframe_list <- lapply(files, read_table_filename)
	combined <- do.call(rbind, dataframe_list)
	above6 <- subset(combined, P < .00001)
	sorted <- arrange(above6, P)
	#print(head(sorted,number))
	print(head(subset(sorted, select=c("SNP", "P", "GeneSymbol", "QT")), number))
}

view_top_snps_linear <- function(location, phenotype, number) {
	pvaldata <- read.table(paste(location,"/PLINK_results_linear_annotated/",phenotype,".assoc.linear", sep=""), sep="\t", header=TRUE, colClasses=c("numeric", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric", "character", "character", "character"), stringsAsFactors=FALSE)
	colnames(pvaldata)[12] <- "QT"
	sorted <- arrange(pvaldata, P)
	print(head(subset(sorted, select=c("CHR", "P", "GeneSymbol", "QT")), number))
}

logistic_assoc_manhattan_plot <- function(location) {
	pvaldata <- read.table(paste(location,"/PLINK_results_logistic_annotated/PLINK_results_logistic.assoc.logistic", sep=""), sep="\t", header=TRUE, colClasses=c("numeric", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric", "character", "character"), stringsAsFactors=FALSE)
	p <- manhattanly(pvaldata, point_size = 3, suggestiveline = -log10(1e-06), genomewideline = -log10(5e-08), snp = "SNP", gene = "GeneSymbol", annotation1 = "GeneLocation")
	embed_notebook(p)
}

view_top_snps_logistic <- function(location, number) {
	pvaldata <- read.table(paste(location,"/PLINK_results_logistic_annotated/PLINK_results_logistic.assoc.logistic", sep=""), sep="\t", header=TRUE, colClasses=c("numeric", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric", "character", "character"), stringsAsFactors=FALSE)
	sorted <- arrange(pvaldata, P)
	print(head(subset(sorted, select=c("CHR", "SNP", "P", "GeneSymbol", "GeneLocation")), number))
}

phe_by_groups <- function(location, phenotype) {
	phedata <- read.table(paste(location,"/final_phe_files/",phenotype,".phe", sep=""), sep=" ", header=FALSE, colClasses=c("numeric", "character", "numeric"), stringsAsFactors=FALSE)
	colnames(phedata) <- c("FID", "PTID", "PHE")
	final_merged <- read.csv(paste(location,"/final_merged.csv", sep=""), header=TRUE, stringsAsFactors=FALSE)
	phe_with_groups <- merge(subset(phedata, select=c("PTID","PHE")), subset(final_merged, select=c("PTID","Screen.Diagnosis")), by=c("PTID"))
	p <- plot_ly(phe_with_groups, y = ~PHE, color = ~Screen.Diagnosis, type = "box")
	embed_notebook(p)	
}


view_external_data <- function(ext_data_loc) {
	comp_data <- read.csv(ext_data_loc, header=TRUE)
	print(comp_data)
}

compare_two <- function(ext_data_loc, location, p_value_name) {
	comp_data <- read.csv(ext_data_loc, header=TRUE, stringsAsFactors=FALSE)
	#print(comp_data)
	files <- list.files(path=paste(location,"/PLINK_results_linear_annotated", sep=""), full.names=T, recursive=F)
	files <- files[which(file.info(files)$size>0)] #get files that don't have size = 0
	dataframe_list <- lapply(files, read_table_filename)
	combined <- do.call(rbind, dataframe_list)
	comp_data_filtered <- subset(comp_data, select=c("SNP", "QT"))
	joined <- left_join(comp_data_filtered, combined, by=c("SNP", "QT"))
	joined_filtered <- subset(joined, select=c("SNP", "QT", "P"))
	merged <- merge(comp_data, joined_filtered, by=c("SNP", "QT"))
	merged <- subset(merged, select=c("SNP", "GeneName", "QT", p_value_name, "P"))
	merged$P <- round(-log10(merged$P),1)
	#final_comparison <- arrange(merged, desc(P_Name))
	final_comparison <- merged[order(-merged[[p_value_name]]),]
	rownames(final_comp_test) <- NULL
	print(final_comparison)
}




