#!/usr/bin/env Rscript
###############################################################################################

# 11/09/2020: Added linear genome visualization.

#Reading in network
suppressMessages(library(igraph))
library(ggplot2)
library(RColorBrewer)
#Load in improved plotting function
source("/home/sam/Documents/Endosymbiosis/MainModel/Xtra_Prokaryotes/scripts/igraphplot2.R")
environment(plot.igraph2) <- asNamespace('igraph')
environment(igraph.Arrows2) <- asNamespace('igraph')
options(scipen=999)

setwd("/home/sam/Documents/Endosymbiosis/MainModel/Projects/")
args <- commandArgs(trailingOnly=T)
readproject=args[1]
readnetwork=args[2]
readannotations=args[3]
readnettype=args[4]	#Can be gene-level network or complete genome with network.
if (length(args) < 4)	print("Error: not enough arguments provided.")
if (args[4] != "genes" && args[4] != "genome" && args[4] != "linear")	print("Error: network type not understood.")

#Tunable parameters
small_circ_radius <- 20	#Radius of the principal circle of the five main TFs.
large_circ <- 50	#Both radius and nr_points of the outer circle
periodicity <- 8
duplicate_space_scale <- 3	#This is used to scale the radius for each additional copy of a gene; for one of the five main TFs, you do (1 + 1/dup_space_scale)*radius. So smaller values give more spacing between paralogs.

links <- read.table(paste(readproject, readnetwork, sep="/"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
nodes <- read.table(paste(readproject, readannotations, sep="/"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (readnettype == "genes") {
	nodes <- rbind(nodes, c( paste("G",5+large_circ+1,".1",sep=""), -3.5), c( paste("G",5+large_circ+2,".1",sep=""), 3.5))	#Add two dummie genes
} else if (readnettype == "linear") {
	nodes <- rbind(nodes, c( paste("H",length(nodes[,1])+1,sep=""), "House", length(nodes[,1])+1, 0.0), c( paste("H",length(nodes[,1])+2,sep=""), "House", length(nodes[,1])+2, 0.0))	#Add two dummie genes
}

###############################################################################################
#Building the network
if (readnettype == "genes") {
	nans <- lapply(list(links$Effect), function(x) ifelse(x==0.0,NaN,x))	#We have to convert zero-interactions to "raw_effect=NaN"s, so that they are excluded just as NaN's are with the "genome" option.
	nans <- unlist(nans, use.names=FALSE)
	Edges <- data.frame(source=links$Source.gene, effect=links$Effect, target=links$Target.gene, raw_effect=nans, stringsAsFactors = FALSE)
	Vertices <- data.frame(name=nodes$Gene, threshold=as.numeric(nodes$Threshold), stringsAsFactors = FALSE)
} else {	#readnettype == "genome"
	Edges <- data.frame(source=links$Source.gene, effect=links$Net.effect, target=links$Target.gene, binding_strength=links$Binding.strength, raw_effect=links$Multiplied.activity, stringsAsFactors = FALSE)
	Vertices <- data.frame(name=nodes$Bead, type=nodes$Type, position=as.numeric(nodes$Genome.position), threshold=as.numeric(nodes$Threshold), stringsAsFactors = FALSE)
}

Network <- graph_from_data_frame(d=Edges[Edges$raw_effect!="NaN",c(1,3)], vertices=Vertices, directed=TRUE)

#Color nodes based on threshold of genes.
threshold_colors <- colorRamp(c("yellow","grey","black"))
minmax_thresh = 4.5

normalised_threshold <- (pmin(pmax(Vertices$threshold,-minmax_thresh), minmax_thresh) + minmax_thresh)/(2*minmax_thresh)
threshold <- threshold_colors( normalised_threshold )
for(i in 1:length(Vertices$name))
{
	if (readnettype == "genes") {
		threshold[i] <- rgb(threshold[i,1], threshold[i,2], threshold[i,3], maxColorValue=255)
	} else {
		if (Vertices$type[i] == "Gene") {
			threshold[i] <- rgb(threshold[i,1], threshold[i,2], threshold[i,3], maxColorValue=255)
		} else if (Vertices$type[i] == "House") {
			threshold[i] <- "#000000"
		} else {
			threshold[i] <- "#FFFFFF"
		}
	}
}

#Color edges based on their effect and set their width according to their effect (readnettype == "genes") or according to the binding probability (readnettype == "genome").
#NOTE: c(4) is binding strength if readnettype == "genome", c(5) is raw effect.
if (readnettype == "genes") {
	interactions_color_map <- colorRamp(c("blue","red"))

	normalised_effects <- (sign(Edges[Edges$raw_effect!="NaN",c(2)]) + 1) / 2
	interactions_colors <- interactions_color_map( normalised_effects )

	for(i in 1:length(Edges[Edges$raw_effect!="NaN",c(2)]))
	{
		interactions_colors[i] <- rgb(interactions_colors[i,1], interactions_colors[i,2], interactions_colors[i,3], maxColorValue=255)
	}
	interactions_widths <- 2*pmin(2*abs(Edges[Edges$raw_effect!="NaN",c(2)]), 10)
	interactions_widths_arrows <- pmax(log(interactions_widths),0.001)

} else {
	interactions_color_ramp <- colorRamp(c("blue","white","red"))
	minmax_effect = 5

	normalised_effects <- (pmin( pmax(Edges[Edges$raw_effect!="NaN",c(5)], -minmax_effect), minmax_effect) + minmax_effect)/(2*minmax_effect)
	interactions_colors <- interactions_color_ramp( normalised_effects )

	for(i in 1:length(Edges[Edges$raw_effect!="NaN",c(5)]))
	{
		interactions_colors[i] <- rgb(interactions_colors[i,1], interactions_colors[i,2], interactions_colors[i,3], maxColorValue=255)
	}
	interactions_widths <- 4*5*Edges[Edges$raw_effect!="NaN",c(4)]/2
	interactions_widths_arrows <- interactions_widths/2/2
}

###############################################################################################
#Generate table for plotting locations of genes.

arc_to_cartesion <- function(point_i, nr_points, radius, center, direction, start_arc)	#arc is in radians, so 2*pi is a full circle, pi is a half circle.
{
	if (start_arc == "top") {
		start_radian <- 0.5*pi
	} else if (start_arc == "bottom") {
		start_radian <- 1.5*pi
	}

	if (direction == "reverse") {
		dir <- -1
	} else if (direction == "forward") {
		dir <- 1
	}

	radian <- start_radian + dir*2*pi*( (point_i-1) / nr_points)	#In R tradition, we say that point_i refers to the first point on the circle.
	xcoord <- center + cos(radian)*radius
	ycoord <- center + sin(radian)*radius
	return(c(xcoord,ycoord))
}

if (readnettype == "genes") {
	GenePosition <- matrix(0, nrow=5+large_circ+2, ncol=2)	#Two dummie genes are added to make sure everything is scaled to the correct size.	0 is the center of the 100x100 field [-50,50] vs [-50,50].
	for(i in 1:length(GenePosition[,1]))
	{
		if(i <= 5) {
		  GenePosition[i,] <- arc_to_cartesion(i, 5, small_circ_radius, 0, "reverse", "top")
		} else if(i > 5 && i < length(GenePosition[,1])-1) {
		  new_index <- ((i-5)%%periodicity)*(large_circ/periodicity) + (i-5)%/%periodicity
		  GenePosition[i,] <- arc_to_cartesion(new_index, large_circ, large_circ, 0, "reverse", "top")
		} else {
		  if(i == length(GenePosition[,1])-1)	{
				GenePosition[i,] <- c(65,65)
		  } else {
				GenePosition[i,] <- c(-65,-65)
			}
		}
	}
} else if (readnettype == "genome") {
	GenePosition <- matrix(0, nrow=length(Vertices$name), ncol=2)
	for(i in 1:length(GenePosition[,1]))
	{
		GenePosition[i,] <- arc_to_cartesion(i, length(Vertices$name), length(Vertices$name), 0, "reverse", "top")
	}
} else {
	GenePosition <- matrix(0, nrow=length(Vertices$name), ncol=2)
	current_position <- 0
	for(i in 1:length(GenePosition[,1]))
	{
		if (i == length(GenePosition[,1])-1) {
			GenePosition[i,] <- c(100,-50)
		} else if (i == length(GenePosition[,1])) {
			GenePosition[i,] <- c(100,50)
		} else {
			GenePosition[i,] <- c(current_position, 0)
			if (i < length(Vertices$type)) {
				if(Vertices$type[i] == "Gene" && Vertices$type[i+1] == "Gene") {
					current_position <- current_position + 10
				} else if(Vertices$type[i] == "TFBS" && Vertices$type[i+1] == "TFBS") {
					current_position <- current_position + 10
				} else if(Vertices$type[i] == "House" && Vertices$type[i+1] == "House") {
					current_position <- current_position + 2
				} else if( (Vertices$type[i] == "Gene" && Vertices$type[i+1] == "TFBS") || (Vertices$type[i] == "TFBS" && Vertices$type[i+1] == "Gene") ) {
					current_position <- current_position + 15
				} else if( (Vertices$type[i] == "Gene" && Vertices$type[i+1] == "House") || (Vertices$type[i] == "House" && Vertices$type[i+1] == "Gene") ) {
					current_position <- current_position + 15
				} else if( (Vertices$type[i] == "TFBS" && Vertices$type[i+1] == "House") || (Vertices$type[i] == "House" && Vertices$type[i+1] == "TFBS") ) {
					current_position <- current_position + 10
				} else {
					print("No option found")
				}
			}
		}
	}
}

###############################################################################################
#Plotting the network

svg(paste(paste(readproject, readnetwork, sep="/"), ".svg", sep=""), width=15, height=15)
#png(paste(paste(readproject, readnetwork, sep="/"), ".png", sep=""), width=1000, height=1000)
par(mar=c(0,0,0,0), bg='#E0E0E0')	#light-grey background to show non-effect but site-occupying interactions!

if (readnettype == "genes") {
	custom_layout <- layout.circle(Network)
	for(i in 1:length(V(Network)))
	{
		name_pieces <- strsplit(V(Network)$name[i],"\\.")
		gene_type <- strtoi(substr(name_pieces[[1]][1],2,10))		#Substring "G15" to "15" (end is 10, so allowing enormous gene types).
		copy_number <- strtoi(name_pieces[[1]][2])
		if (gene_type <= 5) {
			scaling <- (1 + (copy_number-1)/duplicate_space_scale)
		} else {
			scaling <- (1 + (copy_number-1)/((large_circ/small_circ_radius)*duplicate_space_scale))
		}

		custom_layout[i,] <- GenePosition[gene_type,] * scaling	#For the first copy, the second term disappears; for the other copies (x,y) vector is multiplied by scalar.
	}
} else {
	custom_layout <- GenePosition
}

#Distinguish houses, tfbs and genes for the genome-level network.
Shapes <- rep("circle",length(V(Network)))
Labels <- Vertices$name
if (readnettype == "genes") {
	edge_curve <- 0.25
	Sizes <- rep(10,length(V(Network)))
	Sizes[length(V(Network))-1] <- 0
	Sizes[length(V(Network))] <- 0
} else if (readnettype == "genome") {
	edge_curve <- 0.0
	Sizes <- rep(10,length(V(Network)))
	for (i in 1:length(Sizes))
	{
		if (Vertices$type[i] == "House") {
			Sizes[i] <- 5
		} else if (Vertices$type[i] == "TFBS") {
			Sizes[i] <- 5
			Shapes[i] <- "square"
		}
	}
} else {
	edge_curve <- 0.75
	Sizes <- rep(6,length(V(Network)))
	for (i in 1:length(Sizes))
	{
		if (i == length(Sizes)-1) {
			Sizes[i] <- 0
			Labels[i] <- ""
		} else if (i == length(Sizes)) {
			Sizes[i] <- 0
			Labels[i] <- ""
		} else {
			if (Vertices$type[i] == "House") {
				Labels[i] <- ""
				Sizes[i] <- 1
			} else if (Vertices$type[i] == "TFBS") {
				Sizes[i] <- 2
				Shapes[i] <- "square"
				Labels[i] <- ""
			} else {
				Labels[i] <- substr(Labels[i],2,5)
			}
		}
	}
}

#https://igraph.org/r/doc/plot.common.html
plot.igraph2(Network,
  vertex.color = threshold,
  vertex.frame.color = "Black",
  vertex.frame.size = 15,
  vertex.shape = Shapes,
  vertex.size = Sizes,
  vertex.label = Labels,
  vertex.label.dist = 0,
  vertex.label.family = "sans",
  vertex.label.font = 1,
  vertex.label.color = "Black",
  edge.color = interactions_colors,
  edge.width = interactions_widths,
  edge.arrow.size = interactions_widths_arrows,
  edge.arrow.mode = ">",
  edge.curved = edge_curve,
  frame = TRUE,
  layout = custom_layout
  )

legend("topleft", legend=c(-minmax_thresh,-minmax_thresh/2,0,minmax_thresh/2,minmax_thresh), col = c("yellow","#C0C040","grey","404040","black"), bty = "n", pch = 20, pt.cex = 3, cex = 1.5, title="Gene threshold")
if (readnettype == "genome" || readnettype == "linear") {
	legend("bottomleft", legend=c(-minmax_effect,-minmax_effect/2,0,minmax_effect/2,minmax_effect), lty=1, lwd=10, col = c("blue","#8080FF","white","#FF8080","red"), bty = "n", pt.cex = 3, cex = 1.5, title="Interaction effect")
}

#Old threshold legend colors (yellow to blue): c("yellow","#BFBF40","#808080","#4040BF","blue")

garbage <- dev.off()
###############################################################################################
