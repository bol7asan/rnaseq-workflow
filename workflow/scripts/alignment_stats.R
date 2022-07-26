#!/usr/bin/env Rscript

#Load the necessary libraries
library(plotly)
library(readr)
library(htmlwidgets)
library(widgetframe)
library(DESeq2)
source("/media/asangani2/scripts/downstream-analysis/R/functions.R")

#Read feature count assembled table.
df_raw <- CleanTables(snakemake@params[[1]],format = "featurecounts")


##Remove samples with 0 variance
norm.cnt <- t(get_normalized_counts(df_raw))
norm.cnt <- norm.cnt[,which(apply(norm.cnt,2,var) != 0 )]

#Perform PCA
pca <- prcomp(norm.cnt, scale. = TRUE )
pca.data <- data.frame(Sample = rownames(pca$x), X = pca$x[,1], Y = pca$x[,2])
groups <- sub("\\.[0-9]|(_[^_]+$)","",pca.data$Sample)
pca.data["Group"] <- groups


##Configurations for the plot
noax = list(
  title = "",
  zeroline = FALSE,
  showline = FALSE
)

a <- list(
  text = "PCA",
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "center",
  x = 0.5,
  y = 1,
  showarrow = FALSE
)

b <- list(
  text = "Alignment Stats",
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "center",
  x = 0.5,
  y = 1,
  showarrow = FALSE
)

##First plot, PCA
fig <- plot_ly(data = pca.data, x = ~X,y=~Y, type = 'scatter',mode = 'markers',symbol = ~Group,text = ~Sample, marker = list(size=10), hoverinfo = 'text' ) %>% layout(annotations = a,xaxis = noax,yaxis=noax )
fig

#Glob alignment stats for each sample.
workdir <- dirname(snakemake@params[[1]])
files <- Sys.glob(paste0(workdir,"/*"))
samps <- c()
nms <- c()
stats <- c()

##Loop over files and construct the stats data frame.
for (file in files) {
  sample.name <- basename(file)
  if (sample.name != "count-tables" & sample.name != "ucsc"){
    stat.dir <- paste0(file,"/STAR/Log.final.out")
    log <- read.table(stat.dir, sep = "\t", fill = T, skip = 5)
    samps <- c(samps,sample.name)
    nms <- c(nms,as.numeric(sub("%","",log[5,"V2"])))
    stat.labels <- c("Reads","Avg Read Length","Mapped Reads (Unique)","% Mapped Reads (Unique) ","% Mapped Reads (Multimappers)","% Mapped Reads (Many Loci)", "% Unmapped Reads (Many Mismatches)", "% Unmapped Reads (Too Short)", "% Unmapped Reads (Other)")
    stat.labels <- paste0("<b>",stat.labels,"</b>")
    stat.nums <- log$V2[c(1,2,4,5,20,22,25,27,29)]
    stats.combined <- paste(stat.labels,": ",stat.nums)
    stats.combined <- paste(stats.combined, collapse = "<br>")
    stats <- c(stats,stats.combined)
  }
  
}

##Second plot: Bar plot containing alignment stats.
df.bar <- data.frame(sampls = samps, nms = nms,text=stats)

fig1 <- plot_ly(df.bar,
                x=~sampls,
                y=~nms,
                type = "bar",
                hovertext=~text,
                hoverinfo = "text",
                showlegend = FALSE) %>% layout(annotations = b,hoverlabel = list(align="left", bgcolor = "white"),yaxis = list(title = "Unique Alignment %",range = list(0,100)), xaxis = list(title = ""))




##Combine plots and save as html widget
fg <- subplot(fig,fig1, margin = 0.07) %>% layout(autosize = F, width = 1200, height = 500) 
frameableWidget(fg)
saveWidget(fg, snakemake@output[[1]], selfcontained = T)