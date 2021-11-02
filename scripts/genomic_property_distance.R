library(rtracklayer)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(regioneR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

##input: queryRegion: the query regions(GenomicRange), 
#        geno_properties: the properties need to test(List of GenomicRange), 
#        n: number of random locations to be  created (numeric)
#output: A list for each property.
genomic_proerty_distance = function(queryRegion, geno_properties, n=1e6){
  random = createRandomRegions(n, length.mean = 2, length.sd = 1, genome = "hg19", non.overlapping = T)
  res_list = list()
  i = 1
  for(geno_property in geno_properties){
    dis_random = as.data.frame(distanceToNearest(random, geno_property))
    dis_random$type = "rand"
    dis_query = as.data.frame(distanceToNearest(queryRegion, geno_property)) 
    dis_query$type = "query"
    dis = rbind(dis_query, dis_random)
    dis$rank = rank(-dis$distance, ties.method = "random")
    dis$rank_norm = (dis$rank - 1) / (max(dis$rank)-1)
    res = dis[dis$type == "query",] 
    res$med_shift = median(res$rank_norm) -0.5
    res.p = ks.test(res$rank_norm, "punif", 0, 1)$p.value
    res$pvalue = res.p
    res_list[[i]] = res
    i = i + 1
  }
  return(res_list)
}


##gpd_plot: Read the result from above function and draw the plot
gpd_plot = function(gpd_res){
  graph = ggplot(gpd_res, aes(x=rank_norm, fill = med_shift))+
    geom_density()+theme_bw() + 
   xlim(0,1)+
   xlab("")+
   scale_fill_gradient2(low = "red", mid = "white", high = "blue", limits = c(-0.3,0.3), breaks = c(-0.3,0,0.3))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  return(graph)
}
