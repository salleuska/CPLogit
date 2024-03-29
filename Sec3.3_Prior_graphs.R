####################
## Sally Paganin
## Nov 4 2021 
####################
####################
library(CPLogit)
library(ggplot2)
library(grid)
library(gridExtra)
library(BNPmix)
####################
## Check CP process using normalized distance
####################
# Bell numbers - sequence A000110
bell_n = c(1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975, 678570, 4213597, 27644437, 190899322, 1382958545, 10480142147, 82864869804)

####################
## Functions for the CP prior
####################


CPpenalization <- function(distance, psi, returnLog = FALSE) {
  if(returnLog){
    -psi*distance
  } else {
    exp(-psi*distance)
  }
}

####################
## Dirichlet Process EPPF e
####################

dpEPPF <- function(partition, concentration, returnLog = FALSE){
  nObs = length(partition)
  nBlocks = length(unique(partition))
  blockSizes = table(partition)
  
  logProb = nBlocks*log(concentration) + sum(lgamma(blockSizes))
  
  if(returnLog) {
     logProb   
  } else {
     exp(logProb)  
  }
}

###########################
## Finite Dirichlet EPPF
###########################

finiteDirichletEPPF = function(partition, concentration, nClusters, returnLog = FALSE){
  nObs = length(partition)
  nBlocks = length(unique(partition))
  blockSizes = table(partition)
  
  log_c = 0
  if(nBlocks < nClusters) {log_c = - lgamma(nClusters - nBlocks +1)}
  logProb =   lgamma(nClusters +1) + log_c + sum(lgamma(concentration/nClusters + blockSizes) - lgamma(concentration/H))
  ##  + lgamma(concentration) - lgamma(concentration + nObs) # (normalizaoitn constant)
  
  if(returnLog) {
     logProb   
  } else {
       exp(logProb)  
  }
}

##############################
## Pitman-Yor process eppf
##############################
pyEPPF <- function(partition, concentration, discount, returnLog = FALSE){
  n = length(partition)
  nBlocks = length(unique(partition))
  blockSizes = table(partition)
  
  if(nBlocks == 1){
    logProb = -
      lgamma(n + concentration) + lgamma(concentration + 1) + 
      sum(lgamma(discount + blockSizes) - lgamma(1 - discount))   

  } else {
    logProb =  sum(log(concentration + 1:(nBlocks-1)*discount)) - 
      lgamma(n + concentration) + lgamma(concentration + 1) + 
      sum(lgamma(discount + blockSizes) - lgamma(1 - discount))   
  }
  
  
  if(returnLog) {
     logProb   
  } else {
     exp(logProb)  
  }
}

##################
## Plotting
##################
## Set c0 partition color
col2 = "#00B0DA"

## Set a list of c0 partitions
c0_list = list(c(0,0,0,0,0), 
            c(0,0,0,1,1),
            c(0,0,1,1,2),
            c(0,0,1,2,3),
            c(0,1,2,3,4))

##################
## UNIFORM EPPF 
##################
## values of psi parameter (x-axis)
psiVal = c(0,0.5,1,1.5,2,2.5, 3)


for(i in 1:5){
  
  c0 = c0_list[[i]]
  ## Compute distances from c0 and generate all possible partitions
  dd = CPLogit::dist_from(c0, TRUE)
  
  ## Numerator of CP process
  CPNum = sapply(psiVal, function(x) CPpenalization(dd$distances, psi = x))
  ## Compute probabilities (normalization)
  CPProbs = apply(CPNum, 2, function(x) x/sum(x))
  
  nBlocks = apply(dd$partition, 1, function(x) length(unique(x)))
  biggestBlock = apply(dd$partition, 1, function(x) max(table(x)))
  indBiggestBlock = order(as.numeric(paste(nBlocks, biggestBlock, sep = "")), decreasing = TRUE)
  # 
  
  tmp = CPProbs[indBiggestBlock, ]
  df = reshape2::melt(tmp) ## data frame Var1 = partition index Var2 psi index 

  ## Find position of the c0 partition
  c0Pos = which(apply(dd$partitions[indBiggestBlock, ], 1, function(x) all(x == c0)))

  df$isC0 = as.numeric(df$Var1==c0Pos)

  yTicks = c(0, cumsum(tmp[,1])[c(1,16, 41, 51, 52)])
  p = ggplot(df, aes(x = Var2, y = value)) + 
    geom_area(aes(group = Var1, fill = factor(isC0)),col = "grey20" ) +
    scale_fill_manual(values = c('grey90',col2))+
    theme_minimal(18) +
    scale_x_continuous(breaks = 1:7,labels = psiVal, limits = c(1, 7), expand = c(0, 0)) +
    theme(legend.position = 'none',plot.margin=grid::unit(c(1,1,1,1), "cm")) + 
    scale_y_continuous(breaks= yTicks, limits = c(0, 1), expand = c(0, 0),
                       sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                           labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
    theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
          axis.ticks.y.left = element_line(size = 0.6, colour = "grey20", linetype = "dotted"),
          axis.ticks.length.y.left  = unit(2, "cm"),
          axis.ticks.length.x.bottom  = unit(2, "cm"),
          axis.text.y.right = element_text( size = rel(1.6)),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
          axis.title.x=element_text(size=20,face="bold"),
           axis.title.y.right=element_text(size=20)
          ) +
    xlab(expression(psi)) 
  
  pt = p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.01,ymax = 0.01,  xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.2, ymax = 0.2,  xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.55,ymax = 0.55, xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.85,ymax = 0.85,  xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "5 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.99,ymax = 0.99,  xmin = 0,xmax = 0.5)


  gt <- ggplot_gtable(ggplot_build(pt))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid_p = grid.arrange(gt)
  

 ggsave(grid_p, file = paste0("CP_uniform_", i, ".pdf"), device = "pdf", width = 6, height = 5, unit = 'cm', scale = 3, dpi = 300)
  
}

###########################
## DP process ----- 
###########################
alpha <- 1

for(i in 1:5){
  
  c0 = c0_list[[i]]
  ## Compute distances from c0 and generate all possible partitions
  dd = CPLogit::dist_from(c0, TRUE)
  
  ## Numerator of CP process penalization
  CPNum = sapply(psiVal, function(x) CPpenalization(dd$distances, psi = x))
  ## Numerator of DP process EPP
  DPNum = apply(dd$partitions, 1 ,function(x) dpEPPF(partition = x, concentration =  alpha))
  
  # p_dir = apply(dd$partitions, 1 ,function(x) dir_eppf(x, 1/2, 5))
  
  ### ORDER PARTITIONS
  nBlocks = apply(dd$partition, 1, function(x) length(unique(x)))
  biggestBlock = apply(dd$partition, 1, function(x) max(table(x)))
  indBiggestBlock = order(as.numeric(paste(nBlocks, biggestBlock, sep = "")), decreasing = TRUE)
  # 
  
  partitionProbs = apply(CPNum, 2, function(x) (DPNum*x)/sum(DPNum*x))
  
  tmp = partitionProbs[indBiggestBlock, ]
  df = reshape2::melt(tmp) ## data frame Var1 = partition index Var2 psi index 

  ## Find position of the c0 partition
  c0Pos = which(apply(dd$partitions[indBiggestBlock, ], 1, function(x) all(x == c0)))

  df$isC0 = as.numeric(df$Var1==c0Pos)

  yTicks = c(0, cumsum(tmp[,1])[c(1,16, 41, 51, 52)])
  p = ggplot(df, aes(x = Var2, y = value)) + 
    geom_area(aes(group = Var1, fill = factor(isC0)),col = "grey20" ) +
    scale_fill_manual(values = c('grey90',col2))+
    theme_minimal(18) +
    scale_x_continuous(breaks = 1:7,labels = psiVal, limits = c(1, 7), expand = c(0, 0)) +
    theme(legend.position = 'none',plot.margin=grid::unit(c(1,1,1,1), "cm")) + 
    scale_y_continuous(breaks= yTicks, limits = c(0, 1), expand = c(0, 0),
                       sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                           labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
    theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
          axis.ticks.y.left = element_line(size = 0.6, colour = "grey20", linetype = "dotted"),
          axis.ticks.length.y.left  = unit(2, "cm"),
          axis.ticks.length.x.bottom  = unit(2, "cm"),
          axis.text.y.right = element_text( size = rel(1.6)),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
          axis.title.x=element_text(size=20,face="bold"),
           axis.title.y.right=element_text(size=20)
          ) +
    xlab(expression(psi)) 
  
  pt = p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.01,ymax = 0.01,  xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.2, ymax = 0.2,  xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.55,ymax = 0.55, xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.85,ymax = 0.85,  xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "5 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.99,ymax = 0.99,  xmin = 0,xmax = 0.5)


  gt <- ggplot_gtable(ggplot_build(pt))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid_p = grid.arrange(gt)
  

  ggsave(grid_p, file = paste0("CP_dp_", i, ".pdf"), device = "pdf", height = 15, width = 15, unit = "cm", dpi = 300) 
  
}


###########################
## Pitman Yor ----- 
###########################
## Fix discount parameter and pick concentration to 
## match expected number of clusters under the DP
discount = 0.5

## discount = 0.75
## Check the expected number of clusters matches DP (1.6)
## BNPmix::PYcalibrate(log(5), 5, discount = discount)$strength

for(i in 1:5){
  
  c0 = c0_list[[i]]
  ## Compute distances from c0 and generate all possible partitions
  dd = CPLogit::dist_from(c0, TRUE)
  
  ## Numerator of CP process penalization
  CPNum = sapply(psiVal, function(x) CPpenalization(dd$distances, psi = x))

  alpha = BNPmix::PYcalibrate(log(5), 5, discount = discount)$strength

  ## Numerator of PY process EPP
  PYNum = apply(dd$partitions[-1, ], 1 ,function(x) pyEPPF(partition = x, concentration =  alpha, discount = discount))
  PYNum = c(1 - sum(PYNum), PYNum)


  ### ORDER PARTITIONS
  nBlocks = apply(dd$partition, 1, function(x) length(unique(x)))
  biggestBlock = apply(dd$partition, 1, function(x) max(table(x)))
  indBiggestBlock = order(as.numeric(paste(nBlocks, biggestBlock, sep = "")), decreasing = TRUE)
  # 
  
  partitionProbs = apply(CPNum, 2, function(x) (PYNum*x)/sum(PYNum*x))
  
  tmp = partitionProbs[indBiggestBlock, ]
  df = reshape2::melt(tmp) ## data frame Var1 = partition index Var2 psi index 

  ## Find position of the c0 partition
  c0Pos = which(apply(dd$partitions[indBiggestBlock, ], 1, function(x) all(x == c0)))

  df$isC0 = as.numeric(df$Var1==c0Pos)

  yTicks = c(0, cumsum(tmp[,1])[c(1,16, 41, 51, 52)])
  p = ggplot(df, aes(x = Var2, y = value)) + 
    geom_area(aes(group = Var1, fill = factor(isC0)),col = "grey20" ) +
    scale_fill_manual(values = c('grey90',col2))+
    theme_minimal(18) +
    scale_x_continuous(breaks = 1:7,labels = psiVal, limits = c(1, 7), expand = c(0, 0)) +
    theme(legend.position = 'none',plot.margin=grid::unit(c(1,1,1,1), "cm")) + 
    scale_y_continuous(breaks= yTicks, limits = c(0, 1), expand = c(0, 0),
                       sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                           labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
    theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
          axis.ticks.y.left = element_line(size = 0.6, colour = "grey20", linetype = "dotted"),
          axis.ticks.length.y.left  = unit(2, "cm"),
          axis.ticks.length.x.bottom  = unit(2, "cm"),
          axis.text.y.right = element_text( size = rel(1.6)),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
          axis.title.x=element_text(size=20,face="bold"),
           axis.title.y.right=element_text(size=20)
          ) +
    xlab(expression(psi)) 
  
  pt = p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.01,ymax = 0.01,  xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.2, ymax = 0.2,  xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.55,ymax = 0.55, xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.85,ymax = 0.85,  xmin = 0,xmax = 0.5) +
    annotation_custom(grob = textGrob(label = "5 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.99,ymax = 0.99,  xmin = 0,xmax = 0.5)


  gt <- ggplot_gtable(ggplot_build(pt))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid_p = grid.arrange(gt)
  

  ggsave(grid_p, file = paste0("CP_py_", i, "_sigma05_ekDP.pdf"), device = "pdf", height = 15, width = 15, unit = "cm", dpi = 300) 
  
}


