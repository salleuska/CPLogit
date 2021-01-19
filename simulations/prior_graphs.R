bell_n = c(1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975, 678570, 4213597, 27644437, 190899322, 1382958545, 10480142147, 82864869804)
####################à
## functions for CP prior
p_c = function(delta, omega, psi)
{
  norm = sum(omega*exp(-psi*delta))
  exp(-psi*delta)/norm
}
p_clust = function(dist, psi)
{
  exp(-psi*dist)
}

####################
## DP eppf

dp_eppf = function(part, alpha){
  n = length(part)
  n_blocks = length(unique(part))
  block_sizes = table(part)
  
  log_prob = n_blocks*log(alpha) + sum(lgamma(block_sizes))
  ##  + lgamma(alpha) - lgamma(alpha + n) # (cost)
  
  exp(log_prob)
}

## finite dirichlet eppf

dir_eppf = function(part, alpha, H){
  n = length(part)
  n_blocks = length(unique(part))
  block_sizes = table(part)
  
  log_c = 0
  if(n_blocks < H) {log_c = - lgamma(H - n_blocks +1)}
  log_prob =   lgamma(H +1) + log_c + sum(lgamma(alpha/H + block_sizes) - lgamma(alpha/H))
  ##  + lgamma(alpha) - lgamma(alpha + n) # (cost)
  
  exp(log_prob)
}

## Pitman-Yor process eppf
## CORRECTED!!
py_eppf = function(part, alpha, sigma){
  n = length(part)
  n_blocks = length(unique(part))
  block_sizes = table(part)
  
  if(n_blocks == 1){
    log_prob = -
      lgamma(n + alpha) + lgamma(alpha + 1) + 
      sum(lgamma(sigma + block_sizes) - lgamma(1 - sigma))   

  } else {
    log_prob =  sum(log(alpha + 1:(n_blocks-1)*sigma)) - 
      lgamma(n + alpha) + lgamma(alpha + 1) + 
      sum(lgamma(sigma + block_sizes) - lgamma(1 - sigma))   
  }
  
  
  exp(log_prob)
}

# sigma = 0.5
# alpha = BNPmix::PYcalibrate(4, 5, discount = sigma)$strength
# py_eppf(c0_list[[1]], alpha, sigma)
# 
# dd = sharedLogit::dist_from(c0_list[[1]], TRUE)

####################
library(sharedLogit)
library(ggplot2)
library(grid)
library(gridExtra)
library(BNPmix)
library(numbers)
##################

col2 = "lightblue"

c0_list = list(c(0,0,0,0,0), 
            c(0,0,0,1,1),
            c(0,0,1,1,2),
            c(0,0,1,2,3),
            c(0,1,2,3,4))
####################
##################
## UNIFORM -----
##################

for(i in 1:5){
  c0 = c0_list[[i]]
  dd = sharedLogit::dist_from(c0, TRUE)
  tt = table(dd$distances)
  
  ## uniform
  psi_val = c(0,0.5,1,1.5,2,2.5, 3)
  num = sapply(psi_val, function(x) p_clust(dd$distances, psi = x))
  probs = apply(num, 2, function(x) x/sum(x))
  
  n_blocks = apply(dd$partition, 1, function(x) length(unique(x)))
  biggest_block = apply(dd$partition, 1, function(x) max(table(x)))
  ind = order(as.numeric(paste(n_blocks, biggest_block, sep = "")), decreasing = TRUE)
  # n_blocks[ind]
  as.numeric(paste(n_blocks, biggest_block, sep = ""))[ind]
  
  prova = probs[ind, ]
  colnames(prova) = psi_val
  df = reshape2::melt(prova)
  
  c0_pos = 52 - which(apply(dd$partitions[ind, ], 1, function(x) all(x == c0))) + 1
  
  y_ticks = c(0, cumsum(prova[,1])[c(1,16, 41, 51, 52)])
  
  
  p = ggplot(df, aes(x = Var2, y = value)) + geom_area(aes(group = Var1, fill = Var1), col = "grey20", fill = ifelse(df$Var1 == c0_pos, col2, "grey90")) +
    theme_minimal() +
    scale_x_continuous(breaks= psi_val, limits = c(0, 3), expand = c(0, 0)) +
    theme(plot.margin=grid::unit(c(1,1,1,1), "cm")) +
    # scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0)) + 
    scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0),
                       sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                           labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
    theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
          axis.ticks.y.left = element_line(size = 0.6, colour = "grey20", linetype = "dotted"),
          axis.ticks.length.y.left  = unit(2, "cm"),
          axis.ticks.length.x.bottom  = unit(2, "cm"),
          axis.text.y.right = element_text( size = rel(1.6)),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
          axis.title.x=element_text(size=20,face="bold")) +
    xlab(expression(psi)) 
  p
  
  pt = p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.01,ymax = 0.01,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.2, ymax = 0.2,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.55,ymax = 0.55, xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.85,ymax = 0.85,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "5 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.99,ymax = 0.99,  xmin = -1,xmax = 0)
  
  gt <- ggplot_gtable(ggplot_build(pt))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid_p = grid.arrange(gt)
  
  
   ggsave(grid_p, file = paste0("CP_uniform_", i, ".pdf"), device = "pdf", height = 14, width = 16, unit = "cm", dpi = 300)
  
}

###########################
## DP process ----- 
###########################


for(i in 1:5){
  c0 = c0_list[[i]]
  dd = sharedLogit::dist_from(c0, TRUE)
  tt = table(dd$distances)

  
  psi_val = c(0,0.5,1,1.5,2,2.5, 3)
  num = sapply(psi_val, function(x) p_clust(dd$distances, psi = x))
  alpha = 1
  p_dir = apply(dd$partitions, 1 ,function(x) dp_eppf(x, alpha))
  
  # p_dir = apply(dd$partitions, 1 ,function(x) dir_eppf(x, 1/2, 5))
  
  ### ORDER PARTITIONS
  n_blocks = apply(dd$partition, 1, function(x) length(unique(x)))
  biggest_block = apply(dd$partition, 1, function(x) max(table(x)))
  ind = order(as.numeric(paste(n_blocks, biggest_block, sep = "")), decreasing = TRUE)

  probs = apply(num, 2, function(x) (p_dir*x)/sum(p_dir*x))
  
  prova = probs[ind, ]
  colnames(prova) = psi_val
  df = reshape2::melt(prova)

  c0_pos = 52 - which(apply(dd$partitions[ind, ], 1, function(x) all(x == c0))) + 1
  
  y_ticks = c(0, rev(cumsum(rev(prova[,1]))[rev(c(1,16, 41, 51, 52))]))
  
  
  p = ggplot(df, aes(x = Var2, y = value)) + geom_area(aes(group = Var1, fill = Var1), col = "grey20", fill = ifelse(df$Var1 == c0_pos, col2, "grey90")) +
    theme_minimal() +
    scale_x_continuous(breaks= psi_val, limits = c(0, 3), expand = c(0, 0)) +
    theme(plot.margin=grid::unit(c(1,1,1,1), "cm")) +
    # scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0)) + 
    scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0),
                       sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                           labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
    theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
          axis.ticks.y.left = element_line(size = 0.6, colour = "grey20", linetype = "dotted"),
          axis.ticks.length.y.left  = unit(2, "cm"),
          axis.ticks.length.x.bottom  = unit(2, "cm"),
          axis.text.y.right = element_text( size = rel(1.6)),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
          axis.title.x=element_text(size=20,face="bold")) +
    xlab(expression(psi)) 
  
  p
  
  pt = p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.03,ymax = 0.03,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.4, ymax = 0.4,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.8,ymax = 0.8, xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.92,ymax = 0.92,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "5 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.98,ymax = 0.98,  xmin = -1,xmax = 0)
  
  gt <- ggplot_gtable(ggplot_build(pt))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid_p = grid.arrange(gt)
  
  ggsave(grid_p, file = paste0("CP_dp_", i, ".pdf"), device = "pdf", height = 15, width = 15, unit = "cm", dpi = 300) 
}


###########################
## Pitman Yor ----- 
###########################
sigma = 0.5

sigma = 0.75
BNPmix::PYcalibrate(log(5), 5, discount = sigma)$strength
## Expected number of clusters matches DP (1.6)

for(i in 1:5){
  c0 = c0_list[[i]]
  
  dd = sharedLogit::dist_from(c0, TRUE)
  tt = table(dd$distances)
  
  psi_val = c(0,0.5,1,1.5,2,2.5, 3)
  num = sapply(psi_val, function(x) p_clust(dd$distances, psi = x))
  alpha = BNPmix::PYcalibrate(log(5), 5, discount = sigma)$strength
  
  p_dir = apply(dd$partitions[-1, ], 1 ,function(x) py_eppf(x, alpha, sigma))
  p_dir = c(1 - sum(p_dir), p_dir)
  
  ### ORDER PARTITIONS
  n_blocks = apply(dd$partition, 1, function(x) length(unique(x)))
  biggest_block = apply(dd$partition, 1, function(x) max(table(x)))
  ind = order(as.numeric(paste(n_blocks, biggest_block, sep = "")), decreasing = TRUE)

  probs = apply(num, 2, function(x) (p_dir*x)/sum(p_dir*x))
  
  prova = probs[ind, ]
  colnames(prova) = psi_val
  df = reshape2::melt(prova)
  
    c0_pos = 52 - which(apply(dd$partitions[ind, ], 1, function(x) all(x == c0))) + 1
  
  y_ticks = c(0, rev(cumsum(rev(prova[,1]))[rev(c(1,16, 41, 51, 52))]))
  
  p = ggplot(df, aes(x = Var2, y = value)) + geom_area(aes(group = Var1, fill = Var1), col = "grey20", fill = ifelse(df$Var1 == c0_pos, col2, "grey90")) +
    theme_minimal() +
    scale_x_continuous(breaks= psi_val, limits = c(0, 3), expand = c(0, 0)) +
    theme(plot.margin=grid::unit(c(1,1,1,1), "cm")) +
    # scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0)) + 
    scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0),
                       sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                           labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
    theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
          axis.ticks.y.left = element_line(size = 0.6, colour = "grey20", linetype = "dotted"),
          axis.ticks.length.y.left  = unit(2, "cm"),
          axis.ticks.length.x.bottom  = unit(2, "cm"),
          axis.text.y.right = element_text( size = rel(1.6)),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
          axis.title.x=element_text(size=20,face="bold")) +
    xlab(expression(psi)) 
  
  p
  
  pt = p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.4,ymax = 0.4,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.84, ymax = 0.84,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.93,ymax = 0.93,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "4/5 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.99,ymax = 0.99, xmin = -1,xmax = 0) 
  #   annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.92,ymax = 0.92,  xmin = -1,xmax = 0) +
  
  gt <- ggplot_gtable(ggplot_build(pt))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid_p = grid.arrange(gt)

  ggsave(grid_p, file = paste0("CP_py_", i, "_sigma05_ekDP.pdf"), device = "pdf", height = 14, width = 16, unit = "cm", dpi = 300)
  
  
}


## Expected number of clusters equal to 4

for(i in 1:5){
  c0 = c0_list[[i]]
  
  dd = sharedLogit::dist_from(c0, TRUE)
  tt = table(dd$distances)
  
  psi_val = c(0,0.5,1,1.5,2,2.5, 3)
  num = sapply(psi_val, function(x) p_clust(dd$distances, psi = x))
  alpha = BNPmix::PYcalibrate(4, 5, discount = sigma)$strength
  
  p_dir = apply(dd$partitions, 1 ,function(x) py_eppf(x, alpha, sigma))
  
  
  ### ORDER PARTITIONS
  n_blocks = apply(dd$partition, 1, function(x) length(unique(x)))
  biggest_block = apply(dd$partition, 1, function(x) max(table(x)))
  ind = order(as.numeric(paste(n_blocks, biggest_block, sep = "")), decreasing = TRUE)
  
  probs = apply(num, 2, function(x) (p_dir*x)/sum(p_dir*x))
  
  prova = probs[ind, ]
  colnames(prova) = psi_val
  df = reshape2::melt(prova)
  
  c0_pos = 52 - which(apply(dd$partitions[ind, ], 1, function(x) all(x == c0))) + 1
  
  y_ticks = c(0, rev(cumsum(rev(prova[,1]))[rev(c(1,16, 41, 51, 52))]))
  
  p = ggplot(df, aes(x = Var2, y = value)) + geom_area(aes(group = Var1, fill = Var1), col = "grey20", fill = ifelse(df$Var1 == c0_pos, col2, "grey90")) +
    theme_minimal() +
    scale_x_continuous(breaks= psi_val, limits = c(0, 3), expand = c(0, 0)) +
    theme(plot.margin=grid::unit(c(1,1,1,1), "cm")) +
    # scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0)) + 
    scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0),
                       sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                           labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
    theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
          axis.ticks.y.left = element_line(size = 0.6, colour = "grey20", linetype = "dotted"),
          axis.ticks.length.y.left  = unit(2, "cm"),
          axis.ticks.length.x.bottom  = unit(2, "cm"),
          axis.text.y.right = element_text( size = rel(1.6)),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
          axis.title.x=element_text(size=20,face="bold")) +
    xlab(expression(psi)) 
  
  p
  
  pt = p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.1,ymax = 0.1,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.4, ymax = 0.4,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.8,ymax = 0.8, xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.91,ymax = 0.91,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "5 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.98,ymax = 0.98,  xmin = -1,xmax = 0)
  
  gt <- ggplot_gtable(ggplot_build(pt))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid_p = grid.arrange(gt)
  
  ggsave(grid_p, file = paste0("CP_py_", i, "_sigma05_ek4.pdf"), device = "pdf", height = 14, width = 16, unit = "cm", dpi = 300)
  
  
}



###
sigma = 0.75


for(i in 1:5){
  c0 = c0_list[[i]]
  
  dd = sharedLogit::dist_from(c0, TRUE)
  tt = table(dd$distances)
  
  psi_val = c(0,0.5,1,1.5,2,2.5, 3)
  num = sapply(psi_val, function(x) p_clust(dd$distances, psi = x))
  alpha = BNPmix::PYcalibrate(log(5), 5, discount = sigma)$strength
  
  p_dir = apply(dd$partitions, 1 ,function(x) py_eppf(x, alpha, sigma))
  
  
  ### ORDER PARTITIONS
  n_blocks = apply(dd$partition, 1, function(x) length(unique(x)))
  biggest_block = apply(dd$partition, 1, function(x) max(table(x)))
  ind = order(as.numeric(paste(n_blocks, biggest_block, sep = "")), decreasing = TRUE)
  
  probs = apply(num, 2, function(x) (p_dir*x)/sum(p_dir*x))
  
  prova = probs[ind, ]
  colnames(prova) = psi_val
  df = reshape2::melt(prova)
  
  c0_pos = 52 - which(apply(dd$partitions[ind, ], 1, function(x) all(x == c0))) + 1
  
  y_ticks = c(0, rev(cumsum(rev(prova[,1]))[rev(c(1,16, 41, 51, 52))]))
  
  p = ggplot(df, aes(x = Var2, y = value)) + geom_area(aes(group = Var1, fill = Var1), col = "grey20", fill = ifelse(df$Var1 == c0_pos, col2, "grey90")) +
    theme_minimal() +
    scale_x_continuous(breaks= psi_val, limits = c(0, 3), expand = c(0, 0)) +
    theme(plot.margin=grid::unit(c(1,1,1,1), "cm")) +
    # scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0)) + 
    scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0),
                       sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                           labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
    theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
          axis.ticks.y.left = element_line(size = 0.6, colour = "grey20", linetype = "dotted"),
          axis.ticks.length.y.left  = unit(2, "cm"),
          axis.ticks.length.x.bottom  = unit(2, "cm"),
          axis.text.y.right = element_text( size = rel(1.6)),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
          axis.title.x=element_text(size=20,face="bold")) +
    xlab(expression(psi)) 
  
  p
  
  pt = p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.3,ymax = 0.3,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.84, ymax = 0.84,  xmin = -1,xmax = 0) +
    # annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.93,ymax = 0.93,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "3-5 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.99,ymax = 0.99, xmin = -1,xmax = 0) 
  #   annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.92,ymax = 0.92,  xmin = -1,xmax = 0) +
  
  gt <- ggplot_gtable(ggplot_build(pt))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid_p = grid.arrange(gt)
  
  ggsave(grid_p, file = paste0("CP_py_", i, "_sigma075_ekDP.pdf"), device = "pdf", height = 14, width = 16, unit = "cm", dpi = 300)
  
  
}



###
sigma = 0.25

for(i in 1:5){
  c0 = c0_list[[i]]
  
  dd = sharedLogit::dist_from(c0, TRUE)
  tt = table(dd$distances)
  
  psi_val = c(0,0.5,1,1.5,2,2.5, 3)
  num = sapply(psi_val, function(x) p_clust(dd$distances, psi = x))
  alpha = BNPmix::PYcalibrate(log(5), 5, discount = sigma)$strength
  
  p_dir = apply(dd$partitions, 1 ,function(x) py_eppf(x, alpha, sigma))
  
  
  ### ORDER PARTITIONS
  n_blocks = apply(dd$partition, 1, function(x) length(unique(x)))
  biggest_block = apply(dd$partition, 1, function(x) max(table(x)))
  ind = order(as.numeric(paste(n_blocks, biggest_block, sep = "")), decreasing = TRUE)
  
  probs = apply(num, 2, function(x) (p_dir*x)/sum(p_dir*x))
  
  prova = probs[ind, ]
  colnames(prova) = psi_val
  df = reshape2::melt(prova)
  
  c0_pos = 52 - which(apply(dd$partitions[ind, ], 1, function(x) all(x == c0))) + 1
  
  y_ticks = c(0, rev(cumsum(rev(prova[,1]))[rev(c(1,16, 41, 51, 52))]))
  
  p = ggplot(df, aes(x = Var2, y = value)) + geom_area(aes(group = Var1, fill = Var1), col = "grey20", fill = ifelse(df$Var1 == c0_pos, col2, "grey90")) +
    theme_minimal() +
    scale_x_continuous(breaks= psi_val, limits = c(0, 3), expand = c(0, 0)) +
    theme(plot.margin=grid::unit(c(1,1,1,1), "cm")) +
    # scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0)) + 
    scale_y_continuous(breaks= y_ticks, limits = c(0, 1), expand = c(0, 0),
                       sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                           labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
    theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
          axis.ticks.y.left = element_line(size = 0.6, colour = "grey20", linetype = "dotted"),
          axis.ticks.length.y.left  = unit(2, "cm"),
          axis.ticks.length.x.bottom  = unit(2, "cm"),
          axis.text.y.right = element_text( size = rel(1.6)),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
          axis.title.x=element_text(size=20,face="bold")) +
    xlab(expression(psi)) 
  
  p
  
  pt = p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.2,ymax = 0.2,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.6, ymax = 0.6,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.93,ymax = 0.93,  xmin = -1,xmax = 0) +
    annotation_custom(grob = textGrob(label = "4/5 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.99,ymax = 0.99, xmin = -1,xmax = 0) 
  #   annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.85)), ymin = 0.92,ymax = 0.92,  xmin = -1,xmax = 0) +
  
  gt <- ggplot_gtable(ggplot_build(pt))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid_p = grid.arrange(gt)
  
  ggsave(grid_p, file = paste0("CP_py_", i, "_sigma025_ekDP.pdf"), device = "pdf", height = 14, width = 16, unit = "cm", dpi = 300)
  
  
}
