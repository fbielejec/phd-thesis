################
#---PACKAGES---#
################
require("sqldf")
require("ggplot2")
require("grid")

#################
#---FUNCTIONS---#
#################
sqlSort <- function(table) {
  query = paste("SELECT * FROM", table, "ORDER BY seed ASC;")
  sortedTable = sqldf(query)
  return(sortedTable)
}

# true coverage
coverage <- function(data, true) {
  
  coverage = 0
  n = dim(data)[1]
  for(i in 1 : n) {
    if(data$hpdLower[i] < true && data$hpdUpper[i] > true) {
      
      coverage = coverage + 1
      
    }#END: hpd check
  }#END: row loop
  
  coverage = coverage / n
  cat(coverage)
}

# mean signed difference
msd <- function(data, true) {
  
  bias = c()
  n = dim(data)[1]
  for(i in 1 : n) {
    
    bias[i] = data$mean[i] - true
    
  }#END: row loop
  
  cat(mean(bias))
}

# mean value of the squared deviations of the predictions from the true values 
mse <- function(data, true) {
  
  sqError = c()
  n = dim(data)[1]
  for(i in 1 : n) {
    
    sqError[i] = (data$mean[i] - true)^2
    
  }#END: row loop
  
  cat(mean(sqError))
}

#################
#---SORT DATA---#
#################
kappa_1 = read.table("kappa.1", header = TRUE, comment.char = "*")[, -c(10,11)]
kappa_2 = read.table("kappa.2", header = TRUE, comment.char = "*")[, -c(10,11)]
kappa_3 = read.table("kappa.3", header = TRUE, comment.char = "*")[, -c(10,11)]

kappa_1 <- sqlSort("kappa_1")
kappa_2 <- sqlSort("kappa_2")
kappa_3 <- sqlSort("kappa_3")

#################
#---HPD PLOTS---#
#################
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
vplayout <- function(x, y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}

p <- ggplot(kappa_1)
p <- p + aes(x = seq(min(mean), max(mean), length.out = 100) )
p <- p + geom_errorbar(  aes(ymin = hpdLower, ymax = hpdUpper, x = seed), size = .5, width = .75) 
p <- p +  geom_point(aes(y = mean, x = seed), size = .3)
p <- p + scale_x_continuous('Seed')
p <- p + scale_y_continuous('kappa 1')
p <- p + geom_hline(yintercept = 1)
p <- p + coord_flip()
p <- p + theme_bw()
print(p, vp = vplayout(1, 1))

p <- ggplot(kappa_2)
p <- p + aes(x = seq(min(mean), max(mean), length.out = 100) )
p <- p + geom_errorbar(  aes(ymin = hpdLower, ymax = hpdUpper, x = seed), size = .5, width = .75) 
p <- p +  geom_point(aes(y = mean, x = seed), size = .3)
p <- p + scale_x_continuous('Seed')
p <- p + scale_y_continuous('kappa 2')
p <- p + geom_hline(yintercept = 10)
p <- p + coord_flip()
p <- p + theme_bw()
print(p, vp = vplayout(1, 2))

p <- ggplot(kappa_3)
p <- p + aes(x = seq(min(mean), max(mean), length.out = 100) )
p <- p + geom_errorbar(  aes(ymin = hpdLower, ymax = hpdUpper, x = seed), size = .5, width = .75) 
p <- p +  geom_point(aes(y = mean, x = seed), size = .3)
p <- p + scale_x_continuous('Seed')
p <- p + scale_y_continuous('kappa 3')
p <- p + geom_hline(yintercept = 1)
p <- p + coord_flip()
p <- p + theme_bw()
print(p, vp = vplayout(1, 3))

##################
#---STATISTICS---#
##################
coverage(kappa_1, true = 1.0)
coverage(kappa_2, true = 10.0)
coverage(kappa_3, true = 1.0)

msd(kappa_1, true = 1.0)
msd(kappa_2, true = 10.0)
msd(kappa_3, true = 1.0)

mse(kappa_1, true = 1.0)
mse(kappa_2, true = 10.0)
mse(kappa_3, true = 1.0)





