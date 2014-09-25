#################
#---SCRAPBOOK---#
#################


JC69 <- function (lambda, t) {
  
  p0 = 1/4 + 3/4 * exp(-4*lambda*t)
  p1 = 1/4 - 1/4 * exp(-4*lambda*t)
  
  trans_prob = matrix(c(p0, p1, p1, p1,
                        p1, p0, p1, p1,
                        p1, p1, p0, p1,
                        p1, p1, p1, p0
  ), nrow = 4, ncol = 4, byrow = T)
  
  colnames(trans_prob) <- rownames(trans_prob) <- c("T", "C", "A", "G")
  
  return(trans_prob)
}



N <- 100
meanWaitTime <- 3
times = rexp(N, rate = 1 / meanWaitTime)

lambda = 2
states = array("", c(N, 1))
states[1] = "T"

for (i in 2 : N) {
  time = times[i]
  states[i] = sample(c("T", "C", "A", "G"), size = 1, prob = JC69(lambda, time)[state[i - 1], ])
}


fx <- function(x, delta) {
  return( delta * dnorm(x, 7, 0.5^2) + (1 - delta) * dnorm(x, 10, 0.5^2) )
}

rf <- function(n, delta) {
  u = I(runif(n) < delta)
  y = u * rnorm(n, 7, 0.5) + (1 - u) *
    rnorm(n, 10, 0.5)
  
  return(y)
}

N = 200
delta = 0.7
xs = rf(N, delta)
hist(xs, freq = FALSE, ylim = c(0, 1) )
xgrid <- seq(0, 11, 0.001) 
lines(xgrid, fx(xgrid, delta = .7), col = "red")

A = matrix(c(1,2,3,4), nrow = 2, ncol = 2)

I = matrix(c(1,0,0,1), nrow = 2, ncol = 2)

A %*% I


x <- 
  c(
    0.450165759800979 ,
    0.5226732705308421 ,
    0.02263020875055116 ,
    0.002694186745651535 ,
    0.0013598061297214168 ,
    3.2425993395277794E-4 ,
    1.4765164171136013E-4 ,
    7.321189371973964E-7 ,
    1.6324887117420784E-6 ,
    2.4138011602695543E-6 )

median(x)
quantile(x, probs = 0.5)

7.321189371973964E-7 +
  1.6324887117420784E-6+ 
  2.4138011602695543E-6 

################
#---PACKAGES---#
################
require(ggplot2)
require(reshape2)
require(proto)
require(ape)
require(grid)

require("extrafont")
# for pdf
loadfonts()
# for postscript
loadfonts(device = "postscript")
          
###################
###################
##---FUNCTIONS---##
###################
###################

#---FORTIFY.PHYLO---#
fortify.phylo <- function(phylo) {
  
  Ntip     <- length(phylo$tip.label)
  Nnode    <- phylo$Nnode
  Nedge    <- dim(phylo$edge)[1]
  z        <- reorder(phylo, order = "pruningwise")
  yy       <- numeric(Ntip + Nnode)
  TIPS     <- phylo$edge[phylo$edge[, 2] <= Ntip, 2]
  yy[TIPS] <- 1:Ntip
  
  ans1 <- .C("node_height", as.integer(Ntip), as.integer(Nnode), as.integer(z$edge[, 1]), as.integer(z$edge[, 2]), 
             as.integer(Nedge), 
             as.double(yy), DUP = FALSE, PACKAGE = "ape")
  
  yy <- ans1[[6]]
  
  ans2 <- .C("node_depth_edgelength", as.integer(Ntip),
             as.integer(Nnode), as.integer(z$edge[, 1]), as.integer(z$edge[, 2]),
             as.integer(Nedge), as.double(z$edge.length), double(Ntip + Nnode), DUP
             = FALSE, PACKAGE = "ape")
  
  xx <- ans2[[7]]
  
  edge        <- phylo$edge
  nodes       <- (Ntip + 1) : (Ntip + Nnode)
  x0v         <- xx[nodes]
  y0v         <- y1v <- numeric(Nnode)
  NodeInEdge1 <- vector("list", Nnode)
  
  for (i in nodes) {
    ii      <- i - Ntip
    j       <- NodeInEdge1[[ii]] <- which(edge[, 1] == i)
    tmp     <- range(yy[edge[j, 2]])
    y0v[ii] <- tmp[1]
    y1v[ii] <- tmp[2] 
  }
  
  x0h <- xx[edge[, 1]]
  x1h <- xx[edge[, 2]]
  y0h <- yy[edge[, 2]]
  
  lineh <- data.frame(x = x1h, y = y0h, xend = x0h, yend = y0h)
  linev <- data.frame(x = x0v, y = y1v, xend = x0v, yend = y0v)
  x <- rbind(linev, lineh)
  
  return(x)
}

#---LABEL.PHYLO---#
label.phylo <- function(phylo) {
  
  Ntip <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  Nedge <- dim(phylo$edge)[1]
  z <- reorder(phylo, order = "pruningwise")
  yy <- numeric(Ntip + Nnode)
  TIPS <- phylo$edge[phylo$edge[, 2] <= Ntip, 2]
  yy[TIPS] <- 1:Ntip
  
  ans1 <- .C("node_height_clado", as.integer(Ntip), as.integer(Nnode), as.integer(z$edge[, 1]), 
             as.integer(z$edge[,2]), as.integer(Nedge), double(Ntip + Nnode), as.double(yy), 
             DUP = FALSE, PACKAGE = "ape")
  
  yy <- ans1[[7]]
  
  ans2 <- .C("node_depth_edgelength", as.integer(Ntip), as.integer(Nnode), as.integer(z$edge[, 1]), 
             as.integer(z$edge[, 2]), as.integer(Nedge), as.double(z$edge.length), double(Ntip + Nnode), 
             DUP = FALSE, PACKAGE = "ape")
  
  xx <- ans2[[7]]
  
  data.frame(x = xx[1 : Ntip], y = yy[1 : Ntip], label = phylo$tip.label)
}

#---GEOM_SEGMENT2---#
GeomSegment2 <- proto(ggplot2:::GeomSegment, {
  objname <- "geom_segment2"
  
  draw <- function(., data, scales, coordinates, arrow=NULL, ...) {
    if (is.linear(coordinates)) {
      return(with(coord_transform(coordinates, data, scales),
                  
                  segmentsGrob(x, y, xend, yend, default.units="native",
                               gp = gpar(col=alpha(colour, alpha), lwd = size * .pt,
                                         lty=linetype, lineend = "round"),
                               arrow = arrow)
      ))
    }
  }})

geom_segment2 <- function(mapping = NULL, data = NULL, stat =
                            "identity", position = "identity", arrow = NULL, ...)  {
  GeomSegment2$new(mapping = mapping, data = data, stat = stat,
                   position = position, arrow = arrow, ...)
}

#---GET_EXPR---#
get_expr <- function(eq) {
  return(as.character(as.expression(eq)));  
}

###########################
#---POISSON SAMPLE PATH---#
###########################
set.seed(123456)
n = 10
x = c(0, cumsum(rexp(n)))
y = seq(0, n)

theme1 <- theme(
  axis.line = element_line(colour = "black"),
  axis.text = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  legend.position = "none",
  panel.background = element_rect(size = 1, fill = "white", colour = NA),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

p <- ggplot()
p <- p + geom_step(aes(x = x, y = y), color = "blue")
p <- p + scale_x_continuous(
  breaks = round(x[seq(1, length(x), 2)], 2)
  )
p <- p + scale_y_discrete()
p <- p + theme1
p <- p + xlab("t") + ylab("k")
p

###############################
#---POISSON PROCESS ON TREE---#
###############################
set.seed(123456)

b1 <- rexp(1)
b2 <- sum(c(b1, rexp(1)))
b3 <- sum(c(b1, rexp(1)))
b4 <- rexp(1)

label1 <- "realisation 1"
label2 <- "realisation 2"
label3 <- "realisation 3"

pRealisations <- data.frame(x = c( 0, b1, b2, 0, b1, b3, 0, b4 ), 
                           y = c( 0,1,2,0,1,2,0,1 ), 
                           col = c( rep(label1, 3), rep(label2, 3), rep(label3, 2) )
                           )

tree <- data.frame(
  x = c(0, 0, b1, b1, 0, 0, b1, b1, 0, 0),
  y = c(0, 1, 1, 2, 0, 1, 1, 0, 0, -1),
  xend = c(0, b1, b1, b2, 0, b1, b1, b3, 0, b4),
  yend = c(1, 1, 2, 2, 1, 1, 0, 0, -1, -1),
  col = c( rep(label1, 4), rep(label2, 4), rep(label3, 2) ) 
)

theme1 <- theme(
  axis.line = element_line(colour = "black"),
  axis.title.y = element_blank(),
  axis.text.x = element_text(colour = "black"),
  axis.ticks.x = element_line(colour = "black"),
  legend.position = "none",
  panel.background = element_rect(size = 1, fill = "white", colour = NA),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank()
)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
vplayout <- function(x, y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}

p <- ggplot(pRealisations)
p <- p + geom_step(aes(x = x, y = y, color = col, linetype = col))
p <- p + facet_wrap( ~ col, ncol = 3, scales = "fixed")
p <- p + theme1
p <- p + xlab("t")
print(p, vp = vplayout(1, 1))

p <- ggplot(tree)
p <- p + geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = col, linetype = col), alpha = 1, size = 1.1)
# p <- p + scale_x_continuous(
#   breaks = round(c(b1,b2,b3,b4),2)
# )
p <- p + theme1
p <- p + xlab("t")
print(p, vp = vplayout(2, 1))


################
#---MEP TREE---#
################
n2RateTree <- read.nexus("epoch_2rate_true_data_ordered.tree")

phylo = fortify.phylo(n2RateTree)
phyloLabels = label.phylo(n2RateTree)

yrng = range(phylo$y)
xrng = range(phylo$x)
xlim = 20

# phylo$x = max(phylo$x) - phylo$x
# phylo$xend = max(phylo$x) - phylo$xend
# phyloLabels$x = max(phyloLabels$x) - phyloLabels$x

dotSize = 4
lineSize = 1.5

p <- ggplot()
p <- p + geom_segment2(aes(x = x, y = y, xend = xend, yend = yend), colour = "black", size = lineSize, data = phylo)
p <- p + geom_point(aes(x = x, y = y), size = dotSize, data = phyloLabels)

theme <- theme_update(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.x = element_text(colour = "black"),
  axis.text.x = element_text(colour = "black"),
  axis.title.y = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.line.y = element_blank(),
  legend.position = "none",
  panel.background = element_rect(size = 1, fill = "white", colour = NA),
  plot.background = element_rect(colour = NA),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

# p <- p + scale_fill_brewer(palette = "Blues")

p <- p + theme_set(theme)
# p <- p + scale_x_reverse(limits = c(xlim, 0))
p <- p + xlab("Time")
print(p)


##############################################
#---TIME HETEROGENOUS AND ULTRAMETRIC TREE---#
##############################################
tree1        <- read.tree("data/epoch_2rate_true_data_ordered.newick")
phylo1       <- fortify.phylo(tree1)
phylo1$label <- "Heterogenous"
tree2        <- read.tree("data/epoch_2rate_ultrametric.newick")
phylo2       <- fortify.phylo(tree2)
phylo2$label <- "Ultrametric"

phylo <- rbind(phylo1, phylo2)

phylo1Labels       <- label.phylo(tree1)
phylo1Labels$label <- "Heterogenous"
phylo2Labels       <- label.phylo(tree2)
phylo2Labels$label <- "Ultrametric"

phyloLabels <- rbind(phylo1Labels, phylo2Labels)

dotSize = 3
lineSize = 1.5

p <- ggplot()
p <- p + geom_segment2(aes(x = x, y = y, xend = xend, yend = yend), colour = "black", size = lineSize, data = phylo)
p <- p + geom_point(aes(x = x, y = y), size = dotSize, data = phyloLabels)

theme <- theme_update(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.x = element_text(colour = "black"),
  axis.text.x = element_text(colour = "black"),
  axis.title.y = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.line.y = element_blank(),
  legend.position = "none",
  panel.background = element_rect(size = 1, fill = "white", colour = NA),
  plot.background = element_rect(colour = NA),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text.x = element_text(size = 16, colour = "black")
)

p <- p + theme_set(theme)
p <- p + facet_wrap( ~ label, ncol = 2, scales = "fixed")
p <- p + xlab("Time")
print(p)

#######################################
#---POSTERIOR PRIOR LIKELIHOOD PLOT---#
#######################################
PriorLikePost <- function(N, x, a, b, length.out = 1000) {
  
  theta      = seq(0, .25, length.out = length.out)
  prior      = dbeta(theta, a, b)
  likelihood = dbinom(rep(x, length.out), N, theta)
  posterior  = dbeta(theta, a + x, N - x + b)
  
    prior      = prior / sum(prior)
    likelihood = likelihood / sum(likelihood)
    posterior  = posterior / sum(posterior)
  
  data <- data.frame(
    theta      = theta, 
    prior      = prior, 
    likelihood = likelihood, 
    posterior  = posterior
  )
  
  return(data)
}

# x successes in N trials modelled with binomial model and beta prior on probability of success
data <- PriorLikePost(N = 100, x = 10, a = 2, b = 100 )
melt_data <- melt(data, id.vars = "theta")

p <- ggplot(melt_data) 
p <- p + geom_polygon(aes(x = theta, y = value, group = variable, fill = variable), color = I("black") )

theme <- theme_update(
  axis.title = element_text(colour = "black"),
  axis.text = element_text(colour = "black"),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  panel.background = element_rect(size = 1, fill = "white", colour = NA),
  plot.background = element_rect(colour = NA),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

p <- p + theme_set(theme)
p <- p + scale_fill_discrete("")
p <- p + ylab("") + xlab(expression(""*theta*""))
print(p)


##############################
#---METROPOLIS SAMPLE PLOT---#
##############################
# estimate sequence distance theta under the JC69 model
data = data.frame(x = 10, n = 100)

loglikelihood <- function(theta, data) {
  x = data$x
  n = data$n
  # probability from JC69 model
  p = (3/4) * (1 - exp( -(4/3) * theta ))
  like = dbinom(x, size = n, prob = p)

  return(log(like))
}

prior <- function(x) {
  mean = 0.1
  return(log(dexp(x, rate = 1 / mean)))
}

proposal <- function(xt) {
  
  window = 0.1
  
  r.cand = runif(1, min = xt - window, max = xt + window)
  r.cand = ifelse((r.cand >= 0) & (r.cand <= 1), r.cand, 
                  ifelse(r.cand < 0, 1 + r.cand, 
                         ifelse(r.cand > 1,  r.cand - 1, cat("error"))))
  
  # they will be the same, proposal is symmetrical
  d.cand = dunif(r.cand, min = xt - window, max = xt + window)
  d.curr = dunif(xt, min = xt - window, max = xt + window)
  
  return(list(r.cand = r.cand, d.cand = d.cand, d.curr = d.curr))
}

metropolisHastings <- function(loglikelihood, prior, proposal, startvalue, data, Nsim) {
  
  chain = array(dim = c(Nsim, 1))
  chain[1, ] = startvalue
  for (t in 1 : (Nsim - 1)) {
    
    candidate = proposal(chain[t, ])
    r.candidate = candidate$r.cand
    d.candidate = candidate$d.cand
    d.curr = candidate$d.curr
    
    probab = exp( ( loglikelihood(r.candidate, data) + prior(r.candidate) + log(d.candidate) ) - 
                    ( loglikelihood(chain[t, ], data) + prior(chain[t, ]) + log(d.curr) )
    )
    
    if (runif(1) < probab) {
      chain[t + 1, ] = r.candidate
    } else {
      chain[t + 1, ] = chain[t, ]
    }#END: accept check
    
  }#END: iterations loop
  
  return(chain)
}

Nsim = 10^3
startvalue = runif(1, min = 0, max = 1)
chain = metropolisHastings(loglikelihood, prior, proposal, startvalue, data, Nsim)

#---PLOT---#
thetaHat = -(3/4) * log(1 - (4/3) * (data$x/data$n))
plotData <- data.frame(iteration = c(1 : length(chain)), value = chain)

thetaHat
mean(chain)

theme2 <- theme(
  axis.title = element_text(colour = "black"),
  axis.text = element_text(colour = "black"),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  panel.background = element_rect(size = 1, fill = "white", colour = NA),
  plot.background = element_rect(colour = NA),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

labels <- paste("acceptance rate:", round(length(unique(chain))/Nsim, 2), "\n",
                "sample mean: ", round(mean(chain), 2), sep = " ")

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
vplayout <- function(x, y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}

p <- ggplot(plotData)
p <- p + geom_line(aes(x = iteration, y = value))    
p <- p + geom_hline(aes(yintercept = thetaHat), color = "red")
p <- p + xlab("") + ylab("Chain state")
p <- p + theme2
print(p, vp = vplayout(1, 1))

p <- ggplot(plotData)
p <- p + geom_histogram(aes(x = value, y = ..density..), binwidth = 0.01 )    
p <- p + geom_density(aes(x = value), color = "red", alpha = 0.2)
p <- p + xlab("") + ylab("Density")
p <- p + ggtitle(labels)
p <- p + theme2
print(p, vp = vplayout(2, 1))


############################
#---LIKELIHOOD ON A TREE---#
############################
tree <- read.tree(text="(((T1:0.5,T2:2.0):2.5,T3:4.5):2.0,(T4:2.0,T5:4.0):3.0);")

phylo = fortify.phylo(tree)
phyloLabels = label.phylo(tree)
yrng  = range(phylo$y)
xrng  = range(phylo$x)
xlim = 8
phylo$x = max(phylo$x) - phylo$x
phylo$xend = max(phylo$x) - phylo$xend
phyloLabels$x = max(phyloLabels$x) - phyloLabels$x

source("indices.r")
# epochs = data.frame(start = c(min(phylo$x)), end   = c(xlim), epoch = c("Q1") )

nodeSize = 4
lineSize = 1.5
textSize = 5
dotSize = 12

p <- ggplot()

p <- p + geom_segment2(aes(x = x, y = y, xend = xend, yend = yend), colour = "black", size = lineSize, data = phylo)
p <- p + geom_point(aes(x = x, y = y), size = nodeSize, pch = 21,  colour = "black", fill = "black", data = phyloLabels)
p <- p + geom_text(data = phyloLabels, aes(x = x, y = y, label = c("T", "C", "A", "C", "C")), hjust = -1.0, 
#                    family = "Arial Black", 
                   vjust = 0.5, size = textSize)

p <- p + geom_point(aes(x = x, y = y, colour = type, fill = type), size = dotSize, pch = 21, data = indices)
p <- p + geom_text(aes(x = x, y = y, label = get_expr(label), colour = type), size = textSize, data = indices, parse = TRUE, family="Arial Black")

theme <- theme_update(
  axis.line = element_line(colour = "black"),
  axis.title.y = element_blank(),
  axis.text.x = element_text(colour = "black"),
  axis.ticks.x = element_line(colour = "black"),
  legend.position = "none",
  panel.background = element_rect(size = 1, fill = "white", colour = NA),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank()
)

p <- p + theme_set(theme)
p <- p + scale_x_reverse(limits = c(xlim, 0))
p <- p + xlab(NULL)

p <- p + scale_fill_grey()
gs.pal <- colorRampPalette(c("white","black"))
p <- p + scale_colour_manual(values = gs.pal(2))

print(p)




###################################
#--- EXPONENTIAL APPROXIMATION ---#
###################################
p = 0.01
k = c(0:100)
geomProbs = 1 - pgeom(k, prob = p)
expProbs = 1 - pexp(k, rate = p)

p <- ggplot()
p <- p + geom_rect( aes(xmin = k, xmax = k + 1, ymin = 0, ymax = geomProbs), linetype = 2, fill = "white", color = "black" )
p <- p + geom_line(aes(x = k, y = expProbs), color = "red")
p <- p + theme_bw()
p <- p + xlab("k") + ylab("P(X>k)")
p



