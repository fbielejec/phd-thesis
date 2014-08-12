
n = 12
indices = data.frame(x = rep(NA, n), y = rep(NA, n), label = rep(NA, n), type = rep(NA, n))

i = 1
index = 7
indices$x[i] = (phylo$xend[index] - phylo$x[index])/2 + phylo$x[index]
indices$y[i] = phylo$y[index]
indices$label[i] = "t[1]"
indices$type[i] = "branch"

i = 2
index = 8
indices$x[i] = (phylo$xend[index] - phylo$x[index])/2 + phylo$x[index]
indices$y[i] = phylo$y[index]
indices$label[i] = "t[2]"
indices$type[i] = "branch"

i = 3
index = 9
indices$x[i] = (phylo$xend[index] - phylo$x[index])/2 + phylo$x[index]
indices$y[i] = phylo$y[index]
indices$label[i] = "t[3]"
indices$type[i] = "branch"

i = 4
index = 11
indices$x[i] = (phylo$xend[index] - phylo$x[index])/2 + phylo$x[index]
indices$y[i] = phylo$y[index]
indices$label[i] = "t[4]"
indices$type[i] = "branch"

i = 5
index = 12
indices$x[i] = (phylo$xend[index] - phylo$x[index])/2 + phylo$x[index]
indices$y[i] = phylo$y[index]
indices$label[i] = "t[5]"
indices$type[i] = "branch"

i = 6
index = 6
indices$x[i] = (phylo$xend[index] - phylo$x[index])/2 + phylo$x[index]
indices$y[i] = phylo$y[index]
indices$label[i] = "t[6]"
indices$type[i] = "branch"

i = 7
index = 5
indices$x[i] = (phylo$xend[index] - phylo$x[index])/2 + phylo$x[index]
indices$y[i] = phylo$y[index]
indices$label[i] = "t[7]"
indices$type[i] = "branch"

i = 8
index = 10
indices$x[i] = (phylo$xend[index] - phylo$x[index])/2 + phylo$x[index]
indices$y[i] = phylo$y[index]
indices$label[i] = "t[8]"
indices$type[i] = "branch"

i = 9
index = 1
indices$x[i] = phylo$x[index]
indices$y[i] = (phylo$yend[index] - phylo$y[index])/2 + phylo$y[index]
indices$label[i] = "x[0]"
indices$type[i] = "node"

i = 10
xindex = 3
yindex = 6
indices$x[i] = phylo$x[xindex]
indices$y[i] = phylo$y[yindex]
indices$label[i] = "x[6]"
indices$type[i] = "node"

i = 11
xindex = 4
yindex = 10
indices$x[i] = phylo$x[xindex]
indices$y[i] = phylo$y[yindex]
indices$label[i] = "x[7]"
indices$type[i] = "node"

i = 12
xindex = 2
yindex = 5
indices$x[i] = phylo$x[xindex]
indices$y[i] = phylo$y[yindex]
indices$label[i] = "x[8]"
indices$type[i] = "node"

