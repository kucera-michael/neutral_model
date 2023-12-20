n.row = 10 #size of matrix (# of rows and or columns) - J

p.outside = .05 #p of outside pool replacement

p.speciation =  .1 #p of speciation (of specie)
p.ispeciation = .4 #p of speciation (of single idividuals)

p.extinction = 0 #p of extinction (of specie)
p.iextinciton = .7 #p of extinction (of single individuals)
t.extinction = .5 #threshold for species being able to expand when there is extinction

p.death = .1 #p of death

n.generations = 1000 #n of generations
n.runs = 10 #n of runs - r

list.runs = list() #collector for iterations

# add parallel here

for (run in 1:n.runs){
  
  list.distribution = list() #collector for o (order of distribution)
  for (outside.order in c(F, T)){ #order of distribution for inside pool replacement (if T inverse, if F normal)
    n.species = 1 #species - N
    list.table = c() #collector for table(M) (distribution of species)
    list.count = c() #collector for max(M) (number of species that were)
    list.matrix = list() #collector for M (species matrix)
  
    species.matrix = matrix(1, nrow = n.row, ncol = n.row) #generate matrix
  
    #set.seed(run)
    for (generation in 1:n.generations){ #cycle generations
    
      generation.table = as.numeric(table(species.matrix))/(n.row^2) #prob distribution of classes
      generation.classes = sort(unique(as.numeric(species.matrix))) #ordered distinct classes
    
      if (sample(c(0,1), 1, prob = c(1-p.speciation, p.speciation))){ 
        #decide if there is speciation in this generation
        speciation.class = sample(generation.classes, 1)  # select speciated class
        speciation.event = T
        n.species = n.species + 1
      } else {
          speciation.event = F
      }
      
      if (sample(c(0,1), 1, prob = c(1-p.extinction, p.extinction))){ 
        #decide if there is extinction in this generation
        extinction.class = sample(generation.classes, 1)  # select speciated class - s
        extinction.event = T
        if (p.iextinciton<t.extinction){
          generation.table = as.numeric(table(species.matrix))[-which(generation.classes==extinction.class)]/((n.row^2)-as.numeric(table(species.matrix)[which(sort(unique(as.numeric(species.matrix)))==extinction.class)]))
          generation.classes = generation.classes[!generation.classes==extinction.class]
        }
      } else {
        extinction.class = F
        extinction.event = F
      }
      
  
    for (i in 1:nrow(species.matrix)){
      for (j in 1:ncol(species.matrix)){ #cycle cells
        if (speciation.event) { #if there is speciation 
          if (species.matrix[i, j]==speciation.class){ #if the cell has speciation class
            if (sample(c(0,1), 1, prob = c(1-p.ispeciation, p.ispeciation))){ #decide if the individual speciates
              species.matrix[i, j] = n.species
            }
          }
        } else {
          #do nothing
        }
        
        if (sample(c(0, 1), 1, prob = c(1-p.death, p.death)) | (extinction.event & species.matrix[i,j]==extinction.class & sample(c(0,1), 1, prob = c(1-p.iextinciton, p.iextinciton)))){ #decide death (natural and extinction)
          if (sample(c(0, 1), 1, prob = c(1-p.outside, p.outside))){ #decide if the dead individual is replaced from outside pool
            if (speciation.event){
              species.matrix[i, j] = n.species + 1
            } else {
              species.matrix[i, j] = n.species
            }
          } else {
            species.matrix[i, j] = sample(generation.classes[order(generation.table, decreasing = outside.order)], 1, prob = generation.table[order(generation.table)]) #pick the class from inside pool
  
          }
        }
      if (is.na(species.matrix[i, j])){
        species.matrix[i, j] = n.species + 1
      }
      }
    }
      
  if (max(species.matrix) < n.species){
    n.species = n.species + 1 #when a species dies out right after speciation
  } else {
      n.species = max(species.matrix) #???
  }

  
  list.table = append(list.table, length(table(species.matrix))) #save table of class counts (# of present species)
  list.count = append(list.count, n.species)#max(unique(as.numeric(M)))) #save count of all classes that have appeared (all time # of species)
  list.matrix[[generation]] = species.matrix #save matrix
}
    list.distribution[[as.integer(outside.order)+1]] = list(list.matrix, list.table, list.count) # O[[2]] has inverted probability data
}
  list.runs[[run]] = list.distribution
}

for (z in 1:2){
plot(list.runs[[1]][[2]][z+1][[1]], type='l', xlab = 'generation', ylab = '#species', main=c('present # of species','all time # of species')[z], col='red')
for (x in 1:10){
lines(list.runs[[x]][[2]][z+1][[1]], col='red')
lines(list.runs[[x]][[1]][z+1][[1]], col='blue')
}
legend(x = 'topleft', legend = c('inverse', 'count'), fill = c('red', 'blue'))
}

library(animation)

set.seed(5) #need to get rid of white
vis.range = c(1:25)
color.sample = sample(rainbow(max(list.runs[[1]][[1]][3][[1]][max(vis.range)], list.runs[[1]][[2]][3][[1]][max(vis.range)]))) # color palette - R


for (x in 1:2){
  saveGIF({
    for (i in vis.range){
      M = list.runs[[1]][[x]][[1]][[i]] # we are picking the first runs of both inverse and count (x) to the 25th generation
      grid = expand.grid(1:nrow(M),1:ncol(M))
      grid = transform(grid, z = as.numeric(M))
      plot(grid$Var1, grid$Var2, col=color.sample[grid$z], pch=19, cex=2, axes=F, xlab='', ylab='', main=paste0('generation: ', i)) # the visualization is flipped position 1,1 is bottom left (top left in matrix form)
      legend('topleft', horiz=T, legend=sort(unique(grid$z)), inset=c(.025,1.025), xpd=T, fill=color.sample[sort(unique(grid$z))])
      }  
}, paste0(c('count', 'inverse')[x],'_neutral_model.gif')) # 1 is count 2 is inverse
}

