setwd('~/Desktop/FYP/Count')

#Import data
dataDir = c('Aero_1', 'Aero_2', 'Aero_3', 'Ana_1', 'Ana_2', 'Ana_3')

data = read.table('Aero_1count.csv', sep = '\t', header = F)
for(i in 2:6)
{
  tem = read.table(paste(dataDir[i], 'count.csv', sep = ''), sep = '\t', header = F)
  data = cbind(data, i = tem[,2])
  tem = NULL
}
colnames(data) = c('gene', dataDir)

datastat = data[-(1:4436),]
data = data[1:4436,]

#Plot for viewing
colarr = c('red', 'blue', 'yellow', 'green', 'cyan')
plot(data$Aero_1)
for(i in 2:6)
{
  points(data[,i], col = colarr[i])
}
#assign(dataDir[i], read.csv('Aero_1count.csv', sep = '\t', header = F))