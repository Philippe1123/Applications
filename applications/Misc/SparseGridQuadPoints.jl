using SparseGrids
using DelimitedFiles
using FastGaussQuadrature
using Revise

Dim=3
n=4

nodes, weights = sparsegrid(Dim,n,gausslegendre)

pathToNodes="/home/philippe/Desktop/nodesSparse.txt"
pathToWeights="/home/philippe/Desktop/weigthsSparse.txt"


fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")

writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)

close(fileNodes)
close(fileWeigths)


nodes, weights = gausslegendre(29)

pathToNodes="/home/philippe/Desktop/nodesGl.txt"
pathToWeights="/home/philippe/Desktop/weigthsGl.txt"


fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")

writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)

close(fileNodes)
close(fileWeigths)






nodes, weights = sparsegrid(2,1,gausslegendre)
pathToNodes="/home/philippe/Desktop/nodesSparse_1.txt"
pathToWeights="/home/philippe/Desktop/weigthsSparse_1.txt"
fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")
writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)
close(fileNodes)
close(fileWeigths)


nodes, weights = sparsegrid(2,2,gausslegendre)
pathToNodes="/home/philippe/Desktop/nodesSparse_2.txt"
pathToWeights="/home/philippe/Desktop/weigthsSparse_2.txt"
fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")
writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)
close(fileNodes)
close(fileWeigths)


nodes, weights = sparsegrid(2,3,gausslegendre)
pathToNodes="/home/philippe/Desktop/nodesSparse_3.txt"
pathToWeights="/home/philippe/Desktop/weigthsSparse_3.txt"
fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")
writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)
close(fileNodes)
close(fileWeigths)

nodes, weights = sparsegrid(2,4,gausslegendre)
pathToNodes="/home/philippe/Desktop/nodesSparse_4.txt"
pathToWeights="/home/philippe/Desktop/weigthsSparse_4.txt"
fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")
writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)
close(fileNodes)
close(fileWeigths)


nodes, weights = sparsegrid(2,5,gausslegendre)
pathToNodes="/home/philippe/Desktop/nodesSparse_5.txt"
pathToWeights="/home/philippe/Desktop/weigthsSparse_5.txt"
fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")
writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)
close(fileNodes)
close(fileWeigths)


nodes, weights = sparsegrid(2,6,gausslegendre)
pathToNodes="/home/philippe/Desktop/nodesSparse_6.txt"
pathToWeights="/home/philippe/Desktop/weigthsSparse_6.txt"
fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")
writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)
close(fileNodes)
close(fileWeigths)


nodes, weights = sparsegrid(2,7,gausslegendre)
pathToNodes="/home/philippe/Desktop/nodesSparse_7.txt"
pathToWeights="/home/philippe/Desktop/weigthsSparse_7.txt"
fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")
writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)
close(fileNodes)
close(fileWeigths)


nodes, weights = sparsegrid(2,8,gausslegendre)
pathToNodes="/home/philippe/Desktop/nodesSparse_8.txt"
pathToWeights="/home/philippe/Desktop/weigthsSparse_8.txt"
fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")
writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)
close(fileNodes)
close(fileWeigths)


nodes, weights = sparsegrid(2,9,gausslegendre)
pathToNodes="/home/philippe/Desktop/nodesSparse_9.txt"
pathToWeights="/home/philippe/Desktop/weigthsSparse_9.txt"
fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")
writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)
close(fileNodes)
close(fileWeigths)


nodes, weights = sparsegrid(2,10,gausslegendre)
pathToNodes="/home/philippe/Desktop/nodesSparse_10.txt"
pathToWeights="/home/philippe/Desktop/weigthsSparse_10.txt"
fileNodes=open(pathToNodes,"w")
fileWeigths=open(pathToWeights,"w")
writedlm(fileNodes,nodes)
writedlm(fileWeigths,weights)
close(fileNodes)
close(fileWeigths)
