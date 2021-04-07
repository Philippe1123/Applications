using ScatteredInterpolation
using MultilevelEstimators
using Revise
using DelimitedFiles

Nodes=Dict()

pathToInterm=joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE","data")
folder = string(pathToInterm,"/Interm/Beam") # for report
folder_with_elements=string(pathToInterm,"/Mesh/Beam")
println(folder)

type="h"
indices =[Level(0) ;Level(1)]

PathE=string(string(folder,string("/","E_1_1")))
Handle_E=open(PathE)
E=readdlm(Handle_E);
close(Handle_E)
println(E)

for id in indices
println(id)
SizeId=length(id)
Access="_";
for i=1:SizeId
Access=string(Access,id[i])
println(Access)
end







PathNodes=string(string(folder_with_elements,string("/",type,"_refinement/Nodes_L",Access)),".txt")
Handle_Nodes=open(PathNodes)
Node=readdlm(Handle_Nodes);
Node=Node[:,2:3]
Nodes[id]=Node
#println(Nodes[id])
close(Handle_Nodes)
println("h-ref")


end

itp = ScatteredInterpolation.interpolate(NearestNeighbor(), Nodes[Level(1)]', E);
println(Nodes[Level(1)])
Zc=evaluate(itp, Nodes[Level(0)]')
println(Zc)
