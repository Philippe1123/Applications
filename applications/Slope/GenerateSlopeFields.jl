push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","MultilevelEstimators","applications","SPDE")))
push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))

using ScatteredInterpolation
using MultilevelEstimators
using Revise
using DelimitedFiles
using Statistics
using GaussianRandomFields
#using FieldTransformation
using Random


function compute_centers(p,t)
    d = size(p, 2)
    vec_t = vec(t)
    size_t = size(t)

    pts = Array{Float64}(undef, size(t, 1), d)
    @inbounds for i in 1:d
        x = reshape(p[vec_t, i], size_t)
        mean!(view(pts, :, i), x)
    end
    pts
end

function expfield(Z::Vector{T})where {T<:Real}

Mean=log(6000)
σ=0.05

Z=Z.*σ
Z=Z.+Mean
Z=exp.(Z)

return Z
end



pathToInterm=joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE","data")
folder = string(pathToInterm,"/Interm/Beam") # for report
folder_with_elements=string(pathToInterm,"/Mesh/Slope")
#println(folder)

folder_with_elements=string(pathToInterm,"/Mesh/Slope")

#println(folder_with_elements)

p=2
#exp_field = GaussianRandomFields.Exponential(1.0,σ=1,p=2)
#exp_field = GaussianRandomFields.Exponential(1.0,σ=1,p=2)
exp_field = GaussianRandomFields.Matern(1.0,0.6,σ=1,p=2)
#exp_field = GaussianRandomFields.SquaredExponential(5.0,σ=1,p=2)
Nterms=10
lv=1
lvl_1=lv-1
type="p"
if(type=="h")
    typenum=0
    isHigerOrderRefinement=false

else
    typenum=1
    isHigerOrderRefinement=true

end

println(isHigerOrderRefinement)


if(lv==0)
indices =[Level(lv)]
else
indices =[Level(lvl_1) ;Level(lv)]
end
index=Level(lv)
index_1=lv-1
index_1=Level(index_1)


Elements=Dict()
Nodes=Dict()
Centers=Dict()
Nodes_Sorted=Dict()

for id in indices
println(id)
SizeId=length(id)
Access="_";
for i=1:SizeId
Access=string(Access,id[i])
end


PathElement=string(string(folder_with_elements,string("/",type,"_refinement/Elements_L",Access)),".txt")

println(PathElement)

Handle_Elements=open(PathElement)
Element=readdlm(Handle_Elements,Int);
Element=Element[:,5:7]
Elements[id]=Element

println(isHigerOrderRefinement)

if(isHigerOrderRefinement==false)
#PathNodes_Sorted=string(string(folder_with_elements,string("/",type,"_refinement/Nodes_Sorted_L",Access)),".txt")
PathNodes=string(string(folder_with_elements,string("/",type,"_refinement/Nodes_L",Access)),".txt")

PathElementCenters=string(string(folder_with_elements,string("/",type,"_refinement/ElementsCenter_L",Access)),".txt")
Handle_PathElementCenters=open(PathElementCenters,"w")

#Handle_PathElementCenters=open(PathElementCenters)
Handle_Nodes=open(PathNodes)
Node=readdlm(Handle_Nodes);
Node=Node[:,2:3]
Nodes[id]=Node

#Handle_Nodes_Sorted=open(PathNodes_Sorted)
#Node_Sorted=readdlm(Handle_Nodes_Sorted);
#Node_Sorted=Node_Sorted[:,2:3]
#Nodes_Sorted[id]=Node_Sorted
#close(Handle_Nodes_Sorted)


Center=compute_centers(Node,Element)
Centers[id]=Center
Centers[id]=[Centers[id]; Nodes[id]]
#Centers[id]=readdlm(Handle_PathElementCenters);
writedlm(Handle_PathElementCenters,Centers[id])
#Centers[id]=readdlm(Handle_PathElementCenters);
close(Handle_PathElementCenters)
close(Handle_Nodes)
Nodes[id]=Centers[id]################################################################################################################################################3
println("h-ref")
else

    PathNodesElements=string(string(folder_with_elements,string("/",type,"_refinement/Nodes_L",Access)),".txt")
    Handle_Nodes_Elements=open(PathNodesElements)
    Nodes_Elements=readdlm(Handle_Nodes_Elements);
Nodes_Elements=Nodes_Elements[:,2:3]

  PathNodes=string(string(folder_with_elements,string("/",type,"_refinement/GaussPoints_L",Access)),".txt")
  Handle_Nodes=open(PathNodes)
  Node=readdlm(Handle_Nodes);
  Node=Node[:,2:3]
  Nodes[id]=Node
  Nodes[id]=[Nodes[id];Nodes_Elements]
  close(Handle_Nodes)
  println("p-ref and hp-ref")

end


close(Handle_Elements)


end





#exp_field = GaussianRandomFields.Exponential(0.3,σ=1,p=p)
println("P of covar equals")
println(p)

cov = CovarianceFunction(2,exp_field)


# Note fields are generated as vy x vy thus 160 x 40 example the solve transposes this

grfs=Dict()
for index in indices
println(index)
Random.seed!(123456)
#println(Nodes[index])
Nods=Nodes[index];
Nodx=Nods[:,1]
Nody=Nods[:,2]
#println(Nodx)
grfs[index] = GaussianRandomField(cov,KarhunenLoeve(Nterms),Nodes[index],Elements[index],quad=GaussLegendre())
#println(grfs[index])
end

println(folder_with_elements)


field_res=string(folder_with_elements,"/Testing_Fields")

Nodes_Fine=Nodes[index]
Elements_Fine=Elements[index]
Zf = GaussianRandomFields.sample(grfs[index]) # compute GRF
#Zf=FieldTransformation.Transform_L(Zf)
Zf=expfield(Zf)
Zf_res=open(string(field_res,"/Fine"),"w")
#Zf=grfs[index].data.eigenval

writedlm(Zf_res,Zf)
close(Zf_res)

if(lv>0)
Nodes_Coarse=Nodes[index_1]
#println(size(Zf))
#println(size(Nodes_Fine))
itp = ScatteredInterpolation.interpolate(NearestNeighbor(), Nodes_Fine', Zf);
Zc=evaluate(itp, Nodes_Coarse')
Zc=vec(Zc)

Zc_res=open(string(field_res,"/Coarse"),"w")
writedlm(Zc_res,Zc)
close(Zc_res)

data_type=open(string(field_res,"/Data"),"w")
writedlm(data_type,typenum)
close(data_type)
end
data_type=open(string(field_res,"/lvl"),"w")
writedlm(data_type,lv)
close(data_type)
