using GaussianRandomFields, DelimitedFiles, LinearAlgebra, Plots, Statistics, Dierckx, MultilevelEstimators, Random
function expfield(Z::Vector{T})where {T<:Real}

Mean=log(6000)
σ=0.05

Z=Z.*σ
Z=Z.+Mean
Z=exp.(Z)

return Z
end

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


    Elements=Dict()
    Nodes=Dict()
   Centers=Dict()

type="h"

folder_with_elements="/home/philippe/Desktop/Slope_Matlab/Mesh_Steep_8/h_ref/Matlab"

for id=0:3

Access=id




PathElement=string(string(folder_with_elements,string("/","Elements_L_",Access)),".txt")
Handle_Elements=open(PathElement)
Element=readdlm(Handle_Elements,Int);
Element=Element[:,5:7]
Elements[id]=Element



PathNodes=string(string(folder_with_elements,string("/","Nodes_L_",Access)),".txt")
PathElementCenters=string(string(folder_with_elements,string("/","ElementsCenter_L_",Access)),".txt")
Handle_PathElementCenters=open(PathElementCenters,"w")
Handle_Nodes=open(PathNodes)
Node=readdlm(Handle_Nodes);
Node=Node[:,2:3]
Nodes[id]=Node
Center=compute_centers(Node,Element)
Centers[id]=Center
writedlm(Handle_PathElementCenters,Center)
close(Handle_PathElementCenters)
close(Handle_Nodes)
close(Handle_Elements)

end
corr_len=3.0
p=2
exp_field = GaussianRandomFields.SquaredExponential(corr_len,σ=1.0,p=p)
println("P of covar equals")
println(p)
cov = CovarianceFunction(2,exp_field)


# all other levels
grfs=Dict()
nterms=101
for id=0:3
 Random.seed!(1234)
 grfs[id] = GaussianRandomField(cov,KarhunenLoeve(nterms),Nodes[id],Elements[id],quad=GaussLegendre())
 println(grfs[id])
 field=expfield(sample(grfs[id]))
 PathSol=string(string(folder_with_elements,string("/","Res_L_",id)),".txt")
 Handle_PathSol=open(PathSol,"w")
 writedlm(Handle_PathSol,field)
 close(Handle_PathSol)
end
#mat = Matern(0.75,2.0,σ= 1.0)
