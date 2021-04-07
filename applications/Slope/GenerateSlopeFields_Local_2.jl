push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","MultilevelEstimators","applications","SPDE")))
push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))

using ScatteredInterpolation
using MultilevelEstimators
using Revise
using DelimitedFiles
using Statistics
using GaussianRandomFields
using FieldTransformation
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

Mean=log(8000)
σ=0.05

Z=Z.*σ
Z=Z.+Mean
Z=exp.(Z)

return Z
end



pathToInterm=joinpath("/","home","philippe","Desktop","Slope_Matlab")
folder = string(pathToInterm,"/Field_Search") # for report
folder_with_elements=string(pathToInterm,"/Field_Search/Mesh")
#println(folder)


#println(folder_with_elements)

p=2
#exp_field = GaussianRandomFields.Exponential(1.0,σ=1,p=2)
#exp_field = GaussianRandomFields.Exponential(1.0,σ=1,p=2)
exp_field = GaussianRandomFields.Matern(0.3,1.0,σ=1.,p=2)
#exp_field = GaussianRandomFields.SquaredExponential(5.0,σ=1,p=2)
println(exp_field)
Nterms=400
lv=7





#exp_field = GaussianRandomFields.Exponential(0.3,σ=1,p=p)
println("P of covar equals")
println(p)

cov = CovarianceFunction(2,exp_field)


# Note fields are generated as vy x vy thus 160 x 40 example the solve transposes this



println(folder_with_elements)


field_res=string(folder)
PathNodes=string(string(folder_with_elements,string("/","GaussPoints_L_7_M7_S7")),".txt")
Handle_Nodes=open(PathNodes)
Node=readdlm(Handle_Nodes);
Node=Node[:,2:3]
Nodes_Fine=Node

PathElement=string(string(folder_with_elements,string("/","Elements_L","_7")),".txt")

println(PathElement)

Handle_Elements=open(PathElement)
Element=readdlm(Handle_Elements,Int);
Element=Element[:,5:7]

Random.seed!(137899)
grfs = GaussianRandomField(cov,KarhunenLoeve(Nterms),Nodes_Fine,Element,quad=GaussLegendre())
Zf = GaussianRandomFields.sample(grfs) # compute GRF
#Zf=FieldTransformation.Transform_L(Zf)
#Zf=expfield(Zf)

println(minimum(Zf))

Zf_res=open(string(field_res,"/Fine"),"w")
#Zf=grfs[index].data.eigenval

writedlm(Zf_res,Zf)
close(Zf_res)

if(lv>0)
    field_res=string(folder)
    PathNodesC=string(string(folder_with_elements,string("/","GaussPoints_L_6_M7_S6")),".txt")
    Handle_NodesC=open(PathNodesC)
    NodeC=readdlm(Handle_NodesC);
    NodeC=NodeC[:,2:3]
    Nodes_C=NodeC

    PathElementC=string(string(folder_with_elements,string("/","Elements_L","_6")),".txt")

    println(PathElementC)

    Handle_ElementsC=open(PathElementC)
    ElementC=readdlm(Handle_ElementsC,Int);
    ElementC=ElementC[:,5:7]

    Random.seed!(137899)
    grfsC = GaussianRandomField(cov,KarhunenLoeve(Nterms),Nodes_C,ElementC,quad=GaussLegendre())
    Zc = GaussianRandomFields.sample(grfsC) # cSompute GRF
    #Zf=FieldTransformation.Transform_L(Zf)
#    Zf=expfield(Zf)

    println(minimum(Zc))

    Zc_res=open(string(field_res,"/Coarse"),"w")
    #Zf=grfs[index].data.eigenval

    writedlm(Zc_res,Zc)
    close(Zc_res)
end
data_type=open(string(field_res,"/lvl"),"w")
writedlm(data_type,lv)
close(data_type)
