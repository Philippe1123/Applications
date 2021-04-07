
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
folder = string(pathToInterm,"/Mesh_Quad") # for report
folder_with_elements="/home/philippe/Desktop/Slope_Matlab/Mesh_Quad/p_ref/Matlab"
#folder_with_elements="/home/philippe/Desktop/RunFiles/debug/MLMC_2/h_refinement"

#println(folder)


#println(folder_with_elements)

p=2

exp_field = GaussianRandomFields.Matern(1.5,1.0,σ=1.,p=2)
#exp_field = GaussianRandomFields.SquaredExponential(5.0,σ=1,p=2)
println(exp_field)
Nterms=400
lv=1





#exp_field = GaussianRandomFields.Exponential(0.3,σ=1,p=p)
println("P of covar equals")
println(p)

cov = CovarianceFunction(2,exp_field)


# Note fields are generated as vy x vy thus 160 x 40 example the solve transposes this



println(folder_with_elements)


field_res=string(folder)
PathNodes=string(string(folder_with_elements,string("/","GaussPoints_L_0")),".txt")
Handle_Nodes=open(PathNodes)
Node=readdlm(Handle_Nodes);
Node=Node[:,2:3]
Nodes_Fine=Node

PathElement=string(string(folder_with_elements,string("/","Elements_L","_0")),".txt")

println(PathElement)

Handle_Elements=open(PathElement)
Element=readdlm(Handle_Elements,Int);
Element=Element[:,5:7]

Random.seed!(12)
grfs = GaussianRandomField(cov,KarhunenLoeve(Nterms),Nodes_Fine,Element,quad=GaussLegendre())
Zf = GaussianRandomFields.sample(grfs) # compute GRF


println(expfield(Zf))
Zf=expfield(Zf)
Zf_res=open("/home/philippe/Desktop/Slope_Matlab/Mesh_Quad/Field_0","w")

writedlm(Zf_res,Zf)
close(Zf_res)
