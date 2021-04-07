module RF

using Dierckx, DelimitedFiles, GaussianRandomFields, Statistics, Random

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

level=2

smoothness=1.0
corr_len=3
p=2
exp_field = GaussianRandomFields.SquaredExponential(corr_len,σ=smoothness,p=p)
#exp_field = Matern(4,5.0)

println("P of covar equals")
println(p)
cov = CovarianceFunction(2,exp_field)


folder_with_elements="/home/philippe/Desktop/GMSH_Matlab_ReadIn/Mesh_Steep/Matlab"

if level>0

PathElements_Fine=string(string(folder_with_elements,string("/Elements_L_",level)),".txt")
PathNodes_Fine=string(string(folder_with_elements,string("/Nodes_L_",level)),".txt")



Handle_Nodes=open(PathNodes_Fine)
Node=readdlm(Handle_Nodes);
Node_fine=Node[:,2:3]

Handle_Elements=open(PathElements_Fine)
Element=readdlm(Handle_Elements,Int);
Element_fine=Element[:,5:7]



level_1=level-1
println(level_1)
PathElements_Coarse=string(string(folder_with_elements,string("/Elements_L_",level_1)),".txt")
PathNodes_Coarse=string(string(folder_with_elements,string("/Nodes_L_",level_1)),".txt")


Handle_Nodes=open(PathNodes_Coarse)
Node=readdlm(Handle_Nodes);
Node_coarse=Node[:,2:3]

Handle_Elements=open(PathElements_Coarse)
Element=readdlm(Handle_Elements,Int);
Element_coarse=Element[:,5:7]

end
Random.seed!(1234)
grf=GaussianRandomField(cov,KarhunenLoeve(101),Node_fine,Element_fine,quad=GaussLegendre())
Zf = GaussianRandomFields.sample(grf) # compute GRF
Zf = GaussianRandomFields.sample(grf)
Zf = GaussianRandomFields.sample(grf)

#Zf=expfield(Zf)
spline = Spline2D(Node_fine[:,1],Node_fine[:,2],Zf,w=ones(length(Node_fine[:,1])),kx=1,ky=1,s=10)
Zf=expfield(Zf)
#Zf=exp.(Zf)

field_fine=open(string(folder_with_elements,string("/field_L_",string(level,"_Original"))),"w")
writedlm(field_fine,Zf)
close(field_fine)
Random.seed!(1234)
grf_1=GaussianRandomField(cov,KarhunenLoeve(101),Node_coarse,Element_coarse,quad=GaussLegendre())
Zc = GaussianRandomFields.sample(grf_1) # compute GRF
Zc = GaussianRandomFields.sample(grf_1)
Zc = GaussianRandomFields.sample(grf_1)
#Zc=evaluate(spline,Node_coarse[:,1],Node_coarse[:,2])
#Zc=exp.(Zc)
println(length(Zc))
println(length(Node_coarse[:,1]))

field_coarse=open(string(folder_with_elements,string("/field_L_",string(level_1,"_interp"))),"w")
Zc=expfield(Zc)
writedlm(field_coarse,Zc)
close(field_coarse)


end
