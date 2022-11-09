using ScatteredInterpolation
using MultilevelEstimators
using Revise
using DelimitedFiles
using Statistics
using GaussianRandomFields
#using FieldTransformation
using Random
using PyPlot

function main()
#	PathNodes=string("/home/philippe/.julia/dev/Applications/applications/Mesh/Slope/Meshes_paper_mdpi_disp/Original Runs/h_refinement/ElementsCenter_L_4.txt")
    PathNodes=string("/home/philippe/.julia/dev/Applications/applications/Mesh/Slope/Meshes_paper_mdpi_disp/Original Runs/h_refinement/Nodes_L_0.txt")

    elem = readdlm("/home/philippe/.julia/dev/Applications/applications/Mesh/Slope/Meshes_paper_mdpi_disp/Original Runs/h_refinement/Elements_L_0.txt")
    elem=elem[:,5:end]
    nodes = readdlm(PathNodes)
    nodes=nodes[:,2:3]

    matern=Matern(1.,0.6,σ=1.0,p=2)
    cov=CovarianceFunction(2,matern)
    grfs = GaussianRandomField(cov,KarhunenLoeve(210),nodes,Int64.(elem),quad=GaussLegendre())
    a=randn(210,1)
    f = expfield(sample(grfs,xi=a))


    figure()
    t=scatter3D(nodes[:,1],nodes[:,2],f)
    gca()[:view_init](30,-60)
     xlabel("x position")
    ylabel("y position")
    zlabel("Value of the random field")
"""
    figure()
    scatter3D(nodes[:,1],nodes[:,2],f)
    gca()[:view_init](90,-90)
    xlabel("x position")
    ylabel("y position")
"""
	#plt3d= Plots.plot(nodes[:,1],nodes[:,2],expfield(sample(grfs)),
	  #          seriestype=:scatter, markersize = 7)


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

function expfield(Z::Vector{T})where {T<:Real}

Mean=log(6000)
σ=0.05

Z=Z.*σ
Z=Z.+Mean
Z=exp.(Z)

return Z
end


main()