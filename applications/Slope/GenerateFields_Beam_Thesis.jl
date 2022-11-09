using ScatteredInterpolation
using MultilevelEstimators
using Revise
using DelimitedFiles
using Statistics
using GaussianRandomFields
#using FieldTransformation
using Random
using PyPlot
using Distributions
using Interpolations

include("FieldTransformation.jl")

function main()
	PathNodes=string("/home/philippe/.julia/dev/Applications/applications/Mesh/Beam/PaperMDPI/h_refinement/Nodes_L_4.txt")
    elem = readdlm("/home/philippe/.julia/dev/Applications/applications/Mesh/Beam/PaperMDPI/h_refinement/Elements_L_4.txt")
    nodes = readdlm(PathNodes)
    nodes=nodes[:,2:3]
    elem=Int64.(elem[:,5:7])

    Center=compute_centers(nodes,elem)
    matern=Matern(.3,0.6,σ=1.0,p=2)
    cov=CovarianceFunction(2,matern)
    grfs = GaussianRandomField(cov,KarhunenLoeve(587),nodes,elem,quad=GaussLegendre())
    a=randn(587,1)
    f = FieldTransformation.Transform_L(GaussianRandomFields.sample(grfs,xi=a))


    figure()
    t=scatter3D(nodes[:,1],nodes[:,2],f)
    gca()[:view_init](30,-60)
     xlabel("x position")
    ylabel("y position")
    zlabel("Value of the random field")

"""
    grf2 = GaussianRandomField(cov,KarhunenLoeve(587),nodes,elem,quad=GaussLegendre())
    field = FieldTransformation.Transform_L(GaussianRandomFields.sample(grf2,xi=a))

    itp = ScatteredInterpolation.interpolate(NearestNeighbor(), nodes', field);
    field = evaluate(itp, Center')


    figure()
    t=scatter3D(Center[:,1],Center[:,2],field)
    gca()[:view_init](30,-60)
     xlabel("x position")
    ylabel("y position")
    zlabel("Value of the random field")
"""

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

"""
function Transform_L(Points::Array{Float64,1})
    Cov=0.25
    
    beta=30E3*Cov^2/(1-Cov^2);
    alpha=30E3/beta;
    GamDist=Gamma(alpha,beta)
    #@show Cov
    Evaluation_points_Gam=(0:0.01:10)*10.0E4;
    CDF_Gamma=cdf.(GamDist,Evaluation_points_Gam);
    NormDist=Distributions.Normal(0.0,1.0)
    Evaluation_points_Norm=-4.0:5.0^-4:4.0;
    CDF_Norm=cdf.(NormDist,Evaluation_points_Norm)
    cdf_mat=interp1(Evaluation_points_Norm,CDF_Norm,Points)
    GammaField=interp1(CDF_Gamma,Evaluation_points_Gam,cdf_mat)
    
    #@show CDF_Gamma
    #figure()
    #plot(Evaluation_points_Gam,CDF_Gamma)
    
    #figure()
    #plot(Evaluation_points_Norm,CDF_Norm)
    
    return GammaField
    end

    function interp1(xpt, ypt, x)
        y = zeros(size(x,1),size(x,2))
        idx = trues(size(x,1),size(x,2))


        intf = Interpolations.extrapolate(Interpolations.interpolate((xpt,), ypt, Interpolations.Gridded(Interpolations.Linear())),Interpolations.Flat())

         y= intf(x)



    return y
end
"""

main()