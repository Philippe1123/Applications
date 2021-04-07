using GaussianRandomFields, Plots, LinearAlgebra
gr()
λ = 1 # correlation length
ν = .75 # smoothness
n = 150 # number of KL terms
ell = 6 # 'level'

# covariance function
C = CovarianceFunction(2, Matern(λ, ν))

# points
pts = range(0, stop=1, length=2^ell)

# compute gaussian random fields with eig / eigs
#grf = GaussianRandomField(C, KarhunenLoeve(n), pts, pts, eigensolver=EigenSolver())
grf = GaussianRandomField(C, KarhunenLoeve(n), pts, pts, eigensolver=EigsSolver())
println(grf)
# plot first 9 eigenfunctions
p = map(i -> plot_eigenfunction(grf, i), 1:9)
plot(p..., layout=(3,3))
