#module KLEigenvals


using   GaussianRandomFields,Revise, DelimitedFiles

p=2
nterms=586
nq=Int64(ceil(sqrt(3*nterms)))
#exp_field = GaussianRandomFields.Exponential(0.3,σ=1,p=p)

#cov = CovarianceFunction(2,exp_field)
exp_field = GaussianRandomFields.Matern(0.3,0.6,σ=1,p=p)
#exp_field = GaussianRandomFields.Exponential(1.0,σ=1,p=p)

cov = CovarianceFunction(2,exp_field)

vx = 0:0.1:1
vy= 0:0.1:1
grfs = GaussianRandomField(cov,KarhunenLoeve(nterms),vx,vy,quad=GaussLegendre(),nq=nq)
pathToInterm=joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE","data")
folder_with_elements=string(pathToInterm,"/Mesh")


println(grfs)
#println(grfs.data.eigenval.^2)
eigv=grfs.data.eigenval.^2

Zc_res=open(string(folder_with_elements,"/KL"),"w")
writedlm(Zc_res,eigv)
close(Zc_res)

#end
