using GaussianRandomFields, Random, DelimitedFiles

λ = 0.3 # correlation length
ν = 1.0 # smoothness
n = 1400 # number of KL terms
ell = 6 # 'level'

# covariance function
C = CovarianceFunction(2, Matern(λ, ν))

ptsx = range(0, stop=1, length=201)
ptsy = range(0, stop=5, length=201)
ptsz = range(0, stop=1, length=201)


#@time grf_Spec = GaussianRandomField(C, Spectral(), pts, pts,pts)

#@time grf_cholesky = GaussianRandomField(C, Cholesky(), pts, pts,pts)

Random.seed!(1234)
@time grf_CE = GaussianRandomField(C, CirculantEmbedding(), ptsx, ptsy,minpadding=2)
#@time grf_CE = GaussianRandomField(C, KarhunenLoeve(n), ptsx, ptsy,ptsz,eigensolver=EigsSolver())


Random.seed!(1234)
sample_grf=GaussianRandomFields.sample(grf_CE)
println(typeof(sample_grf))
#sample_grf=GaussianRandomFields.sample(grf_CE,xi=vec(collect(1:randdim(grf_CE))))

stp=open("/home/philippe/Desktop/3d_Matlab/FieldInstances/3dfield.txt","w")
 writedlm(stp,sample_grf)
 close(stp)

 stp=open("/home/philippe/Desktop/3d_Matlab/FieldInstances/3dfield_pts_x.txt","w")
  writedlm(stp,grf_CE.pts[1])
  close(stp)

  stp=open("/home/philippe/Desktop/3d_Matlab/FieldInstances/3dfield_pts_y.txt","w")
   writedlm(stp,grf_CE.pts[2])
   close(stp)

   stp=open("/home/philippe/Desktop/3d_Matlab/FieldInstances/3dfield_pts_z.txt","w")
    writedlm(stp,grf_CE.pts[3])
    close(stp)
