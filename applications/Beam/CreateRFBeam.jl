module CreateRFBeam

using  DelimitedFiles, GaussianRandomFields, Statistics, MultilevelEstimators,Plots,LinearAlgebra,Random

get_max_index_set(index_set, Int64) = get_index_set(index_set, Int64)

#plotly()
pyplot()
nterms=150
he_start=250/4000
Lx = 2.5;
Ly = 0.25;
he = he_start;

nelx=Lx/he
nely=Ly/he
max_level=5

indices = get_max_index_set(ML(), max_level)



distributions = [MultilevelEstimators.Normal() for i in 1:nterms]
corr_len=1.0
smoothness=1.0
p=2
   exp_field = GaussianRandomFields.Exponential(corr_len,σ=smoothness,p=p)
   println("P of covar equals")
   println(p)

   cov = CovarianceFunction(2,exp_field)
   let
   grfs=Dict()

   i=0.0
   coarse_dof=1
   for index in indices
       j=i
       m = coarse_dof*2^i
       n = coarse_dof*2^j

#       vx = 0:(he/m)/Lx:0.99999
#       vy= 0:(he/n)/Ly:0.99999
       vx=1/80/(2^i):1/40/(2^i):1-1/80/(2^i)
       vy=1/8/(2^i):1/4/(2^i):1-1/8/(2^i)

   #    println(vx)
   #    println(vy)
       Random.seed!(1234)
   #    println(nterms)
      grfs[index] = GaussianRandomField(cov,KarhunenLoeve(nterms),vx,vy,eigensolver=EigenSolver())
      Random.seed!(1234)

      i=i+1
   end

println(grfs)

#figure()
Random.seed!(1234)
Field0=GaussianRandomFields.sample(grfs[Index(0)],xi=randn(nterms,1))
Random.seed!(1234)
Field1=GaussianRandomFields.sample(grfs[Index(1)],xi=randn(nterms,1))
Random.seed!(1234)
Field2=GaussianRandomFields.sample(grfs[Index(2)],xi=randn(nterms,1))
Random.seed!(1234)
Field3=GaussianRandomFields.sample(grfs[Index(3)],xi=randn(nterms,1))
Random.seed!(1234)
Field4=GaussianRandomFields.sample(grfs[Index(4)],xi=randn(nterms,1))
Random.seed!(1234)
Field5=GaussianRandomFields.sample(grfs[Index(5)],xi=randn(nterms,1))

println("Sampled")

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,1],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,1],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,1],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,1],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,1],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,1],1280,128),fill=true)
p1=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))

display(p1)


c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,2],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,2],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,2],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,2],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,2],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,2],1280,128),fill=true)
p12=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))

display(p12)


c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,3],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,3],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,3],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,3],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,3],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,3],1280,128),fill=true)
p13=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p13)

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,4],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,4],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,4],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,4],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,4],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,4],1280,128),fill=true)
p14=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p14)

c0=plot(grfs[Index(0)].data.eigenfunc[:,4])
c1=plot(grfs[Index(1)].data.eigenfunc[:,4])
c2=plot(grfs[Index(2)].data.eigenfunc[:,4])
c3=plot(grfs[Index(3)].data.eigenfunc[:,4])
c4=plot(grfs[Index(4)].data.eigenfunc[:,4])
c5=plot(grfs[Index(5)].data.eigenfunc[:,4])
p144=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p144)

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,5],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,5],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,5],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,5],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,5],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,5],1280,128),fill=true)
p15=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p15)

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,6],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,6],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,6],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,6],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,6],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,6],1280,128),fill=true)
p16=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p16)

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,7],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,7],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,7],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,7],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,7],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,7],1280,128),fill=true)
p17=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p17)

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,8],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,8],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,8],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,8],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,8],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,8],1280,128),fill=true)
p18=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p18)

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,9],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,9],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,9],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,9],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,9],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,9],1280,128),fill=true)
p19=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p19)

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,10],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,10],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,10],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,10],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,10],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,10],1280,128),fill=true)
p20=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p20)

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,11],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,11],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,11],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,11],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,11],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,11],1280,128),fill=true)
p21=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p21)

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,12],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,12],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,12],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,12],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,12],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,12],1280,128),fill=true)
p22=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p22)

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,12],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,12],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,12],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,12],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,12],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,12],1280,128),fill=true)
p22=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p22)

c0=contour(reshape(grfs[Index(0)].data.eigenfunc[:,13],40,4),fill=true)
c1=contour(reshape(grfs[Index(1)].data.eigenfunc[:,13],80,8),fill=true)
c2=contour(reshape(grfs[Index(2)].data.eigenfunc[:,13],160,16),fill=true)
c3=contour(reshape(grfs[Index(3)].data.eigenfunc[:,13],320,32),fill=true)
c4=contour(reshape(grfs[Index(4)].data.eigenfunc[:,13],640,64),fill=true)
c5=contour(reshape(grfs[Index(5)].data.eigenfunc[:,13],1280,128),fill=true)
p23=plot(c0,c1,c2,c3,c4,c5,reuse=false,layout=(2,3))
display(p23)


#display(plot(grfs[Index(0)].pts,Field0))
f0=(contour(Field0,reuse = false,fill=true))
f1=(contour(Field1,reuse = false,fill=true))
f2=(contour(Field2,reuse = false,fill=true))
f3=(contour(Field3,reuse = false,fill=true))
f4=(contour(Field4,reuse = false,fill=true))
f5=(contour(Field5,reuse = false,fill=true))
p2=plot(f0,f1,f2,f3,f4,f5,reuse=false,layout=(2,3))
display(p2)




Field1_1 = Array(view(Field1, 2:2:size(Field1, 1), 2:2:size(Field1, 2)))
f0=contour(Field0,fill=true)
f1=contour(Field1_1,fill=true)
f2=contour(Field2,fill=true)

p3=plot(f0,f1,f2,reuse=false,layout=(1,3))
display(p3)









#display(plot(grfs[Index(0)].pts,Field0))
#gui(plot(grfs[Index(0)].pts,Field0))
#view(p1,p2,layout=(1,2))
#gui(contour(grfs[Index(1)].pts,Field1))
#s0=plot(Field0,reuse = false,st=:surface)
#s1=plot(Field1,reuse = false,st=:surface)
#s2=plot(Field2,reuse = false,st=:surface)
#s3=plot(Field3,reuse = false,st=:surface)
#s4=plot(Field4,reuse = false,st=:surface)
#s5=plot(Field5,reuse = false,st=:surface)
#display(plot(s0,s1,s2,s3,s4,s5,layput=(2,3)))

#figure()
#plot(Statistics.cov(Field0))
#plot(Statistics.cor(Field0,dims=2))

#plot(Statistics.cor(Field0,dims=1))
#figure()
#plot(grfs[Index(0)].data.eigenval)
#figure()
#plot(view(grfs[Index(0)].data.eigenfunc, :, 1:nterms))

#Random.seed!(1234)
#Field1=GaussianRandomFields.sample(grfs[Index(1)],xi=randn(nterms,1))
#plot(Field1)
#Field1_1 = Array(view(Field1, 2:2:size(Field1, 1), 2:2:size(Field1, 2)))
#figure()
#surf(Field1_1)
#figure()
#plot(Statistics.cov(Field1))
#figure()
#plot(grfs[Index(1)].data.eigenval)
#figure()
#plot(view(grfs[Index(1)].data.eigenfunc, :, 1:nterms))

#Random.seed!(1234)
#Field2=GaussianRandomFields.sample(grfs[Index(2)],xi=randn(nterms,1))
#plot(Field2)
#Field2_1 = Array(view(Field2, 2:2:size(Field2, 1), 2:2:size(Field2, 2)))
#figure()
#surf(Field2_1)
#figure()
#plot(Statistics.cov(Field2))
#figure()
#plot(grfs[Index(2)].data.eigenval)
#figure()
#plot(view(grfs[Index(2)].data.eigenfunc, :, 1:nterms))

#Random.seed!(1234)
#Field3=GaussianRandomFields.sample(grfs[Index(3)],xi=randn(nterms,1))
#plot(Field3)
#Field3_1 = Array(view(Field3, 2:2:size(Field3, 1), 2:2:size(Field3, 2)))
#figure()
#surf(Field3_1)
#figure()
#plot(Statistics.cov(Field3))

#figure()
#plot(grfs[Index(3)].data.eigenval)
#figure()
#plot(view(grfs[Index(3)].data.eigenfunc, :, 1:nterms))



#Random.seed!(1234)
#Field4=GaussianRandomFields.sample(grfs[Index(4)],xi=randn(nterms,1))
#plot(Field4)
#Field4_1 = Array(view(Field4, 2:2:size(Field4, 1), 2:2:size(Field4, 2)))
#figure()
#surf(Field4_1)
#println(Statistics.cor(Field4,dims=2))
#plot(Statistics.cor(Field4,dims=1))
#plot(Statistics.cor(Field4,dims=2))
#figure()
#plot(view(grfs[Index(4)].data.eigenfunc, :, 1:nterms))

#figure()
#surf(grfs[Index(4)].data.eigenfunc)

#Field3_2 = Array(view(Field4, 2:2:size(Field4, 1), 2:2:size(Field4, 2)))
#figure()
#surf(Field3_2)
#figure()
#plot(Statistics.cov(Field3_2))

#Random.seed!(1234)
#Field5=GaussianRandomFields.sample(grfs[Index(5)],xi=randn(nterms,1))
#plot(Field5)
#plot(Statistics.cor(Field5,dims=1))
#plot(Statistics.cor(Field5,dims=2))
#println(Statistics.cov((Field)))
#println(Statistics.cor((Field)))
#figure()
#surf(Statistics.cor((Field)))

end
end
