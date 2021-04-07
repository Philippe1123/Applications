using Dierckx, DelimitedFiles, GaussianRandomFields, Statistics


PathElements_Fine="/home/philippe/Desktop/GMSH_Matlab_ReadIn/Elements_Fine.txt"
PathNodes_Fine="/home/philippe/Desktop/GMSH_Matlab_ReadIn/Points_Fine.txt"
PathRes_Fine="/home/philippe/Desktop/GMSH_Matlab_ReadIn/field_Fine.txt"

PathElements_Coarse="/home/philippe/Desktop/GMSH_Matlab_ReadIn/Elements_Coarse.txt"
PathNodes_Coarse="/home/philippe/Desktop/GMSH_Matlab_ReadIn/Points_Coarse.txt"
PathRes_Coarse="/home/philippe/Desktop/GMSH_Matlab_ReadIn/field_Coarse.txt"

PathRes_Fine_EigFun="/home/philippe/Desktop/GMSH_Matlab_ReadIn/Eig_Fine_KL.txt"
EigFun_Fine=open(PathRes_Fine_EigFun,"w")


Elements_Fine=open(PathElements_Fine)
Elements_Fine=readdlm(Elements_Fine,Int);

Nodes_Fine=open(PathNodes_Fine)
Nodes_Fine=readdlm(Nodes_Fine);

Nodes_Coarse=open(PathNodes_Coarse)
Nodes_Coarse=readdlm(Nodes_Coarse);

mat = Matern(1.2,2.0)
cov = CovarianceFunction(2,mat)
grf = GaussianRandomField(cov,KarhunenLoeve(100),Nodes_Fine,Elements_Fine)

println(length(view(grf.data.eigenfunc, :, 1)))
writedlm(EigFun_Fine,view(grf.data.eigenfunc, :, 1))
close(EigFun_Fine)

field_Fine=sample(grf)



println(length(Nodes_Fine[:,2]))

spline = Spline2D(Nodes_Fine[:,1],Nodes_Fine[:,2],field_Fine,kx=3,ky=3,s=10)
println("Spline")
field_Coarse=evaluate(spline,Nodes_Coarse[:,1],Nodes_Coarse[:,2])
println("coarse")

Samp_Fine=open(PathRes_Fine,"w")
Samp_Coarse=open(PathRes_Coarse,"w")

writedlm(Samp_Coarse,field_Coarse)
close(Samp_Coarse)

writedlm(Samp_Fine,field_Fine)
close(Samp_Fine)
