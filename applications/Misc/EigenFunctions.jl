using DelimitedFiles, GaussianRandomFields


PathElements_Fine="/home/philippe/Desktop/GMSH_Matlab_ReadIn/Elements_Fine.txt"
PathNodes_Fine="/home/philippe/Desktop/GMSH_Matlab_ReadIn/Points_Fine.txt"
PathRes_Fine="/home/philippe/Desktop/GMSH_Matlab_ReadIn/field_Fine.txt"



PathRes_Fine_EigFun_KL="/home/philippe/Desktop/GMSH_Matlab_ReadIn/Eig_Fine_KL.txt"
EigFun_Fine_KL=open(PathRes_Fine_EigFun_KL,"w")


Elements_Fine=open(PathElements_Fine)
Elements_Fine=readdlm(Elements_Fine,Int);

Nodes_Fine=open(PathNodes_Fine)
Nodes_Fine=readdlm(Nodes_Fine);



mat = Matern(0.75,2.0)
cov = CovarianceFunction(2,mat)
grf_KL = GaussianRandomField(cov,KarhunenLoeve(200),Nodes_Fine,Elements_Fine)

writedlm(EigFun_Fine_KL,view(grf_KL.data.eigenfunc, :, 1:10))
close(EigFun_Fine_KL)


PathRes_Fine_EigFun_Spec="/home/philippe/Desktop/GMSH_Matlab_ReadIn/Eig_Fine_Spec.txt"
EigFun_Fine_Spec=open(PathRes_Fine_EigFun_Spec,"w")
grf_Spec = GaussianRandomField(cov,Spectral(),Nodes_Fine,Elements_Fine,n=200)

writedlm(EigFun_Fine_Spec,view(grf_Spec.data.eigenfunc, :, 1:10))
close(EigFun_Fine_Spec)
