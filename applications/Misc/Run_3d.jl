push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))
push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
push!(LOAD_PATH,"/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE")


using Distributed
using MultilevelEstimators
using Coupling_3d
using Pkg
using MATLAB
using DelimitedFiles
using Populate_B_3d

using GaussianRandomFields
using FFTW


#./julia /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/Run.jl 2>&1 | tee /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/10112018_Run_Log.txt
numberOfProcs=24

addprocs(numberOfProcs)
println("Procs added")
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
@everywhere push!(LOAD_PATH,"/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE")

@everywhere using Coupling_3d
@everywhere using MATLAB
@everywhere using GaussianRandomFields
@everywhere using FFTW

#pathToInterm=joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE","data")

#pathToInterm=joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE","data")

#pathToInterm="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/packages/MultilevelEstimators/KoRDo/applications/SPDE/data"


#folder = string("/home/philippe/Desktop/3d_Matlab","/Interm/3d") # for report
#folder_with_elements=string("/home/philippe/Desktop/3d_Matlab/GenerateMeshes_1_p_dev_beam/p_ref/Matlab","/Mesh/3d")




folder = string("/scratch/leuven/330/vsc33032","/Interm/3d")
folder_with_elements=string("/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE/data/Mesh/3d_beam")

if(isdir(folder)==false)
mkdir(folder)
else
    rm(folder,recursive=true)
    mkdir(folder)
end




@everywhere MatlabSampler() = eval_string(string("Solver_3d_MATLAB(", myid(), ")"))



#--------------------
#MLMC
#------------------
maximum_level=4
#Populate_B_3d.init(maximum_level,folder_with_elements,folder)
init_Beam_MC_L_Het_Single=Coupling_3d.init_3d(ML(),true,false,false,true,true,startlevel=0,MatlabSampler,folder,folder_with_elements,false,max_level=maximum_level,numberoftol=31,nb_of_warm_up_samples=2,nshifts=10,corr_len=0.3,smoothness=2.0)
estimator=init_Beam_MC_L_Het_Single
history = run(estimator,2.e-06)
