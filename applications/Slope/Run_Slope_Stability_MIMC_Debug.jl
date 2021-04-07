
#module Run
push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))
push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
push!(LOAD_PATH,"/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/packages/MultilevelEstimators/KoRDo/applications/SPDE")
using Distributed
using MultilevelEstimators
using Coupling_Slope_MIMC_DEBUG
using Pkg
using MATLAB
using DelimitedFiles
using Revise
using DigitalNets
using LatticeRules


#./julia /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/Run.jl 2>&1 | tee /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/10112018_Run_Log.txt
numberOfProcs=24

#addprocs(numberOfProcs)
#println("Procs added")
#println(numberOfProcs)

#@everywhere using Coupling
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
@everywhere push!(LOAD_PATH,"/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE")


@everywhere using Coupling_Slope_MIMC_DEBUG
@everywhere using MATLAB


#pathToInterm=joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE","data")

#pathToInterm=joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE","data")

#pathToInterm="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/packages/MultilevelEstimators/KoRDo/applications/SPDE/data"

pathToInterm=joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE","data")

folder = string(pathToInterm,"/Interm/Slope") # for report
folder_with_elements=string(pathToInterm,"/Mesh/Slope")
if(isdir(folder)==false)
mkdir(folder)
else
    rm(folder,recursive=true)
    mkdir(folder)
end



@everywhere MatlabSampler() = eval_string(string("Slope_Slover_JULIA_MATLAB(", myid(), ")"))



#--------------------
#MIMC
#------------------
nterms=400
pt=LatticeRule("/home/philippe/.julia/dev/Applications/applications/SPDE/lattice-32001-1024-1048576.3600.txt",nterms)

init_Beam_MC_L_Het_Single=Coupling_Slope_MIMC_DEBUG.init_Slope(SL(),true,false,false,true,true,startlevel=0,max_level=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=1,nb_of_warm_up_samples=2,nshifts=10,nterms=nterms,correlateOnlyDiffs=true,corr_len=1.0,point_generator=pt,adaptive_MIMC=true)
estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,1.75e-04)
