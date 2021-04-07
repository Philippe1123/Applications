#__precompile__()

#module Run
push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","MultilevelEstimators","applications","SPDE")))
push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
push!(LOAD_PATH,"/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE")

#push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))

using Distributed
using MultilevelEstimators
using Coupling_Beam
using Pkg
using MATLAB
using DelimitedFiles
using Random
using LatticeRules
using DigitalNets
#using SharedArrays
#using DistributedArrays
#push!(LOAD_PATH,(joinpath("/","home","philippe","JuliaRuns","MultilevelEstimators","applications","SPDE")))
#push!(LOAD_PATH,Pkg.dir(joinpath("MultilevelEstimators","applications","SPDE")))
#push!(LOAD_PATH,(joinpath("/","home","philippe","JuliaRuns","MultilevelEstimators","src")))
#push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","v0.6","MultilevelEstimators","applications","SPDE")))

#./julia /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/Run.jl 2>&1 | tee /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/10112018_Run_Log.txt
numberOfProcs=22

addprocs(numberOfProcs)
println("Procs added")
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
@everywhere push!(LOAD_PATH,"/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE")

@everywhere using Coupling_Beam
@everywhere using MATLAB

#pathToInterm=joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE","data")

#pathToInterm=joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE","data")

pathToInterm="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE/data"

folder = string(pathToInterm,"/Interm9/Beam") # for report
folder_with_elements=string(pathToInterm,"/Mesh/Beam")

if(isdir(folder)==false)
mkdir(folder)
else
    rm(folder,recursive=true)
    mkdir(folder)
end



nterms=100
pt=DigitalNet64_2(nterms) # changed to sobol3
@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN_SOBOL3(", myid(), ")"))
init_Beam_mlqmc_L_Het_Single_HigherOrder=Coupling_Beam.init_Beam(SL(),true,false,false,false,true,true,nb_of_warm_up_samples=2,max_level=2,startlevel=2,nshifts=8,MatlabSampler,folder,folder_with_elements,true,nterms=nterms,NQoI=1,numberoftol=200,smoothness=6.0,ptgen=pt,corr_len=2.0)
estimator=init_Beam_mlqmc_L_Het_Single_HigherOrder
if(estimator.internals.sample_method_internals.generators[1][1] isa AbstractDigitalNets)
if(Int64.(DigitalNets.reversebits((estimator.internals.sample_method_internals.generators[1][1].digital_net.C)[1,1]))==1)
println("Used point set is: Sobol_CS")
elseif(Int64.(DigitalNets.reversebits((estimator.internals.sample_method_internals.generators[1][1].digital_net.C)[1,1]))==3)
println("Used point set is: Sobol_CS_2")
elseif(Int64.(DigitalNets.reversebits((estimator.internals.sample_method_internals.generators[1][1].digital_net.C)[1,1]))==7)
println("Used point set is: Sobol_CS_3")
end
elseif(estimator.internals.sample_method_internals.generators[1][1] isa AbstractLatticeRule)
println("Lattice")
#println(Int64.(estimator.internals.sample_method_internals.generators[1][1].lattice_rule.z))
end
history = run(estimator,2.0e-12)