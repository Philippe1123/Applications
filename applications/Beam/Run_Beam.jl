#__precompile__()

#module Run
push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))
push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
push!(LOAD_PATH,"/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/packages/MultilevelEstimators/KoRDo/applications/SPDE")

#push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))

using Distributed
using MultilevelEstimators
using Coupling_Beam
using Pkg
using MATLAB
using DelimitedFiles

#using SharedArrays
#using DistributedArrays
#push!(LOAD_PATH,(joinpath("/","home","philippe","JuliaRuns","MultilevelEstimators","applications","SPDE")))
#push!(LOAD_PATH,Pkg.dir(joinpath("MultilevelEstimators","applications","SPDE")))
#push!(LOAD_PATH,(joinpath("/","home","philippe","JuliaRuns","MultilevelEstimators","src")))
#push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","v0.6","MultilevelEstimators","applications","SPDE")))

#./julia /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/Run.jl 2>&1 | tee /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/10112018_Run_Log.txt
numberOfProcs=24

addprocs(numberOfProcs)
println("Procs added")
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
@everywhere push!(LOAD_PATH,"/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/packages/MultilevelEstimators/KoRDo/applications/SPDE")

@everywhere using Coupling_Beam
@everywhere using MATLAB

#pathToInterm=joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE","data")

#pathToInterm=joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE","data")

pathToInterm="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/packages/MultilevelEstimators/KoRDo/applications/SPDE/data"

folder = string(pathToInterm,"/Interm/Beam") # for report
folder_with_elements=string(pathToInterm,"/Mesh/Beam")

if(isdir(folder)==false)
mkdir(folder)
else
    rm(folder,recursive=true)
    mkdir(folder)
end



##-------------------------------------------------------------
## NON LINEAR
##----------------------------------------------------
#-------------
# Higer order elements
#------------
#In test phase on 02/04/2019

#Heterogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlmc_NL_Het_Single_HigherOrder= Coupling_Beam.init_Beam(ML(),false,false,false,true,true,true,250/4000,nb_of_warm_up_samples=40,max_level=3,MatlabSampler,folder,nterms=101)
#estimator=init_Beam_mlmc_NL_Het_Single_HigherOrder
#history = run(estimator,2.5e-6)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_NL_Het_Single_HigherOrder=Coupling_Beam.init_Beam(ML(),true,false,false,true,true,true,250/4000,nb_of_warm_up_samples=2,max_level=3,nshifts=10,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_mlqmc_NL_Het_Single_HigherOrder
#history = run(estimator,2.5e-6)

#Homogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlmc_NL_Hom_Single_HigherOrder= Coupling_Beam.init_Beam(ML(),false,false,false,true,false,true,250/4000,nb_of_warm_up_samples=40,max_level=3,MatlabSampler,folder,nterms=101)
#estimator=init_Beam_mlmc_NL_Hom_Single_HigherOrder
#history = run(estimator,2.5e-6)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_NL_Hom_Single_HigherOrder=Coupling_Beam.init_Beam(ML(),true,false,false,true,false,true,250/4000,nb_of_warm_up_samples=2,max_level=3,nshifts=10,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_mlqmc_NL_Hom_Single_HigherOrder
#history = run(estimator,2.5e-6)

##-------------
#-------------
# Refinement
#------------
#finised and used for Korea Paper
#In test fase on 02/04/2019
#Heterogeneous


#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_NL_Het_Single= Coupling_Beam.init_Beam(ML(),true,false,false,true,true,false,250/4000,nb_of_warm_up_samples=2,max_level=3,nshifts=10,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_mlqmc_NL_Het_Single
#history = run(estimator,2.5e-6)

#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_NL_Het_Single= Coupling_Beam.init_Beam(ML(),false,false,false,true,true,false,250/4000,max_level=3,MatlabSampler,folder,nterms=101)
#estimator=init_Beam_mlmc_NL_Het_Single
#history = run(estimator,2.5e-6)



#Homogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_NL_Hom_Single= Coupling_Beam.init_Beam(ML(),false,false,false,true,false,false,250/4000,max_level=3,MatlabSampler,folder,nterms=101)
#estimator=init_Beam_mlmc_NL_Hom_Single
#history = run(estimator,2.5e-6)

#In test fase on 02/04/2019
#mlqmc`
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_NL_Hom_Single= Coupling_Beam.init_Beam(ML(),true,false,false,true,false,false,250/4000,nb_of_warm_up_samples=2,max_level=3,nshifts=10,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_mlqmc_NL_Hom_Single
#history = run(estimator,2.5e-6)


##-------------------------------------------------------------
##-------------------------------------------------------------
##  LINEAR Testing on 06/04/2019
##----------------------------------------------------
#-------------
# Higer order elements
#------------
#Heterogeneous
##mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Het_Single_HigherOrder= Coupling_Beam.init_Beam(ML(),false,false,false,false,true,true,250/4000,nb_of_warm_up_samples=40,max_level=4,MatlabSampler,folder,nterms=101)
#estimator=init_Beam_mlmc_L_Het_Single_HigherOrder
#history = run(estimator,8.0e-5)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_L_Het_Single_HigherOrder=Coupling_Beam.init_Beam(ML(),true,false,false,false,true,true,250/4000,nb_of_warm_up_samples=2,max_level=4,nshifts=10,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_mlqmc_L_Het_Single_HigherOrder
#history = run(estimator,8.0e-5)

#Homogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Hom_Single_HigherOrder= Coupling_Beam.init_Beam(ML(),false,false,false,false,false,true,250/4000,nb_of_warm_up_samples=40,max_level=5,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_mlmc_L_Hom_Single_HigherOrder
#history = run(estimator,3.861e-05)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_L_Hom_Single_HigherOrder=Coupling_Beam.init_Beam(ML(),true,false,false,false,false,true,250/4000,nb_of_warm_up_samples=2,max_level=5,nshifts=10,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_mlqmc_L_Hom_Single_HigherOrder
#history = run(estimator,3.861e-05)



#Homogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Hom_Single= Coupling_Beam.init_Beam(ML(),false,false,false,false,false,false,250/4000,max_level=5,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_mlmc_L_Hom_Single
#history = run(estimator,3.861e-05)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_L_Hom_Single= Coupling_Beam.init_Beam(ML(),true,false,false,false,false,false,250/4000,nb_of_warm_up_samples=2,max_level=5,nshifts=10,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_mlqmc_L_Hom_Single
#history = run(estimator,3.861e-05)

#############################################################
################ MONTE CARLO
######################################################
#LINEAR

#HETEROGENEOUS

#HIGH
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Het_Single_high=Coupling_Beam.init_Beam(SL(),false,false,false,false,true,true,250/16000,startlevel=3,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_MC_L_Het_Single_high
#history = run(estimator,3.861e-05)
#REF WORKING
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Het_Single=Coupling_Beam.init_Beam(SL(),false,false,false,false,true,false,startlevel=3,numberoftol=10,MatlabSampler,folder,folder_with_elements,false,nterms=101)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,3.861e-05)



############
#HOMOGENEOUS

#HIGH
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Hom_Single_high=Coupling_Beam.init_Beam(SL(),false,false,false,false,false,true,250/4000,startlevel=3,numberoftol=10,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_MC_L_Hom_Single_high
#history = run(estimator,3.861e-05)
#REF
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Hom_Single=Coupling_Beam.init_Beam(SL(),false,false,false,false,false,false,250/4000,startlevel=3,numberoftol=6,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_MC_L_Hom_Single
#history = run(estimator,0.000110274021000000)
#######

#NONLINEAR

#HETEROGENEOUS

#HIGH
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_MC_NL_Het_Single_high=Coupling_Beam.init_Beam(SL(),false,false,false,true,true,true,250/4000,startlevel=3,numberoftol=8,MatlabSampler,folder,nterms=101)
#estimator=init_Beam_MC_NL_Het_Single_high
#history = run(estimator,4.22500000000000e-06)
#REF
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_NL_Het_Single=Coupling_Beam.init_Beam(SL(),false,false,false,true,true,false,250/4000,startlevel=3,numberoftol=7,MatlabSampler,folder,false,nterms=101)
#estimator=init_Beam_MC_NL_Het_Single
#history = run(estimator,5.49250000000000e-06)
############
#HOMOGENEOUS

#HIGH
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_MC_NL_Hom_Single_high=Coupling_Beam.init_Beam(SL(),false,false,false,true,false,true,250/4000,startlevel=3,numberoftol=6,MatlabSampler,folder,nterms=101)
#estimator=init_Beam_MC_NL_Hom_Single_high
#history = run(estimator,7.14025000000000e-06)
#REF
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_NL_Hom_Single=Coupling_Beam.init_Beam(SL(),false,false,false,true,false,false,250/4000,startlevel=3,numberoftol=3,MatlabSampler,folder,nterms=101)
#estimator=init_Beam_MC_NL_Hom_Single
#history = run(estimator,1.56871292500000e-05)
#######


####### MC

#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN_GP(", myid(), ")"))
#init_Beam_MC_L_Het_Single_high=Coupling_Beam.init_Beam(SL(),false,false,false,false,true,true,250/4000,startlevel=3,MatlabSampler,folder,true,nterms=101)
#estimator=init_Beam_MC_L_Het_Single_high
#history = run(estimator,3.861e-05)

#-------------------------------------------------------------
## NON LINEAR
##----------------------------------------------------
#-------------
# Higer order elements Gausspoints
#------------
#In test phase on

#Heterogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN_GP(", myid(), ")"))
#init_Beam_mlmc_NL_Het_Single_HigherOrder= Coupling_Beam.init_Beam(ML(),false,false,false,true,true,true,250/4000,nb_of_warm_up_samples=40,max_level=3,MatlabSampler,folder,true,nterms=101)
#estimator=init_Beam_mlmc_NL_Het_Single_HigherOrder
#history = run(estimator,2.5e-6)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN_GP(", myid(), ")"))
#init_Beam_mlqmc_NL_Het_Single_HigherOrder=Coupling_Beam.init_Beam(ML(),true,false,false,true,true,true,250/4000,nb_of_warm_up_samples=2,max_level=3,nshifts=10,MatlabSampler,folder,true,nterms=101)
#estimator=init_Beam_mlqmc_NL_Het_Single_HigherOrder
#history = run(estimator,2.5e-6)

####MC

#HIGH
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN_GP(", myid(), ")"))
#init_Beam_MC_NL_Het_Single_high=Coupling_Beam.init_Beam(SL(),false,false,false,true,true,true,250/4000,startlevel=3,numberoftol=8,MatlabSampler,folder,true,nterms=101)
#estimator=init_Beam_MC_NL_Het_Single_high
#history = run(estimator,4.22500000000000e-06)

#end








########################################################################## RUN DONE
##################3
######## LINEAR GAUSS POINTS
#-------------
# Higer order elements
#------------
#Heterogeneous
##mlmc Working
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Het_Single_HigherOrder= Coupling_Beam.init_Beam(ML(),false,false,false,false,true,true,nb_of_warm_up_samples=40,max_level=3,MatlabSampler,folder,folder_with_elements,true,nterms=78,NQoI=1)
#estimator=init_Beam_mlmc_L_Het_Single_HigherOrder
#history = run(estimator,2.0e-05)

#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Het_Single_HigherOrder= Coupling_Beam.init_Beam(ML(),false,false,false,false,true,true,nb_of_warm_up_samples=40,max_level=3,MatlabSampler,folder,folder_with_elements,true,nterms=586,NQoI=1)
#estimator=init_Beam_mlmc_L_Het_Single_HigherOrder
#history = run(estimator,2.0e-05)

#mlqm working
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_L_Het_Single_HigherOrder=Coupling_Beam.init_Beam(ML(),true,false,false,false,true,true,nb_of_warm_up_samples=2,max_level=3,nshifts=10,MatlabSampler,folder,folder_with_elements,true,nterms=78,NQoI=1)
#estimator=init_Beam_mlqmc_L_Het_Single_HigherOrder
#history = run(estimator,2.0e-05)

#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_L_Het_Single_HigherOrder=Coupling_Beam.init_Beam(ML(),true,false,false,false,true,true,nb_of_warm_up_samples=2,max_level=5,nshifts=10,MatlabSampler,folder,folder_with_elements,true,nterms=586,NQoI=1,numberoftol=30)
#estimator=init_Beam_mlqmc_L_Het_Single_HigherOrder
#history = run(estimator,2.0e-07)






#-------------
# Refinement
#------------
#Heterogeneous
#mlmc WORKING
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Het_Single= Coupling_Beam.init_Beam(ML(),false,false,false,false,true,false,max_level=5,MatlabSampler,folder,folder_with_elements,false,nterms=78,numberoftol=10,NQoI=1)
#estimator=init_Beam_mlmc_L_Het_Single
#history = run(estimator,2.0e-05)


#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Het_Single= Coupling_Beam.init_Beam(ML(),false,false,false,false,true,false,max_level=5,MatlabSampler,folder,folder_with_elements,false,nterms=586,numberoftol=10,NQoI=1)
#estimator=init_Beam_mlmc_L_Het_Single
#history = run(estimator,2.0e-05)


#mlqmc WORKING
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_L_Het_Single= Coupling_Beam.init_Beam(ML(),true,false,false,false,true,false,nb_of_warm_up_samples=2,max_level=5,nshifts=10,MatlabSampler,folder,folder_with_elements,false,nterms=78,numberoftol=10,NQoI=1)
#estimator=init_Beam_mlqmc_L_Het_Single
#history = run(estimator,2.0e-05)

#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_L_Het_Single= Coupling_Beam.init_Beam(ML(),true,false,false,false,true,false,nb_of_warm_up_samples=2,max_level=5,nshifts=10,MatlabSampler,folder,folder_with_elements,false,nterms=586,numberoftol=10,NQoI=1)
#estimator=init_Beam_mlqmc_L_Het_Single
#history = run(estimator,2.0e-05)




#-------------
# Refinement
#------------
#Heterogeneous
#mlmc WORKING
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Het_Single= Coupling_Beam.init_Beam(ML(),false,false,false,false,true,false,max_level=5,MatlabSampler,folder,folder_with_elements,false,nterms=78,numberoftol=9,NQoI=641)
#estimator=init_Beam_mlmc_L_Het_Single
#history = run(estimator,2.6e-05)
#S= MultilevelEstimators.samples_diff(estimator)
#name = estimator[:name]
#sf_dir= string(name[1:end-5])
#isdir(joinpath(pathToInterm,sf_dir,"Samples")) || mkdir(joinpath(pathToInterm,sf_dir,"Samples"))
#for idx in CartesianIndices(S)
#    idx_dir = joinpath(pathToInterm,sf_dir,"Samples", join(idx.I,"_"))
#    isdir(idx_dir) || mkdir(idx_dir)
#    for k in keys(S[idx])
#        writedlm(joinpath(idx_dir,string("samples_level_",k[1]-1,".txt")),S[idx][k])
#    end
#end

#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Het_Single= Coupling_Beam.init_Beam(ML(),false,false,false,false,true,false,max_level=5,MatlabSampler,folder,folder_with_elements,false,nterms=586,numberoftol=6,NQoI=641)
#estimator=init_Beam_mlmc_L_Het_Single
#history = run(estimator,5.71220000000000e-05)
#S= MultilevelEstimators.samples_diff(estimator)
#name = estimator[:name]
#sf_dir= string(name[1:end-5])
#isdir(joinpath(pathToInterm,sf_dir,"Samples")) || mkdir(joinpath(pathToInterm,sf_dir,"Samples"))
#for idx in CartesianIndices(S)
#    idx_dir = joinpath(pathToInterm,sf_dir,"Samples", join(idx.I,"_"))
#    isdir(idx_dir) || mkdir(idx_dir)
#    for k in keys(S[idx])
#        writedlm(joinpath(idx_dir,string("samples_level_",k[1]-1,".txt")),S[idx][k])
#    end
#end

#HIGH
##MC all


#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Het_Single_high=Coupling_Beam.init_Beam(SL(),false,false,false,false,true,true,startlevel=3,numberoftol=10,MatlabSampler,folder,folder_with_elements,false,nterms=78,NQoI=1)
#estimator=init_Beam_MC_L_Het_Single_high
#history = run(estimator,2.0e-05)

#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Het_Single_high=Coupling_Beam.init_Beam(SL(),false,false,false,false,true,true,startlevel=3,numberoftol=10,MatlabSampler,folder,folder_with_elements,false,nterms=586,NQoI=1)
#estimator=init_Beam_MC_L_Het_Single_high
#history = run(estimator,2.0e-05)

#######QMC



@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
init_Beam_MC_L_Het_Single_high=Coupling_Beam.init_Beam(SL(),true,false,false,false,true,true,nb_of_warm_up_samples=2,nshifts=10,max_level=3,startlevel=3,numberoftol=10,MatlabSampler,folder,folder_with_elements,false,nterms=78,NQoI=1)
estimator=init_Beam_MC_L_Het_Single_high
history = run(estimator,2.0e-05)

@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
init_Beam_MC_L_Het_Single_high=Coupling_Beam.init_Beam(SL(),true,false,false,false,true,true,nb_of_warm_up_samples=2,nshifts=10,maxlevel=3,startlevel=3,numberoftol=10,MatlabSampler,folder,folder_with_elements,false,nterms=586,NQoI=1)
estimator=init_Beam_MC_L_Het_Single_high
history = run(estimator,2.0e-05)

################################################



#REF WORKING

#########MC
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Het_Single=Coupling_Beam.init_Beam(SL(),false,false,false,false,true,false,startlevel=3,numberoftol=10,MatlabSampler,folder,folder_with_elements,false,nterms=78,NQoI=1)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,2.0e-05)

#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Het_Single=Coupling_Beam.init_Beam(SL(),false,false,false,false,true,false,startlevel=3,numberoftol=7,MatlabSampler,folder,folder_with_elements,false,nterms=586,NQoI=1)
#estimator1=init_Beam_MC_L_Het_Single
#history = run(estimator1,4.39400000000000e-05)

#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Het_Single=Coupling_Beam.init_Beam(SL(),false,false,false,false,true,false,max_level=4,numberoftol=10,MatlabSampler,folder,folder_with_elements,false,nterms=586,NQoI=1)
#estimator2=init_Beam_MC_L_Het_Single
#history = run(estimator2,2.0e-05)

######QMC

#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Het_Single=Coupling_Beam.init_Beam(SL(),true,false,false,false,true,false,nb_of_warm_up_samples=2,nshifts=10,maxlevel=3,startlevel=3,numberoftol=10,MatlabSampler,folder,folder_with_elements,false,nterms=78,NQoI=1)
#estimator2=init_Beam_MC_L_Het_Single
#history = run(estimator2,2.0e-05)


#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Het_Single=Coupling_Beam.init_Beam(SL(),true,false,false,false,true,false,nb_of_warm_up_samples=2,nshifts=10,maxlevel=3,startlevel=3,numberoftol=7,MatlabSampler,folder,folder_with_elements,false,nterms=586,NQoI=1)
#estimator1=init_Beam_MC_L_Het_Single
#history = run(estimator1,4.39400000000000e-05)

@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
init_Beam_MC_L_Het_Single=Coupling_Beam.init_Beam(SL(),true,false,false,false,true,false,nb_of_warm_up_samples=2,nshifts=10,max_level=4,startlevel=4,numberoftol=10,MatlabSampler,folder,folder_with_elements,false,nterms=586,NQoI=1)
estimator2=init_Beam_MC_L_Het_Single
#history = run(estimator2,2.0e-05)

##############QMC


##########################################################################################################
