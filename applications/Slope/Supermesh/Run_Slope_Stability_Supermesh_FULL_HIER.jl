
#module Run
push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))
push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
push!(LOAD_PATH,"/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE")
using Distributed
using MultilevelEstimators
using Coupling_Slope_Supermesh_FULL_HIER
using Pkg
using MATLAB
using DelimitedFiles
using Revise
using LatticeRules

using DigitalNets


#./julia /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/Run.jl 2>&1 | tee /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/10112018_Run_Log.txt
numberOfProcs=36

addprocs(numberOfProcs)
println("Procs added")
println(numberOfProcs)

#@everywhere using Coupling
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE")))
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
@everywhere push!(LOAD_PATH,"/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE")

@everywhere using Coupling_Slope_Supermesh_FULL_HIER
@everywhere using MATLAB


#pathToInterm=joinpath("/","home","philippe",".julia","dev","Applications","applications","SPDE","data")

#pathToInterm=joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE","data")

pathToInterm="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE/data"

folder = string(pathToInterm,"/Interm10/Slope") # for report
folder_with_elements=string(pathToInterm,"/Mesh/Slope")
if(isdir(folder)==false)
mkdir(folder)
else
    rm(folder,recursive=true)
    mkdir(folder)
end



@everywhere MatlabSampler() = eval_string(string("Slope_Slover_JULIA_MATLAB_Supermesh(", myid(), ")"))

#---------------------------------------------------
#p refinement MLMC
#---------------


#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(ML(),false,false,false,true,true,startlevel=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=10,nterms=73)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,5e-05)

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

#







#H refinement MLMC
#-------------
#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(ML(),false,false,false,true,false,startlevel=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=10,nterms=101)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,1.75e-04)



#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(ML(),false,false,false,true,false,startlevel=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=10,nterms=73)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,5e-05)
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




#---------------------------------------------------
#p refinement MC
#-------------
#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(SL(),false,false,false,true,true,max_level=4,MatlabSampler,folder,folder_with_elements,false,numberoftol=10)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,4e-04)
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

#---------------------------------------------------
# h refinement MC
#-------------
#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(SL(),false,false,false,true,false,max_level=3,MatlabSampler,folder,folder_with_elements,false,numberoftol=4)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,0.0008446915750000001)
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


#--------------------
#MIMC
#------------------
#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(FT(2),false,false,false,true,false,startlevel=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=10,nb_of_warm_up_samples=10,nshifts=10)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,1.75e-04)

#H refinement MLMC
#-------------


#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(ML(),false,false,false,true,false,startlevel=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=10,nterms=10,max_level=4)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,5e-05)

#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(ML(),false,false,false,true,false,startlevel=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=10,nterms=210,max_level=4)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,5e-05)






#---------------------------------------------------
#p refinement MLMC
#---------------

#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(ML(),false,false,false,true,true,startlevel=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=10,nterms=210)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,5e-05)
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

#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(ML(),false,false,false,true,true,startlevel=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=10,nterms=10)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,5e-05)
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


#H refinement MLQMC
#-------------
#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(ML(),true,false,false,true,false,startlevel=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=10,nb_of_warm_up_samples=2,nshifts=10,nterms=210)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,5e-05)

#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(ML(),true,false,false,true,false,startlevel=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=10,nb_of_warm_up_samples=2,nshifts=10,nterms=10)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,5e-05)

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


#---------------------------------------------------
#p refinement MLQMC
#---------------




nterms=400
pt=LatticeRule(nterms)
#pt=LatticeRule(nterms)
#pt=DigitalNet64(nterms)
init_Beam_MC_L_Het_Single=Coupling_Slope_Supermesh_FULL_HIER.init_Slope(ML(),true,false,false,true,true,startlevel=0,max_level=6,MatlabSampler,folder,folder_with_elements,false,numberoftol=34,nb_of_warm_up_samples=2,nshifts=10,nterms=nterms,correlateOnlyDiffs=false,corr_len=0.2,smoothness=2.0,point_generator=pt)
estimator=init_Beam_MC_L_Het_Single
#println(Int64.(estimator.internals.sample_method_internals.generators[1][1].lattice_rule.z))
history = run(estimator,5e-07)

#init_Beam_MC_L_Het_Single=Coupling_Slope.init_Slope(ML(),true,false,false,true,true,startlevel=0,MatlabSampler,folder,folder_with_elements,false,numberoftol=10,nb_of_warm_up_samples=2,nshifts=10,nterms=10)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,5e-05)
#S= MultilevelEstimators.samples_diff(estimator)
#name = estimator[:name]
#sf_dir= string(name[1:end-5])
#isdir(joinpath("/vsc-hard-mounts/leuven-user/330/vsc33032/RunReports",sf_dir,"Samples")) || mkdir(joinpath("/vsc-hard-mounts/leuven-user/330/vsc33032/RunReports",sf_dir,"Samples"))
#for idx in CartesianIndices(S)
#    idx_dir = joinpath("/vsc-hard-mounts/leuven-user/330/vsc33032/RunReports",sf_dir,"Samples", join(idx.I,"_"))
#    isdir(idx_dir) || mkdir(idx_dir)
#    for k in keys(S[idx])
#        writedlm(joinpath(idx_dir,string("samples_level_",k[1]-1,".txt")),S[idx][k])
#    end
#end

#---------------------------------------------------
