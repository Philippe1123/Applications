module Coupling_Slope_Tent


using GaussianRandomFields,Interact,MultilevelEstimators,Interpolations,MATLAB,Pkg
using Distributions
using Distributed
using DelimitedFiles
#using Dierckx
using LinearAlgebra
using Statistics
using Random
using Dates
using ScatteredInterpolation
using Revise
using DigitalNets
using LatticeRules


macro get_arg(key_name, default_value)
    @eval get_arg(args::Dict{Symbol, Any}, ::Val{$key_name}) = haskey(args, $key_name) ? args[$key_name] : $default_value
end

#@get_arg :max_index_set_param 3

@get_arg :minpadding index->0

get_arg(args::Dict{Symbol,Any}, arg::Symbol) = get_arg(args, Val(arg))

get_arg(args::Dict{Symbol,Any}, arg::Val{T}) where T = throw(ArgumentError(string("in init_lognormal, invalid key ", T, " found")))

get_max_index_set(index_set, Int64) = get_index_set(index_set, Int64)

get_max_index_set(::SL, args) = [Level(args)]

get_max_index_set(::Union{AD, U}, args) = get_index_set(get_arg(args, :max_search_space), get_arg(args, :max_index_set_param))


function init_Slope(index_set::AbstractIndexSet, is_qmc::Bool, is_multiple_qoi::Bool, is_analyse::Bool,isField::Bool,isHigerOrderRefinement::Bool,MatlabSampler::Function,folder_Interm::String,folder_with_elements::String,GaussPoints::Bool; corr_len::T=1.0, smoothness::T=1.0, nterms::N=101,
     max_level::N=3, nshifts::N=1,nb_of_warm_up_samples::N=40,continuate::Bool=true,numberoftol::N=10,NQoI::N=1,correlateOnlyDiffs::Bool=false,point_generator::Random.AbstractRNG,kwargs...) where{T<:AbstractFloat,N<:Integer,V<:MSession}
    # println(do_regression)
    # Dimensions of beam
#isHigerOrderRefinement=true
println("ishigherorderRef ",isHigerOrderRefinement)
  println("number Of TOlerances ",numberoftol)
  println("Max number of levels ",max_level)

  nb_qoi=1

    max_index_set_param=max_level
    args = Dict{Symbol,Any}(kwargs)
    args[:index_set] = index_set
    indices = get_max_index_set(index_set, max_level)
    minpadding = get_arg(args, :minpadding)
    ## Gaussian random fields ##


    println(indices)

    Elements=Dict()
    Nodes=Dict()


    if(isa(index_set,SL)||isa(index_set,ML))

      if(isHigerOrderRefinement==false)
        Centers=Dict()
        isElementRefinement=true
      type="h"
      elseif(isHigerOrderRefinement==true)
      type="p"
      isElementRefinement=false
      end

    else
        type="hp"
        isElementRefinement=true
        isHigerOrderRefinement=true

       end

println("ElementRefinement ",isElementRefinement)
println("isHigerOrderRefinement ",isHigerOrderRefinement)

    NumberOfDicEntries=Dict()
    StringForAcces=Dict()
    ij=1
    LevelDict=Dict()
    Positions=Dict()

    for id in indices
      IntermPositions=Dict()
      IntermLevelDict=Dict()
      LevelDict[id]=IntermLevelDict
      Positions[id]=IntermPositions

    end


    for id in indices
    SizeId=length(id)
    Access="_";
    for i=1:SizeId
    Access=string(Access,id[i])
    end

   NumberOfDicEntries[id]=0
#%%%% written only for MLMC
    if(correlateOnlyDiffs==true)
    NumberOfDicEntries[id]=NumberOfDicEntries[id]+1

    Access_test=string("_",id,"_M",id,'_',"S",id)


    LevelDict[id][id]=Access_test

      for (key,value) in diff(id)
      #IntermLevelDict=LevelDict[key]
            #   println(IntermLevelDict)



      NumberOfDicEntries[key]=NumberOfDicEntries[key]+1
      Access_test=string("_",key,"_M",id,'_',"S",key)
      LevelDict[key][id]=Access_test
      end




    end



#println(LevelDict)
#println(NumberOfDicEntries)
#println(StringForAcces)
#println("-------------------------------------------------------------------------")

    PathElement=string(string(folder_with_elements,string("/",type,"_refinement/Elements_L",Access)),".txt")
    Handle_Elements=open(PathElement)
    Element=readdlm(Handle_Elements,Int);
    Element=Element[:,5:7]
    Elements[id]=Element



    if((isa(index_set,SL)||isa(index_set,ML))&&isHigerOrderRefinement==false)
    PathNodes=string(string(folder_with_elements,string("/",type,"_refinement/Nodes_L",Access)),".txt")
    PathElementCenters=string(string(folder_with_elements,string("/",type,"_refinement/ElementsCenter_L",Access)),".txt")
    Handle_PathElementCenters=open(PathElementCenters,"w")
    Handle_Nodes=open(PathNodes)
    Node=readdlm(Handle_Nodes);
    Node=Node[:,2:3]
    Nodes[id]=Node
    Center=compute_centers(Node,Element)
    Centers[id]=Center
    #Nodes[id]=[Centers[id];Nodes[id]] ###############################################################################################use centers to generate
    writedlm(Handle_PathElementCenters,Center)
    close(Handle_PathElementCenters)
    close(Handle_Nodes)
    println("h-ref")
    else
      PathNodes_elements=string(string(folder_with_elements,string("/",type,"_refinement/Nodes_L",Access)),".txt")
      Nodes_elements=readdlm(PathNodes_elements);
      Nodes_elements=Nodes_elements[:,2:3]


        if(correlateOnlyDiffs==false)

          PathNodes=string(string(folder_with_elements,string("/",type,"_refinement/GaussPoints_L",Access)),".txt")
          Handle_Nodes=open(PathNodes)
          Node=readdlm(Handle_Nodes);
          close(Handle_Nodes)
          Node=Node[:,2:3]
      #Nodes[id]=[Node;Nodes_elements]############update
          Nodes[id]=Node

        end




      println("p-ref and hp-ref")




    end


    close(Handle_Elements)

    end

  #  println(LevelDict)
  #  println(LevelDict[Level(0)])
  #  println(LevelDict[Level(1)])
  #  println(LevelDict[Level(2)])
  #  println(LevelDict[Level(3)])


  #  for id in indices
  #    println("Slaves of", id )
  #    for (key,value) in diff(id)
  #     println(LevelDict[key][id])
  #    end
  #  end



    if(correlateOnlyDiffs==true)
      Nodes=Dict()

      for id in indices
        Node_Dict=Dict()
        Nodes[id]=Node_Dict=Dict()
      end


      for id in indices

        PathNodes=string(string(folder_with_elements,string("/",type,"_refinement/GaussPoints_L",LevelDict[id][id])),".txt")
        Handle_Nodes=open(PathNodes)
        Node=readdlm(Handle_Nodes);
        close(Handle_Nodes)
        Node=Node[:,2:3]
        Nodes[id][id]=Node


        for (key,value) in diff(id)
       PathNodes=string(string(folder_with_elements,string("/",type,"_refinement/GaussPoints_L",LevelDict[key][id])),".txt")
       Handle_Nodes=open(PathNodes)
       Node=readdlm(Handle_Nodes);
       close(Handle_Nodes)
       Node=Node[:,2:3]
       Nodes[key][id]=Node


    end
      end
    end






    #   distributions = [MultilevelEstimators.Normal() for i in 1:nterms]
    distributions = [MultilevelEstimators.Uniform() for i in 1:nterms]
println(distributions[1])
       p=2
      # exp_field = GaussianRandomFields.SquaredExponential(corr_len,σ=smoothness,p=p)
     #  exp_field = GaussianRandomFields.Exponential(corr_len,σ=smoothness,p=p)
       exp_field = GaussianRandomFields.Matern(corr_len,smoothness,σ=1.0,p=p)

       println("P of covar equals")
       println(p)
       cov = CovarianceFunction(2,exp_field)


         # all other levels
       grfs=Dict()
       if(correlateOnlyDiffs==false)
            for index in indices
                println(index)
                Random.seed!(1234)
                grfs[index] = GaussianRandomField(cov,KarhunenLoeve(nterms),Nodes[index],Elements[index],quad=GaussLegendre())
                println(grfs[index])
              end


              for index in indices
                if(index>Level(0))
                  index_1=index[1]-1
                  index_1=Level(index_1)
                  println(index)

                  if(grfs[index].data.eigenfunc[1,1]==grfs[index_1].data.eigenfunc[1,1])
                    println("Level ", index," and Level ",index_1, "same eigenfunc")
  #                Random.seed!(7328)
    #              println(GaussianRandomFields.sample(grfs[index])[13])
  else
    println("Warning: Level ", index," and Level ",index_1, "DIFFERENT eigenfunc")

                end
                end
                end
        else


          for id in indices
            grfs_Dict=Dict()
            grfs[id]=grfs_Dict
          end


          for index in indices
            println(string("Master Field for Level: ",index))
            Random.seed!(1234)
            grfs[index][index]=GaussianRandomField(cov,KarhunenLoeve(nterms),Nodes[index][index],Elements[index],quad=GaussLegendre(),eigensolver=EigenSolver())
            Random.seed!(1234)
            println(grfs[index][index])
            for (key,value) in diff(index)
              println(string("Slave Fields for Level: ",index))
          #    println(Nodes)
          Random.seed!(1234)
              grfs[key][index] = GaussianRandomField(cov,KarhunenLoeve(nterms),Nodes[key][index],Elements[index],quad=GaussLegendre(),eigensolver=EigenSolver())
              Random.seed!(1234)
#              println("here")
           println(grfs[key][index])
           if(grfs[key][index].data.eigenfunc[1,1]==grfs[index][index].data.eigenfunc[1,1])
             println("Eigenfunctions Matching")
           else
             println("WARNING: Eigenfunctions Not Matching, Program will switch to interpolate")

           end

            end
            println("---------------------------------------------")

            end

        end

        println(grfs)



      if(isa(index_set,SL))
      increment = max_level
      else
      increment=0
      end





         if is_qmc
           sample_method=QMC()
         else
             sample_method=MC()
         end

          # name
        name = "Slope "
        name = is_analyse ? string(name,"analyse ") : name
        name = isa(index_set,AD) ? string(name,"A") : name
        name = isa(index_set,ML) ? string(name,"ML") : MultilevelEstimators.ndims(index_set) > 1 ? string(name,"MI") : name
        name = is_qmc ? string(name,"Q") : name
        name = string(name,"MC")
        name = isField ? string(name,"_Het") : string(name,"_Hom")
        name = isHigerOrderRefinement ? string(name,"_High") : string(name,"_Ref")
        name = GaussPoints ? string(name,"_GP") :
        name = is_multiple_qoi ? string(name," (multiple)") : name
        name = string(name,"_")
        timenow = Dates.now()
        timenow = Dates.format(timenow, "dd-mm-yyyy-T:HH:MM:SS")
        name = string(name,timenow)
        rd=rand(1:1000000,1,1)
        name = string(name,rd)
        #nb_of_qoi = is_multiple_qoi ? Int(Lx/he*2^(max_level-1)+1) : 1
        sample_function = correlateOnlyDiffs ? (index, ξ) -> Slope_cslo(index, ξ, grfs,Nodes,Elements,MatlabSampler,folder_Interm,folder_with_elements,isHigerOrderRefinement,isElementRefinement,increment,correlateOnlyDiffs,LevelDict) : (index, ξ) -> Slope(index, ξ, grfs,Nodes,Elements,MatlabSampler,folder_Interm,folder_with_elements,isHigerOrderRefinement,isElementRefinement,increment,correlateOnlyDiffs)
        isa(index_set,SL) ? println("All samples taken on level ", max_level) : println()



#   folder = string(joinpath(Pkg.dir("MultilevelEstimators"),"applications","SPDE","data",name)) # for report
#    folder = string("/home/philippe/.julia/dev/Applications/applications/SPDE/data/",name)
folder=string("/vsc-hard-mounts/leuven-user/330/vsc33032/RunReports/",name)


namefix=name
if(isdir(folder)==false)
mkdir(folder)
else
  ic=0
  while(isdir(folder)==true)
   name = string(namefix,ic)
   folder=string("/vsc-hard-mounts/leuven-user/330/vsc33032/RunReports/",name)
   ic=ic+1
  end
    mkdir(folder)
end


        ## Estimator ##
        #println(folder)
    if(!isa(index_set,SL))


        if(!isHigerOrderRefinement)
        γ = 2
        else
        γ = 1.5
        end
        println("Gamma is")
        println(γ)

    if(is_qmc)


        MultilevelEstimators.Estimator(
        index_set, # index_set: ML, SL, TD...
        sample_method,
        sample_function,
        distributions,
        name = name, # estimator name
        folder = folder, # for report
        do_mse_splitting=false,
        nb_of_shifts=nshifts,
        nb_of_warm_up_samples=nb_of_warm_up_samples,
        max_index_set_param=max_index_set_param,
        continuate=continuate,
        nb_of_tols=numberoftol,
        save_samples = false,
         nb_of_qoi = NQoI,
         point_generator=point_generator,
         sample_mul_factor=2.0,
          #user_data = grfs, # GRF's
         #verbose = true, # display information
         #nb_of_qoi = nb_of_qoi, # number of qoi
         #cost_model = (index) -> geometric_cost_model(4,1.5,index), # cost model
         #sample_multiplication_factor = sample_multiplication_factor, # qmc multiplication factor
         #store_samples=false,
         )
    else

        MultilevelEstimators.Estimator(
        index_set, # index_set: ML, SL, TD...
        sample_method,
        sample_function,
        distributions,
        name = name, # estimator name
        folder = folder, # for report
        do_mse_splitting=false,
        nb_of_warm_up_samples=nb_of_warm_up_samples,
        max_index_set_param=max_index_set_param,
        continuate=continuate,
        nb_of_tols=numberoftol,
        save_samples = false,
        nb_of_qoi = NQoI,
         #user_data = grfs, # GRF's
         #verbose = true, # display information
         #nb_of_qoi = nb_of_qoi, # number of qoi
         #cost_model = (index) -> geometric_cost_model(4,1.5,index), # cost model
         #sample_multiplication_factor = sample_multiplication_factor, # qmc multiplication factor
         #store_samples=false,
         )

    end
    else
      checkPtgen= @isdefined point_generator

    if(checkPtgen==true)
      MultilevelEstimators.Estimator(
      index_set, # index_set: ML, SL, TD...
      sample_method,
      sample_function,
      distributions,
      name = name, # estimator name
      folder = folder, # for report
      nb_of_warm_up_samples=nb_of_warm_up_samples,
  #    max_index_set_param=max_index_set_param,
      continuate=continuate,
      nb_of_tols=numberoftol,
      nb_of_qoi = NQoI,
      point_generator=point_generator,
      sample_mul_factor=2.0,
      nb_of_shifts=nshifts,
  # number of qoi
       )

   else
       MultilevelEstimators.Estimator(
       index_set, # index_set: ML, SL, TD...
       sample_method,
       sample_function,
       distributions,
       name = name, # estimator name
       folder = folder, # for report
       nb_of_warm_up_samples=nb_of_warm_up_samples,
   #    max_index_set_param=max_index_set_param,
       continuate=continuate,
       nb_of_tols=numberoftol,
       nb_of_qoi = NQoI,
        )
   end
   end




end

## user data ##
struct Field_Data{V}
    fields::V
end





function Slope(index::Index, ξ::Vector{T} where {T<:Real}, grf::Dict, Nodes::Dict, Elements::Dict,MatlabSampler::Function,folder::String,folder_with_elements::String,isHigerOrderRefinement::Bool,isElementRefinement::Bool,increment::Int64,correlateOnlyDiffs::Bool)

    StepChange=0
   if(length(increment)==length(index))
    index=index+Level(increment)
      if(index[1]>=7)
      StepChange=1

      end
    end

    Nodes_Fine=Nodes[index]
    Elements_Fine=Elements[index]
    Zf = GaussianRandomFields.sample(grf[index],xi=ξ) # compute GRF
    Zf=expfield(Zf)

    Qf = Slope_Sample(Zf,MatlabSampler,folder,index,false,isHigerOrderRefinement,isElementRefinement,StepChange)

    dQ = Qf
    if(increment==0)
    for (key,value) in diff(index)

    #    index_1=index[1]-1
        index_1=key
        Nodes_Coarse=Nodes[index_1]
        Elements_Coarse=Elements[index_1]
#        itp = ScatteredInterpolation.interpolate(NearestNeighbor(), Nodes_Fine', Zf);
#        Zc=evaluate(itp, Nodes_Coarse')
#        Zc=vec(Zc)
   if(grf[index_1].data.eigenfunc[1,1]==grf[index].data.eigenfunc[1,1])
         Zc=GaussianRandomFields.sample(grf[index_1],xi=ξ)
        Zc=expfield(Zc)
      else
        itp = ScatteredInterpolation.interpolate(NearestNeighbor(), Nodes_Fine', Zf);
                Zc=evaluate(itp, Nodes_Coarse')
                Zc=vec(Zc)

      end
        Qc = Slope_Sample(Zc,MatlabSampler,folder,index_1,true,isHigerOrderRefinement,isElementRefinement,StepChange)
        dQ += value*Qc
  #      println(float(Qf)," break ",float(Qc)," break ",float(dQ)," break ",value)
    end
  end
    return (dQ,Qf)
end


## sample functions ##
function Slope_Sample(Z::Vector{T},MatlabSampler::Function,folder::String,index::Index,doRestriction::Bool,isHigerOrderRefinement::Bool,isElementRefinement::Bool,StepChange::Int64) where {T<:Real}

    ind="_"
    arr=[]
for len=1:length(index)
     append!(arr,index[len])
     ind=string(ind,index[len])
end
ind=string(ind,"_")

       if(isHigerOrderRefinement==true&&isElementRefinement==false)
      isHigerOrderRefinement_num=1
  elseif(isHigerOrderRefinement==false&&isElementRefinement==true)
      isHigerOrderRefinement_num=0
  elseif(isHigerOrderRefinement==true&&isElementRefinement==true)
      isHigerOrderRefinement_num=3
    end

        Res=open(string(folder,string("/Res_",myid())),"w")
        close(Res)


        stp=open(string(folder,string("/StepChange_",myid())),"w")
         writedlm(stp,StepChange)
         close(stp)

        indx=open(string(folder,string("/Index_",myid())),"w")
         writedlm(indx,arr)
         close(indx)

#        println(folder)
#        println(ind)

         ZFile=open(string(folder,"/Cohesion",ind,myid()),"w")
          writedlm(ZFile,Z)
          close(ZFile)




          HighOrder=open(string(folder,string("/HighOrder_",myid())),"w")
          writedlm(HighOrder,isHigerOrderRefinement_num)
          close(HighOrder)



          if doRestriction==true
              numres=1
          else
              numres=0
          end

           RestHandle=open(string(folder,string("/doRestriction_",myid())),"w")
           writedlm(RestHandle,numres)
           close(RestHandle)



           t=@elapsed     MatlabSampler()
            Res=open(string(folder,string("/Res_",myid())));
            #println("level ",index[1]," ","time",t)



            u=readdlm(Res);
            close(Res)



        #    u=u[1]


     return u
end








#cslo= correlate only successive levels
function Slope_cslo(index::Index, ξ::Vector{T} where {T<:Real}, grfs::Dict, Nodes::Dict, Elements::Dict,MatlabSampler::Function,folder::String,folder_with_elements::String,isHigerOrderRefinement::Bool,isElementRefinement::Bool,increment::Int64,correlateOnlyDiffs::Bool,Extension::Dict)

#println("Slope_cslo")
    StepChange=0
   if(length(increment)==length(index))
    index=index+Level(increment)
      if(index[1]>=7)
      StepChange=1

      end
    end


    Nodes_Fine=Nodes[index][index]

    Elements_Fine=Elements[index]
    grf=grfs[index][index]



    ω=ones(length(ξ),1)
	x=ones(length(ξ),1)-abs.(2*ξ.-1)
	for id=1:length(x)
		ω[id]=transform(MultilevelEstimators.TruncatedNormal(0,1,-2,2),x[id])
	end




    Zf = GaussianRandomFields.sample(grf,xi=ω) # compute GRF
    Zf=expfield(Zf)
    StringExtension=Extension[index][index][4:end]
    Qf = Slope_Sample_cslo(Zf,MatlabSampler,folder,index,false,isHigerOrderRefinement,isElementRefinement,StepChange,StringExtension)

    dQ = Qf
    if(increment==0)
    for (key,value) in diff(index)

    #    index_1=index[1]-1
        index_1=key
        Nodes_Coarse=Nodes[index_1][index]
        Elements_Coarse=Elements[index_1]
        StringExtension=Extension[index_1][index][4:end]

  if(grfs[index_1][index].data.eigenfunc[1,1]==grfs[index][index].data.eigenfunc[1,1])
        Zc=GaussianRandomFields.sample(grfs[index_1][index],xi=ω)
       Zc=expfield(Zc)
  #      println("success")
      else
  #      println("fail")
       itp = ScatteredInterpolation.interpolate(NearestNeighbor(), Nodes_Fine', Zf);
        Zc=evaluate(itp, Nodes_Coarse')
        Zc=vec(Zc)
        end


        Qc = Slope_Sample_cslo(Zc,MatlabSampler,folder,index_1,true,isHigerOrderRefinement,isElementRefinement,StepChange,StringExtension)
        dQ += value*Qc
  #      println(float(Qf)," break ",float(Qc)," break ",float(dQ)," break ",value)
#  println("loop end ")

    end
  end
    return (dQ,Qf)
end




## sample functions ##
function Slope_Sample_cslo(Z::Vector{T},MatlabSampler::Function,folder::String,index::Index,doRestriction::Bool,isHigerOrderRefinement::Bool,isElementRefinement::Bool,StepChange::Int64,StringExtension::String) where {T<:Real}

    ind="_"
    arr=[]
for len=1:length(index)
     append!(arr,index[len])
     ind=string(ind,index[len])
end
ind=string(ind,StringExtension,"_")

       if(isHigerOrderRefinement==true&&isElementRefinement==false)
      isHigerOrderRefinement_num=1
  elseif(isHigerOrderRefinement==false&&isElementRefinement==true)
      isHigerOrderRefinement_num=0
  elseif(isHigerOrderRefinement==true&&isElementRefinement==true)
      isHigerOrderRefinement_num=3
    end

  #   println(StringExtension)
      Ext=open(string(folder,string("/Extension_",myid())),"w")
      write(Ext,StringExtension)
      close(Ext)

        Res=open(string(folder,string("/Res_",myid())),"w")
        close(Res)


        stp=open(string(folder,string("/StepChange_",myid())),"w")
         writedlm(stp,StepChange)
         close(stp)

        indx=open(string(folder,string("/Index_",myid())),"w")
         writedlm(indx,arr)
         close(indx)

#        println(folder)
#        println(ind)

         ZFile=open(string(folder,"/Cohesion",ind,myid()),"w")
          writedlm(ZFile,Z)
          close(ZFile)




          HighOrder=open(string(folder,string("/HighOrder_",myid())),"w")
          writedlm(HighOrder,isHigerOrderRefinement_num)
          close(HighOrder)



          if doRestriction==true
              numres=1
          else
              numres=0
          end

           RestHandle=open(string(folder,string("/doRestriction_",myid())),"w")
           writedlm(RestHandle,numres)
           close(RestHandle)



           t=@elapsed     MatlabSampler()
            Res=open(string(folder,string("/Res_",myid())));
            #println("level ",index[1]," ","time",t)



            u=readdlm(Res);
            close(Res)
        #    u=u[1]


     return u
end








function compute_centers(p,t)
    d = size(p, 2)
    vec_t = vec(t)
    size_t = size(t)

    pts = Array{Float64}(undef, size(t, 1), d)
    @inbounds for i in 1:d
        x = reshape(p[vec_t, i], size_t)
        mean!(view(pts, :, i), x)
    end
    pts
end


function expfield(Z::Vector{T})where {T<:Real}

Mean=log(8000)
σ=0.05

Z=Z.*σ
Z=Z.+Mean
Z=exp.(Z)

return Z
end
end
