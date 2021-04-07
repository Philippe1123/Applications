__precompile__()

 module Coupling_Beam

 #sample_function = isa(index_set,SL) ? isHigerOrderRefinement ? (index,ξ)->Beam_NL_Het_MC_High(index, ξ,Lx,Ly, grfs[index],startlevel) : (index,ξ)->Beam_NL_Het_MC_Ref(index, ξ,Lx,Ly, grfs[index]) : isNL ? isField ? isHigerOrderRefinement ? (index,ξ) -> Beam_NL_Het_High(index, ξ,Lx,Ly, grfs[index]) :  (index,ξ) -> Beam_NL_Het(index, ξ,Lx,Ly, grfs[index]) :  isHigerOrderRefinement ? (index, ξ) -> Beam_NL_Hom_High(index, ξ,Lx,Ly, grfs[index]) : (index, ξ) -> Beam_NL_Hom(index, ξ,Lx,Ly, grfs[index])  : isField ? is_multiple_qoi ? (index,ξ,data)->Beam_L_Het_multiple(index,ξ,data,nb_of_qoi) : isHigerOrderRefinement ? (index, ξ) -> Beam_L_Het_High(index, ξ,Lx,Ly, grfs[index]) : (index, ξ) -> Beam_L_Het(index, ξ,Lx,Ly, grfs[index]) : is_multiple_qoi ? (index,ξ,data)->Beam_L_Het_multiple(index,ξ,data,nb_of_qoi) : isHigerOrderRefinement ? (index, ξ) -> Beam_L_Hom_High(index, ξ,Lx,Ly, grfs[index]) : (index, ξ) -> Beam_L_Hom(index, ξ,Lx,Ly, grfs[index])



#push!(LOAD_PATH,Pkg.dir(joinpath("Coupling","src")))
#push!(LOAD_PATH,Pkg.dir(joinpath("MultilevelEstimators","applications","SPDE")))
  using   GaussianRandomFields,Interact, FieldTransformation,MultilevelEstimators,Interpolations,Solver_L,HomogeneousNormalMatrixGen,MATLAB,Pkg
  using Distributions
  using Distributed
  using DelimitedFiles
  using LinearAlgebra
  using Statistics
  using Random
  using Dates
using ScatteredInterpolation
using DigitalNets
using LatticeRules


 # @reexport using MultilevelEstimators, GaussianRandomFields,HomogeneousNormalMatrixGen

#Distributions,
include("FieldTransformation.jl")

## import statements ##
 #import Base.getindex

##export statements
#export init_Beam_Test
## Continuation == true, samples on level 0 1 2 and then extrap
macro get_arg(key_name, default_value)
    @eval get_arg(args::Dict{Symbol, Any}, ::Val{$key_name}) = haskey(args, $key_name) ? args[$key_name] : $default_value
end

#@get_arg :max_index_set_param 6

@get_arg :minpadding index->0

get_arg(args::Dict{Symbol,Any}, arg::Symbol) = get_arg(args, Val(arg))

get_arg(args::Dict{Symbol,Any}, arg::Val{T}) where T = throw(ArgumentError(string("in init_lognormal, invalid key ", T, " found")))

get_max_index_set(index_set, Int64) = get_index_set(index_set, Int64)

get_max_index_set(::SL, args) = [Level(args)]

get_max_index_set(::Union{AD, U}, args) = get_index_set(get_arg(args, :max_search_space), get_arg(args, :max_index_set_param))






function init_Beam(index_set::AbstractIndexSet, is_qmc::Bool, is_multiple_qoi::Bool, is_analyse::Bool,isNL::Bool,isField::Bool,isHigerOrderRefinement::Bool,MatlabSampler::Function,folder_Interm::String,folder_with_elements::String,GaussPoints::Bool; corr_len::T=0.3, smoothness::T=1.0, nterms::N=1,
     max_level::N=3, nshifts::N=1,nb_of_warm_up_samples::N=40,continuate::Bool=true,startlevel::N=0,numberoftol::N=10,NQoI::N=641,ptgen::Random.AbstractRNG,kwargs...) where{T<:AbstractFloat,N<:Integer,V<:MSession}
    # println(do_regression)
    # Dimensions of beam

     println(GaussPoints)
     println(numberoftol)

#     Lx = 2.5;
#     Ly = 0.25;

#     nelx=Lx/he
#     nely=Ly/he

lvlM=open(string(folder_Interm,"/MaxLVL"),"w")
writedlm(lvlM,max_level)
close(lvlM)

NQ=open(string(folder_Interm,"/NQoI"),"w")
writedlm(NQ,NQoI)
close(NQ)

     max_index_set_param=max_level
     args = Dict{Symbol,Any}(kwargs)
     args[:index_set] = index_set
     indices = get_max_index_set(index_set, max_level)
     minpadding = get_arg(args, :minpadding)
     ## Gaussian random fields ##


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


     for id in indices

     SizeId=length(id)
     Access="_";
     for i=1:SizeId
     Access=string(Access,id[i])
     end



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
#     Nodes[id]=Center ###############################################################################################use centers to generate RF
     writedlm(Handle_PathElementCenters,Center)
     close(Handle_PathElementCenters)
     close(Handle_Nodes)
     println("h-ref")
     else
       PathNodes=string(string(folder_with_elements,string("/",type,"_refinement/GaussPoints_L",Access)),".txt")
       Handle_Nodes=open(PathNodes)
       Node=readdlm(Handle_Nodes);
       Node=Node[:,2:3]
       Nodes[id]=Node
       close(Handle_Nodes)
       println("p-ref and hp-ref")

     end


     close(Handle_Elements)

     end








         coarse_dof = 1
#         distributions = [MultilevelEstimators.Normal() for i in 1:nterms]
         distributions = [MultilevelEstimators.TruncatedNormal(0,1,-2,2) for i in 1:nterms]
println(distributions[1])

           if(isField)

             p=2
        #     exp_field = GaussianRandomFields.Exponential(corr_len,σ=smoothness,p=p)
             exp_field = GaussianRandomFields.Matern(corr_len,smoothness,σ=1,p=p)

             println("P of covar equals")
             println(p)

            cov = CovarianceFunction(2,exp_field)


# Note fields are generated as vy x vy thus 160 x 40 example the solve transposes this

grfs=Dict()
for index in indices
 println(index)
 #Random.seed!(1234)
 grfs[index] = GaussianRandomField(cov,KarhunenLoeve(nterms),Nodes[index],Elements[index],quad=GaussLegendre())
 println(grfs[index])
end



#println(isField)
elseif(!isField)
############## for loop is working implement it everywhere
 if(isHigerOrderRefinement==false)
             grfs=Dict()
             i=startlevel
             for index in indices


                 j=i
                 m = coarse_dof*2^i
                 n = coarse_dof*2^j

                 vx = 0:(he/m)/Lx:0.99999
                 vy= 0:(he/n)/Ly:0.99999
                 grfs[index] = HomogeneousNormalMatrixGen.HomogeneousNormalMatrix(vx,vy,(vx,vy))
                 i=i+1
             end

else

    grfs=Dict()
    i=0
    for index in indices


        j=i
        m = coarse_dof*2^i
        n = coarse_dof*2^j

        vx = 0:(he/m)/Lx:0.99999
        vy= 0:(he/n)/Ly:0.99999
        grfs[index] = HomogeneousNormalMatrixGen.HomogeneousNormalMatrix(vx,vy,(vx,vy))
    end


end

end


if(isa(index_set,ML))
increment = 0
elseif(isa(index_set,SL))
increment=max_level
end


     if is_qmc
       sample_method=QMC()
     else
         sample_method=MC()
     end

      # name
    name = "Beam "
    name = is_analyse ? string(name,"analyse ") : name
    name = isa(index_set,AD) ? string(name,"A") : name
    name = isa(index_set,ML) ? string(name,"ML") : MultilevelEstimators.ndims(index_set) > 1 ? string(name,"MI") : name
    name = is_qmc ? string(name,"Q") : name
    name = string(name,"MC")
    name = isField ? string(name,"_Het") : string(name,"_Hom")
    name = isNL ? string(name,"_NonLin") : string(name,"_Lin")
    name = isHigerOrderRefinement ? string(name,"_High") : string(name,"_Ref")
    name = GaussPoints ? string(name,"_GP") :

    name = is_multiple_qoi ? string(name," (multiple)") : name

    timenow = Dates.now()
    timenow = Dates.format(timenow, "dd-mm-yyyy-T:HH:MM:SS")
    name = string(name,timenow)
    nb_of_qoi = is_multiple_qoi ? Int(Lx/he*2^(max_level-1)+1) : 1
    sample_function = isNL ? isField ? (index,ξ)->Beam_NL_Het_ALL(index, ξ,Lx,Ly , grfs[index],MatlabSampler,folder_Interm,isHigerOrderRefinement,increment) : (index,ξ)->Beam_NL_Hom_ALL(index, ξ,Lx,Ly , grfs[index],MatlabSampler,folder_Interm,isHigerOrderRefinement,increment) : isField ? (index,ξ)->Beam_L_Het_ALL(index, ξ, grfs,Nodes,Elements,MatlabSampler,folder_Interm,folder_with_elements,isHigerOrderRefinement,isElementRefinement,increment) : (index,ξ)->Beam_L_Hom_ALL(index, ξ,Lx,Ly , grfs[index],MatlabSampler,folder_Interm,isHigerOrderRefinement,increment)


#println(startlevel)
#println(sample_function)


#folder = string(joinpath(Pkg.dir("MultilevelEstimators"),"applications","SPDE","data",name)) # for report
folder=string("/vsc-hard-mounts/leuven-user/330/vsc33032/RunReports/",name)
if(isdir(folder)==false)
mkdir(folder)
else
    rm(folder,recursive=true)
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
    #cost_model = level -> 2^(γ * level[1]),
     #user_data = grfs, # GRF's
     #verbose = true, # display information
     nb_of_qoi = NQoI, # number of qoi
     point_generator=ptgen,
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
    #cost_model = level -> 2^(γ * level[1]),
     #user_data = grfs, # GRF's
     #verbose = true, # display information
     nb_of_qoi = NQoI, # number of qoi
     #cost_model = (index) -> geometric_cost_model(4,1.5,index), # cost model
     #sample_multiplication_factor = sample_multiplication_factor, # qmc multiplication factor
     #store_samples=false,
     )

end
else

    checkPtgen= @isdefined ptgen

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
    point_generator=ptgen,
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





function Beam_NL_Het_ALL(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::GaussianRandomField,MatlabSampler::Function,folder::String,isHigerOrderRefinement::Bool,increment::Int64)

    index=index+Level(increment)

    Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
    Zf=FieldTransformation.Transform_NL(Zf)
    he=Lx/size(Zf,1)
    Qf = Beam_sample_NL(Zf,Lx,Ly,he,index[1],MatlabSampler,folder)
    dQ = Qf
    if(increment==0)
    for (key,value) in diff(index)
        step = (index - key).I .+ 1
        Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
        he=Lx/size(Zc,1)
        Qc = Beam_sample_NL(Zc,Lx,Ly,he,(index[1]-1),MatlabSampler,folder)
        dQ += value*Qc
    end
    end
    return (dQ,Qf)

end

function Beam_NL_Hom_ALL(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::HomogeneousNormalMatrixGen.HomogeneousNormalMatrix,MatlabSampler::Function,folder::String,isHigerOrderRefinement::Bool,increment::Int64)

    ξ_number=ξ[1]

    index=index+Level(increment)

    Zf = HomogeneousNormalMatrixGen.sample(grf,ξ_number)
    Zf=FieldTransformation.Transform_NL(Zf)
    he=Lx/size(Zf,1)
    Qf = Beam_sample_NL(Zf,Lx,Ly,he,index[1],MatlabSampler,folder)
    dQ = Qf
    if(increment==0)
    for (key,value) in diff(index)
        step = (index - key).I .+ 1

        Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
        he=Lx/size(Zc,1)
        Qc = Beam_sample_NL(Zc,Lx,Ly,he,(index[1]-1),MatlabSampler,folder)
        dQ += value*Qc
    end
    end
    return (dQ,Qf)

end

function Beam_L_Het_ALL(index::Index, ξ::Vector{T} where {T<:Real}, grf::Dict, Nodes::Dict, Elements::Dict,MatlabSampler::Function,folder::String,folder_with_elements::String,isHigerOrderRefinement::Bool,isElementRefinement::Bool,increment::Int64)
    index=index+Level(increment)
#    println("Finelevel ",index)
    Nodes_Fine=Nodes[index]
    Elements_Fine=Elements[index]
    Zf = GaussianRandomFields.sample(grf[index],xi=ξ[1:randdim(grf[index])]) # compute GRF


    #Zf=FieldTransformation.Transform_L(Zf)
    # Zf=30E3.+6E3*Zf


    #Z=sample(grfs, xi=ω)
	Mean=log(30E3)
	σ=0.1
	Zf=Zf.*σ
	Zf=Zf.+Mean
	Zf=exp.(Zf)


    Zf=vec(Zf)
    Qf = Beam_sample_L(Zf,MatlabSampler,folder,index,false,isHigerOrderRefinement,isElementRefinement)

    dQ = Qf
    if(increment==0)
    for (key,value) in diff(index)
        index_1=key
        Nodes_Coarse=Nodes[index_1]
        #println(size(Zf))
        #println(size(Nodes_Fine))
        itp = ScatteredInterpolation.interpolate(NearestNeighbor(), Nodes_Fine', Zf);
        Zc=evaluate(itp, Nodes_Coarse')
        Zc=vec(Zc)
        #step = (index - key).I .+ 1
        #Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
        Qc = Beam_sample_L(Zc,MatlabSampler,folder,index_1,true,isHigerOrderRefinement,isElementRefinement)

        dQ += value*Qc
    end
end

GC.gc()

return (dQ,Qf)

end

function Beam_L_Hom_ALL(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::HomogeneousNormalMatrixGen.HomogeneousNormalMatrix,MatlabSampler::Function,folder::String,isHigerOrderRefinement::Bool,increment::Int64)
    ξ_number=ξ[1]
    Zf = HomogeneousNormalMatrixGen.sample(grf,ξ_number) # compute GRF
    Zf=FieldTransformation.Transform_L(Zf)
    he=Lx/size(Zf,1)
    Qf = Beam_sample_L(Zf,Lx,Ly,he,index[1],MatlabSampler,folder)
    dQ = Qf
    if(increment==0)
    for (key,value) in diff(index)
        step = (index - key).I .+ 1

        Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
        he=Lx/size(Zc,1)
        Qc = Beam_sample_L(Zc,Lx,Ly,he,(index[1]-1),MatlabSampler,folder)
        dQ += value*Qc
    end
end

GC.gc()

return (dQ,Qf)


end




function Beam_sample_NL(Z::Matrix{T},Lx::T,Ly::T,he::T,level::Int64,MatlabSampler::Function,folder::String) where {T<:Real}
    nu = 0.25;  #move to matlab
    fy = 240.0; #move to matlab
    t = 1.0;  #move to matlab
    Lx=Lx*1000
    Ly=Ly*1000
    he=he*1000

    Res=open(string(folder,string("/Res_",myid())),"w")
    close(Res)

     LxFile=open(string(folder,string("/Lx_",myid())),"w")
     LyFile=open(string(folder,string("/Ly_",myid())),"w")
     heFile=open(string(folder,string("/he_",myid())),"w")
     ZFile=open(string(folder,string("/E_",myid())),"w")
     LevelFile=open(string(folder,string("/Lev_",myid())),"w")

      writedlm(LxFile,Lx)
      writedlm(LyFile,Ly)
      writedlm(heFile,he)
      writedlm(ZFile,Z)
      writedlm(LevelFile,level)

      close(LxFile)
      close(LyFile)
      close(heFile)
      close(ZFile)
      close(LevelFile)



   t=@elapsed     MatlabSampler()
   #println(t)
        Res=open(string(folder,string("/Res_",myid())));
        u=readdlm(Res);


   u=abs.(u/1000);
   return u
end









     function Beam_sample_L(Z::Vector{T},MatlabSampler::Function,folder::String,index::Index,doRestriction::Bool,isHigerOrderRefinement::Bool,isElementRefinement::Bool) where {T<:Real}

     nu=0.15;
     t = 1000.0;


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

         indx=open(string(folder,string("/Index_",myid())),"w")
          writedlm(indx,arr)
          close(indx)

 #        println(folder)
 #        println(ind)

          ZFile=open(string(folder,"/E",ind,myid()),"w")
           writedlm(ZFile,Z)
           close(ZFile)




           HighOrder=open(string(folder,string("/HighOrder_",myid())),"w")
           writedlm(HighOrder,isHigerOrderRefinement_num)
           close(HighOrder)


     Res=open(string(folder,string("/Res_",myid())),"w")
     close(Res)

     if doRestriction==true
         numres=1
     else
         numres=0
     end

      RestHandle=open(string(folder,string("/doRestriction_",myid())),"w")
      writedlm(RestHandle,numres)
      close(RestHandle)


t=@elapsed     MatlabSampler()
#println(t)
         Res=open(string(folder,string("/Res_",myid())));
         u=readdlm(Res);
GC.gc()

 u=abs.(u/1000);

  return u
     end





function Extract_matrix_Displacements_From_Vector(Points::Array{Float64,1},Lx::Float64,Ly::Float64,he::Float64)
nelx=Int64(Lx/he)
nely=Int64(Ly/he)
u=Points
#matxDir::Array{Float64,2}
#matyDir::Array{Float64,2}
#Matx::Array{Float64,2}
#Maty::Array{Float64,2}
        matxDir=zeros(nely+1,nelx+1);
        matyDir=zeros(nely+1,nelx+1);
        x=1;  y=1;  z=2;  t=1;v=1;
        while(x<=(nelx+1))
            y=1;
            while(y<=(nely+1))
                    matxDir[y,x]=u[t];
                    matyDir[y,x]=u[z];
                    z=z+2;  t=t+2;

                v=v+1;
                y=y+1;
            end
            x=x+1;
        end
        Matx=matxDir;
        Maty=matyDir;

        return Matx,Maty

end

function Beam_L_Het_High_GP(index::Index, ξ::Vector{T} where {T<:Real}, Lx::Float64,Ly::Float64 , grf::GaussianRandomField,MatlabSampler::Function,folder::String,nelx::Int64,nely::Int64,he::Float64)
    Error=true

    while(Error ==true)
    #    println("Start sampling")
       o(Lx,Ly,index,grf,Error)=try
    #       println("Start try catch")
        Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
        Zf=FieldTransformation.Transform_L(Zf)
        Qf = Beam_sample_L_High(Zf,Lx,Ly,he,index[1],MatlabSampler,folder)
        dQ = Qf
        for (key,value) in diff(index)
            step = (index - key).I .+ 1

#            println(index[1]+1)
#            println(index[1])
#            println(nely)
#            println(size(Zf))
            Ax=collect(0:1/(index[1]+2):nelx-1/(index[1]+2))
            Ay=collect(0:1/(index[1]+2):nely-1/(index[1]+2))
#          println(size(Ax))
#          println((Ax))

            knots=(Ax,Ay,)
            itp = Interpolations.interpolate(knots, Zf, Gridded(Linear()))

            if(index[1]==1)
#             Ax_1=collect(0:1:nelx-1)
#             Ax_2=collect(2/3:1:nelx)

#             Ax_1=collect(0.2:1:nelx)
#             Ax_2=collect(0.6:1:nelx)

            Ax_1=collect(0.2236:1:nelx)
             Ax_2=collect(0.6:1:nelx)

             Ax_restricted=vcat(Ax_1,Ax_2)
             Ax_restricted=sort(vec(Ax_restricted))


#             Ay_1=collect(0:1:nely-1)
#             Ay_2=collect(2/3:1:nely)

#             Ay_1=collect(0.2:1:nely)
#             Ay_2=collect(0.6:1:nely)

             Ay_1=collect(0.2236:1:nely)
             Ay_2=collect(0.6:1:nely)

             Ay_restricted=vcat(Ay_1,Ay_2)
             Ay_restricted=sort(vec(Ay_restricted))

         elseif(index[1]==2)

    #         Ay_1=collect(0:1:nely-1)
    #         Ay_2=collect(0.375:1:nely)
    #         Ay_3=collect(0.75:1:nely)

             Ay_1=collect(0.0392:1:nely)
             Ay_2=collect(0.375:1:nely)
             Ay_3=collect(0.71:1:nely)

             Ay_restricted=vcat(Ay_1,Ay_2)
             Ay_restricted=vcat(Ay_restricted,Ay_3)
             Ay_restricted=sort(vec(Ay_restricted))


             Ax_1=collect(0.0392:1:nelx)
             Ax_2=collect(0.375:1:nelx)
             Ax_3=collect(0.71:1:nelx)
             Ax_restricted=vcat(Ax_1,Ax_2)
             Ax_restricted=vcat(Ax_restricted,Ax_3)
             Ax_restricted=sort(vec(Ax_restricted))


         elseif(index[1]==3)

             Ax_1=collect(0:1:nelx-1)
             Ax_2=collect(1/3:1:nelx)
             Ax_3=collect(0.5:1:nelx)
             Ax_4=collect(0.8:1:nelx)
             Ax_restricted=vcat(Ax_1,Ax_2)
             Ax_restricted=vcat(Ax_restricted,Ax_3)
             Ax_restricted=vcat(Ax_restricted,Ax_4)
             Ax_restricted=sort(vec(Ax_restricted))


             Ay_1=collect(0:1:nely-1)
             Ay_2=collect(1/3:1:nely)
             Ay_3=collect(0.5:1:nely)
             Ay_4=collect(0.8:1:nely)
             Ay_restricted=vcat(Ay_1,Ay_2)
             Ay_restricted=vcat(Ay_restricted,Ay_3)
             Ay_restricted=vcat(Ay_restricted,Ay_4)
             Ay_restricted=sort(vec(Ay_restricted))
         elseif(index[1]==4)
             Ax_1=collect(0.1:1:nelx)
             Ax_2=collect(0.2:1:nelx)
             Ax_3=collect(0.5:1:nelx)
             Ax_4=collect(0.6:1:nelx)
             Ax_5=collect(0.7:1:nelx)
             Ax_restricted=vcat(Ax_1,Ax_2)
             Ax_restricted=vcat(Ax_restricted,Ax_3)
             Ax_restricted=vcat(Ax_restricted,Ax_4)
             Ax_restricted=vcat(Ax_restricted,Ax_5)
             Ax_restricted=sort(vec(Ax_restricted))


             Ay_1=collect(0.1:1:nely)
             Ay_2=collect(0.2:1:nely)
             Ay_3=collect(0.5:1:nely)
             Ay_4=collect(0.6:1:nely)
             Ay_5=collect(0.7:1:nely)
             Ay_restricted=vcat(Ay_1,Ay_2)
             Ay_restricted=vcat(Ay_restricted,Ay_3)
             Ay_restricted=vcat(Ay_restricted,Ay_4)
             Ay_restricted=vcat(Ay_restricted,Ay_5)
             Ay_restricted=sort(vec(Ay_restricted))
        end



#            println(size(Ax_restricted))
#           println((Ax_restricted))
#           println((Ax))

            Zc=itp(Ax_restricted,Ay_restricted)
    #        if(index[1]==2)
#figure()
#surf(Zf)
#figure()
#surf(Zc)
#sleep(40)
#            end
            Qc = Beam_sample_L_High(Zc,Lx,Ly,he,(index[1]-1),MatlabSampler,folder)
#            if(index[1]==2)

#            sleep(40)
#        end
            dQ += value*Qc
        end
        Error=false
    #    println("return no error")
        return (dQ,Qf,Error)
    catch ex
            println(ex)
            println("Error in return map caught")
            Zf=0
        Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
        Zf=FieldTransformation.Transform_L(Zf)
        Qf = Beam_sample_L_High(Zf,Lx,Ly,he,index[1],MatlabSampler,folder)
        println("Succesfully sampled at initial level")
        dQ = Qf
        println("Start previous sample")
        for (key,value) in diff(index)
            step = (index - key).I .+ 1

            Zc = Zf
            he=Lx/size(Zc,1)
            Qc = Beam_sample_L_High(Zc,Lx,Ly,he,(index[1]-1),MatlabSampler,folder)

            println("Succesfully sampled at previous level")
            dQ += value*Qc
        end
        Error=false
        println("Return map error ")
        return (dQ,Qf,Error)
        end
        (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
        Error=Error
        dQ=dQ
        Qf=Qf

    end
    Error=Error
    dQ=dQ
    Qf=Qf

        #@show Qf

        # compute difference

        return (dQ,Qf)
end


function Beam_NL_Het_High_GP(index::Index, ξ::Vector{T} where {T<:Real}, Lx::Float64,Ly::Float64 , grf::GaussianRandomField,MatlabSampler::Function,folder::String,nelx::Int64,nely::Int64,he::Float64)
    Error=true

    while(Error ==true)
    #    println("Start sampling")
       o(Lx,Ly,index,grf,Error)=try
    #       println("Start try catch")
        Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
        Zf=FieldTransformation.Transform_NL(Zf)
        Qf = Beam_sample_NL_High(Zf,Lx,Ly,he,index[1],MatlabSampler,folder)
        dQ = Qf
        for (key,value) in diff(index)
            step = (index - key).I .+ 1

#            println(index[1]+1)
#            println(index[1])
#            println(nely)
#            println(size(Zf))
            Ax=collect(0:1/(index[1]+2):nelx-1/(index[1]+2))
            Ay=collect(0:1/(index[1]+2):nely-1/(index[1]+2))
#          println(size(Ax))
#          println((Ax))

            knots=(Ax,Ay,)
            itp = Interpolations.interpolate(knots, Zf, Gridded(Linear()))

            if(index[1]==1)
#             Ax_1=collect(0:1:nelx-1)
#             Ax_2=collect(2/3:1:nelx)

#             Ax_1=collect(0.2:1:nelx)
#             Ax_2=collect(0.6:1:nelx)

            Ax_1=collect(0.2236:1:nelx)
             Ax_2=collect(0.6:1:nelx)

             Ax_restricted=vcat(Ax_1,Ax_2)
             Ax_restricted=sort(vec(Ax_restricted))


#             Ay_1=collect(0:1:nely-1)
#             Ay_2=collect(2/3:1:nely)

#             Ay_1=collect(0.2:1:nely)
#             Ay_2=collect(0.6:1:nely)

             Ay_1=collect(0.2236:1:nely)
             Ay_2=collect(0.6:1:nely)

             Ay_restricted=vcat(Ay_1,Ay_2)
             Ay_restricted=sort(vec(Ay_restricted))

         elseif(index[1]==2)

    #         Ay_1=collect(0:1:nely-1)
    #         Ay_2=collect(0.375:1:nely)
    #         Ay_3=collect(0.75:1:nely)

             Ay_1=collect(0.0392:1:nely)
             Ay_2=collect(0.375:1:nely)
             Ay_3=collect(0.71:1:nely)

             Ay_restricted=vcat(Ay_1,Ay_2)
             Ay_restricted=vcat(Ay_restricted,Ay_3)
             Ay_restricted=sort(vec(Ay_restricted))


             Ax_1=collect(0.0392:1:nelx)
             Ax_2=collect(0.375:1:nelx)
             Ax_3=collect(0.71:1:nelx)
             Ax_restricted=vcat(Ax_1,Ax_2)
             Ax_restricted=vcat(Ax_restricted,Ax_3)
             Ax_restricted=sort(vec(Ax_restricted))


         elseif(index[1]==3)

             Ax_1=collect(0:1:nelx-1)
             Ax_2=collect(1/3:1:nelx)
             Ax_3=collect(0.5:1:nelx)
             Ax_4=collect(0.8:1:nelx)
             Ax_restricted=vcat(Ax_1,Ax_2)
             Ax_restricted=vcat(Ax_restricted,Ax_3)
             Ax_restricted=vcat(Ax_restricted,Ax_4)
             Ax_restricted=sort(vec(Ax_restricted))


             Ay_1=collect(0:1:nely-1)
             Ay_2=collect(1/3:1:nely)
             Ay_3=collect(0.5:1:nely)
             Ay_4=collect(0.8:1:nely)
             Ay_restricted=vcat(Ay_1,Ay_2)
             Ay_restricted=vcat(Ay_restricted,Ay_3)
             Ay_restricted=vcat(Ay_restricted,Ay_4)
             Ay_restricted=sort(vec(Ay_restricted))
         elseif(index[1]==4)
             println("Not Written")

        end



#            println(size(Ax_restricted))
#           println((Ax_restricted))
#           println((Ax))

            Zc=itp(Ax_restricted,Ay_restricted)
    #        if(index[1]==2)
#figure()
#surf(Zf)
#figure()
#surf(Zc)
#sleep(40)
#            end
            Qc = Beam_sample_NL_High(Zc,Lx,Ly,he,(index[1]-1),MatlabSampler,folder)
#            if(index[1]==2)

#            sleep(40)
#        end
            dQ += value*Qc
        end
        Error=false
    #    println("return no error")
        return (dQ,Qf,Error)
    catch ex
            println(ex)
            println("Error in return map caught")
            Zf=0
        Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
        Zf=FieldTransformation.Transform_NL(Zf)
        Qf = Beam_sample_NL_High(Zf,Lx,Ly,he,index[1],MatlabSampler,folder)
        println("Succesfully sampled at initial level")
        dQ = Qf
        println("Start previous sample")
        for (key,value) in diff(index)
            step = (index - key).I .+ 1

            Zc = Zf
            he=Lx/size(Zc,1)
            Qc = Beam_sample_NL_High(Zc,Lx,Ly,he,(index[1]-1),MatlabSampler,folder)

            println("Succesfully sampled at previous level")
            dQ += value*Qc
        end
        Error=false
        println("Return map error ")
        return (dQ,Qf,Error)
        end
        (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
        Error=Error
        dQ=dQ
        Qf=Qf

    end
    Error=Error
    dQ=dQ
    Qf=Qf

        #@show Qf

        # compute difference

        return (dQ,Qf)
end



function interpolate_field(pts_fine,pts_coarse,Z::Matrix{T}) where {T<:Real}
    itp = Interpolations.interpolate(pts_fine, Z, Gridded(Linear()))
    itp[pts_coarse[1],pts_coarse[2]]
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


end
