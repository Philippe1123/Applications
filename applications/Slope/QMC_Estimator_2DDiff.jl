using DelimitedFiles

using LinearAlgebra

using GaussianRandomFields

using DigitalNets, LatticeRules

#using Plots

using PyPlot

using Printf

using Statistics

using LaTeXStrings

using SpecialFunctions

using SimpleMultigrid

using MultilevelEstimators

const shift_file = joinpath(@__DIR__(),"Shifts")
const shiftvec = readdlm(shift_file,Float64)



function Run()
factor=2
dim=1
nshift=8
nsampleStart=1
tolerance=10^-4


npt = 30
h = 1/(npt-1)
# Gaussian random field
mat = Matern(1, 3)
covariance_function = CovarianceFunction(2, mat)
pts = [i/(npt-1) - h/2 for i in 1:npt-1]
grf = GaussianRandomField(covariance_function, KarhunenLoeve(dim), pts,pts)
println(grf)


f(x,nsampletobetaken) = begin
 ω=ones(length(x),1)
	for id=1:length(x)
#	ω[id]=transform(TruncatedNormal(0,1,-1-nsampletobetaken/1024,1+nsampletobetaken/1024),x[id])
	ω[id]=transform(TruncatedNormal(0,1,-1,1),x[id])

	end
	a = exp.(sample(grf, xi=ω))
#	println(ω[1][1])
#	println(a[1][1])
#	println("----")
	u = solve(a)
	u = u[length(u) ÷ 2]
end

varest=10^5
#dig=DigitalNet64(dim)
lat=LatticeRule(dim)
timestart=0;
TimeElasped=[];
VarEstVec=[];
SampleVec=[];
shift=Dict()
for id=1:nshift
#	println(shiftvec[6,id])
#shift[id]=DigitalShiftedDigitalNets64(dig,shiftvec[1:dim,id])
shift[id]=ShiftedLatticeRule(lat,shiftvec[1:dim,id])
end
nsampletobetaken=nsampleStart
nsamplestaken=0
Result=[]
ittstart=nsamplestaken
ittstop=nsampletobetaken
while varest>tolerance


#println("taking ",nsampletobetaken, "samples")
#println(nsampletobetaken)
#println(ittstart)
#println(ittstop-1+ittstart)
sampleMatrix=zeros(nsampletobetaken,nshift,dim)
#println(shift[1][1][1])
idx=1
for i=1:nshift
	idy=1
	for j=ittstart:ittstop-1+ittstart
sampleMatrix[idy,idx,:]=shift[i][j]
idy=idy+1
    end
	idx=idx+1
end

#println(sampleMatrix)
#println(nsamplestaken)
#println(nsampletobetaken)

out=zeros(nsampletobetaken,nshift)
t =@elapsed begin
for i=1:nshift
	for j=1:nsampletobetaken
out[j,i]=f(sampleMatrix[j,i,:],nsampletobetaken)
end
end

end
timestart=timestart+t
append!(TimeElasped,timestart)
append!(SampleVec,nsamplestaken)



nsamplestaken=nsamplestaken+nsampletobetaken
nsampletobetaken=Int64(ceil(nsamplestaken*factor))-nsamplestaken
ittstart=ittstop+ittstart
ittstop=nsampletobetaken
#println(out)

if(isempty(Result))
Result=out
else
Result=[Result ;out]
end

mn=ComputeExpectedVal(Result)
varest=ComputeVarEstimator(Result)
append!(VarEstVec,varest)
println("mean = ", mn, " with var est = ", varest, " with samples = ",nsamplestaken)
#println(varest)


#varest=-1
end
figure()

loglog(SampleVec,VarEstVec,"-o")
vec=[10; 10^3]

loglog(vec,(vec.^(-1/2))./10,label=1/2)
loglog(vec,(vec.^(-1))./10,label=1/2)
loglog(vec,(vec.^(-1.5))./10,label=1/2)
loglog(vec,vec.^(-2),label=1/2)
legend(("sol","O(n^-1/2)","O(n^-1)","O(n^-1.5)","O(n^-2)"))


figure()
loglog(VarEstVec,TimeElasped,"-o")
vec=[10^-5; 10^-3]
loglog(vec,vec.^(-2),label=1/2)
loglog(vec,(vec.^(-1.5))./10,label=1/2)
loglog(vec,(vec.^(-1))./10,label=1/2)
loglog(vec,(vec.^(-1/2))./10,label=1/2)
legend(("sol","O(n^-1/2)","O(n^-1)","O(n^-1.5)","O(n^-2)"))

end


function solve(af)
	Af = elliptic2d(af)
    bf = fill(one(eltype(Af)), size(Af, 1))
    uf = Af\bf
	return uf
end



function ComputeExpectedVal(sol::Array{Float64,2})

	m=mean(mean(sol))

	return m

end


function ComputeSampleVariance(sol::Array{Float64,2})




end


function ComputeVarEstimator(sol::Array{Float64,2})

nshift=size(sol,2)
nsample=size(sol,1)
meanOverSamples=mean(sol,dims=1)
out=std(meanOverSamples)/sqrt(nshift)


return out

end


Run()
