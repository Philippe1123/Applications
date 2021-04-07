#### COMMENT : Shut down Julia , run the piece of code in comment , then execute file
using Pkg
using DigitalNets
using Reporter
using FileIO
using PyPlot
using DelimitedFiles

function Run()
#data_dir = string("/home/philippe/.julia/dev/DigitalNets/Example")
data_dir = string("/home/philippe/Desktop/RunFiles/Results_UNCECOMP/nu2")
filepath = joinpath(data_dir,"Slope QMC_Het_High_01-03-2021-T:09:32:22")
filepath_2 = string(filepath,"/Slope QMC_Het_High_01-03-2021-T:09:32:22.jld2")
history=load(filepath_2,"history")
report(history,filepath,include_preamble=true)
#println("Complete")


num=1
num2=num+1
num3=num2+1


figure()
figure()
#nu=2
vecstr=["Slope QMC_Het_High_01-03-2021-T:09:32:55","Slope QMC_Het_High_01-03-2021-T:09:32:22","Slope QMC_Het_High_01-03-2021-T:09:33:30","Slope QMC_Het_High_05-03-2021-T:14:52:20"]
# nu =14
#data_dir = string("/home/philippe/Desktop/RunFiles/Results_UNCECOMP/nu14")
#vecstr=["Slope QMC_Het_High_27-02-2021-T:13:21:41","Slope QMC_Het_High_27-02-2021-T:13:22:46","Slope QMC_Het_High_27-02-2021-T:13:21:11","Slope QMC_Het_High_05-03-2021-T:14:39:43"]

filepath = joinpath(data_dir,vecstr[1])
filepath_2 = string(filepath,"/",vecstr[1],".jld2")
history=load(filepath_2,"history")
meanabs=history.data[end][:mean]
println(meanabs)


for ln = 1:1
    filepath = joinpath(data_dir,vecstr[ln])
    filepath_2 = string(filepath,"/",vecstr[ln],".jld2")
    history=load(filepath_2,"history")

VarEst=zeros(length(history.data),1)
Nsamples=zeros(length(history.data),1)
meandiff=zeros(length(history.data),1)
mean=zeros(length(history.data),1)

for id = 1:length(history.data)
VarEst[id] = history.data[id][:varest]
Nsamples[id]=history.data[id][:nb_of_samples][1]
meandiff[id]=abs.(history.data[id][:mean]-meanabs)
mean[id]=history.data[id][:mean]

end

figure(num)
loglog(Nsamples,VarEst,"-o")
#vec=[10; 10^3]

#loglog(vec,(vec.^(-1/2))./10,label=1/2)
#loglog(vec,(vec.^(-1))./10,label=1/2)
#loglog(vec,(vec.^(-1.5))./10,label=1/2)
#loglog(vec,vec.^(-2),label=1/2)

figure(num2)
loglog(Nsamples,meandiff,"-o")


figure(num3)
semilogx(Nsamples,mean,"-o")

a=[Nsamples mean]
handl=open("/home/philippe/Paper_and_slides/Paper Uncecomp 2021/results/out","w")
writedlm(handl,a)
close(handl)

end
end

Run()
