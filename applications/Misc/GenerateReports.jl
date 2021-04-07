#### COMMENT : Shut down Julia , run the piece of code in comment , then execute file
using Pkg
using DigitalNets
using Reporter
using FileIO
using PyPlot

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

#vecstr=["Simple_Diffusion_Problem15-02-2021-T:10:21:12Lattice","Simple_Diffusion_Problem15-02-2021-T:10:23:46Sobol","Simple_Diffusion_Problem15-02-2021-T:10:25:04Sobol2","Simple_Diffusion_Problem15-02-2021-T:10:25:12Sobol3"]
#vecstr=["Simple_Diffusion_Problem15-02-2021-T:11:23:13Lattice","Simple_Diffusion_Problem15-02-2021-T:11:38:10Sobol","Simple_Diffusion_Problem15-02-2021-T:11:54:30Sobol2","Simple_Diffusion_Problem15-02-2021-T:12:01:28Sobol3"]

#vecstr=["Beam QMC_Het_Lin_High_GP15-02-2021-T:15:14:49","Beam QMC_Het_Lin_High_GP15-02-2021-T:15:17:10","Beam QMC_Het_Lin_High_GP15-02-2021-T:15:12:37","Beam QMC_Het_Lin_High_GP15-02-2021-T:15:08:30"]
#legend(("Sobol","Lattice","Sobol3","Sobol2"))


#vecstr=["Beam QMC_Het_Lin_High_GP15-02-2021-T:19:21:22","Beam QMC_Het_Lin_High_GP15-02-2021-T:19:21:29","Beam QMC_Het_Lin_High_GP15-02-2021-T:19:21:03","Beam QMC_Het_Lin_High_GP15-02-2021-T:19:15:26"]
#vecstr=["Beam QMC_Het_Lin_High_GP15-02-2021-T:22:48:53","Beam QMC_Het_Lin_High_GP15-02-2021-T:22:48:56","Beam QMC_Het_Lin_High_GP15-02-2021-T:22:49:00","Beam QMC_Het_Lin_High_GP15-02-2021-T:22:48:25"]
#vecstr=["Beam QMC_Het_Lin_High_GP16-02-2021-T:09:34:26","Beam QMC_Het_Lin_High_GP16-02-2021-T:09:34:27","Beam QMC_Het_Lin_High_GP16-02-2021-T:09:33:24","Beam QMC_Het_Lin_High_GP16-02-2021-T:09:33:58"]
#vecstr=["Beam QMC_Het_Lin_High_GP16-02-2021-T:11:56:11","Beam QMC_Het_Lin_High_GP16-02-2021-T:11:56:40","Beam QMC_Het_Lin_High_GP16-02-2021-T:11:55:06","Beam QMC_Het_Lin_High_GP16-02-2021-T:11:55:43"]
#vecstr=["Beam QMC_Het_Lin_High_GP16-02-2021-T:15:20:49","Beam QMC_Het_Lin_High_GP16-02-2021-T:15:20:20","Beam QMC_Het_Lin_High_GP16-02-2021-T:15:17:31","Beam QMC_Het_Lin_High_GP16-02-2021-T:15:19:43"]
#vecstr=["Beam QMC_Het_Lin_High_GP17-02-2021-T:13:57:29","Beam QMC_Het_Lin_High_GP17-02-2021-T:13:58:00","Beam QMC_Het_Lin_High_GP17-02-2021-T:13:54:42","Beam QMC_Het_Lin_High_GP17-02-2021-T:13:55:32"]
#vecstr=["Beam QMC_Het_Lin_High_GP23-02-2021-T:21:25:55","Beam QMC_Het_Lin_High_GP23-02-2021-T:21:27:39","Beam QMC_Het_Lin_High_GP23-02-2021-T:21:24:32","Beam QMC_Het_Lin_High_GP23-02-2021-T:21:24:54","Beam QMC_Het_Lin_High_GP24-02-2021-T:14:23:01"]
vecstr=["Slope QMC_Het_High_26-02-2021-T:12:12:57","Slope QMC_Het_High_26-02-2021-T:12:12:26","Slope QMC_Het_High_26-02-2021-T:12:13:29" ]
#vecstr=["Slope QMC_Het_High_26-02-2021-T:12:11:01","Slope QMC_Het_High_26-02-2021-T:12:05:44","Slope QMC_Het_High_26-02-2021-T:12:11:53" ]
vecstr=["Slope QMC_Het_High_27-02-2021-T:13:21:41","Slope QMC_Het_High_27-02-2021-T:13:22:46","Slope QMC_Het_High_27-02-2021-T:13:21:11" ]

for ln = 1:length(vecstr)
    filepath = joinpath(data_dir,vecstr[ln])
    filepath_2 = string(filepath,"/",vecstr[ln],".jld2")
    history=load(filepath_2,"history")

VarEst=zeros(length(history.data),1)
Nsamples=zeros(length(history.data),1)
elapsd=zeros(length(history.data),1)
el=0
for id = 1:length(history.data)
VarEst[id] = history.data[id][:varest]
Nsamples[id]=history.data[id][:nb_of_samples][1]
elapsd[id]=history.data[id][:elapsed]+el
el=elapsd[id]
end
println(VarEst)
println(Nsamples)
figure(num)
loglog(Nsamples,VarEst,"-o")
#vec=[10; 10^3]

#loglog(vec,(vec.^(-1/2))./10,label=1/2)
#loglog(vec,(vec.^(-1))./10,label=1/2)
#loglog(vec,(vec.^(-1.5))./10,label=1/2)
#loglog(vec,vec.^(-2),label=1/2)

figure(num2)
loglog(Nsamples,sqrt.(VarEst),"-o")


figure(num3)
loglog(VarEst,elapsd,"-o")




end
figure(num)
xlabel("Number of samples n")
ylabel("var est")
legend(("Sobol2","Lattice","lattent"))

figure(num2)

xlabel("Number of samples n")
ylabel("std est")
legend(("Sobol2","Lattice","lattent"))
title("Beam order 2 elements with truncated Normal distributed [-2,2]")

vec=[10; 10^3]

loglog(vec,(vec.^(-1/2))./1000)
loglog(vec,(vec.^(-1))./1000)
loglog(vec,(vec.^(-1.5))./250)
loglog(vec,vec.^(-2)./100)

figure(num3)
xlabel("var est")
ylabel("elapsed")
legend(("Sobol2","Lattice","lattent"))

end

Run()
