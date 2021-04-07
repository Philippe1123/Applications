using GaussianRandomFields, Plots
using Random, Printf
pyplot()
mat = Matern(0.5,2.0)
cov = CovarianceFunction(2,mat)



pts = 0:0.02:1
grf = GaussianRandomField(cov,CirculantEmbedding(),pts,pts,minpadding=250)
p1=contour(sample(grf), fill=(true, :balance), colorbar=false, aspect_ratio=1, size=(1200,1200))
#gr()
display(p1)
v = grf.data[1]
ev = sort(vec(v), rev=true).^2 .* length(v)
p2=plot(1:length(ev), ev,xaxis=:log, yaxis=:log,reuse=false)
display(p2)
