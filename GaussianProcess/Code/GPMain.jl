cd("/home/mv/Dropbox/Teaching/AdvBayesLearn/VT2014/GaussianProcesses/Code/")

include("GPlib.jl")

nSim = 10 
xGrid = linspace(-5, 5, 100)
fSim = SimGP(xGrid, "sin", SquaredExpKernel, [0.6,0.5], nSim)
p = plot(x = xGrid, y = fSim[:,1], Guide.XLabel("x"), Guide.YLabel("f(x)"), Geom.line)
draw(PNG("SimulatedGPs.png", 600px, 400px), p)
