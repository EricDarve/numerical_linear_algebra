Pkg.add("HDF5")

Pkg.add("PyPlot")
Pkg.add("PlotlyJS")
Pkg.add("PGFPlots")
Pkg.add("Plots")

# Pkg.update()

using HDF5
n = 4
A = rand(n, n)
h5write("test.h5", "A", A)
A0 = h5read("test.h5", "A")
@assert A == A0

using Plots

pyplot()
plot(rand(10)) 

plotlyjs()
plot(rand(10)) 

pgfplots()
plot(rand(10)) 

println("Installation completed successfully!")
