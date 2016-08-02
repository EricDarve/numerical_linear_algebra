Pkg.add("PyPlot")
# Pkg.add("PlotlyJS")
# Pkg.add("PGFPlots")
Pkg.add("Plots")

# Pkg.update()

using Plots

pyplot()
plot(rand(10)) 

# plotlyjs()
# plot(rand(10)) 

# pgfplots()
# plot(rand(10)) 

println("Installation completed successfully!")
