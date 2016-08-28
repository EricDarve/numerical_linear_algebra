if length(Base.LOAD_CACHE_PATH) >= 3
    splice!(Base.LOAD_CACHE_PATH, 3)
end

println("[install_all_packages.jl] Installing HDF5")

Pkg.add("HDF5")

println("[install_all_packages.jl] Testing HDF5")

using HDF5
n = 4
A = rand(n, n)
h5open("test.h5", "w") do file
    write(file, "A", A)
end
A0 = h5read("test.h5", "A")
@assert A == A0

println("[install_all_packages.jl] Installation of HDF5 completed successfully!")

# This completes the installation of HDF5

println("[install_all_packages.jl] Installing Plots")

Pkg.add("Plots")
Pkg.add("PlotlyJS")

using Plots
plotlyjs()
plot(rand(10)) 
# This completes the installation of Plots

println("[install_all_packages.jl] Installation completed successfully!")
