# Julia code for the book Numerical Linear Algebra

[![Build Status](https://travis-ci.org/EricDarve/numerical_linear_algebra.svg?branch=master)](https://travis-ci.org/EricDarve/numerical_linear_algebra)

To use this code, please download and install Julia from [julialang.org/downloads](https://julialang.org/downloads/).

### Julia installation

Step 1: if you have installed Julia on your computer, you can download this repository using the green `Clone or download` button above. Once you have Julia installed, you can run all the Julia codes in this repository (files with extension `.jl`). 

To run the Julia notebooks (files with extension `.ipynb`), you need to install [IJulia](https://github.com/JuliaLang/IJulia.jl) and [Jupyter](https://jupyter.readthedocs.io/en/latest/). 

Step 2: IJulia is a Julia package. To install it, start the Julia application. At the `julia>` prompt, type:

    using Pkg
    Pkg.add("IJulia")
    
See [IJulia](https://github.com/JuliaLang/IJulia.jl) for more detailed instructions.    

Step 3: Install Jupyter following these [instructions](https://jupyter.readthedocs.io/en/latest/install.html). 

Step 4: type `jupyter notebook` in a [terminal window](https://jupyter.readthedocs.io/en/latest/running.html). Make sure you are in the directory containing your notebook. You should be able to open the notebook from the jupyter window inside your web browser. 

You can run the example notebook `Demo.ipynb` contained in this repository. After opening `Demo.ipynb`, wait for the kernel to be ready (check the top right corner of the window), then click on `Cell -> Run All` to update the plot. You may have to change the Julia kernel to match the one that is installed. Click on `Kernel -> change kernel` for this.

Read the notebook. You should see a plot at the end.

You are now ready to go!
