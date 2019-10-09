# Julia code for the book Numerical Linear Algebra

[![Build Status](https://travis-ci.org/EricDarve/numerical_linear_algebra.svg?branch=master)](https://travis-ci.org/EricDarve/numerical_linear_algebra)

To use this code, please go to [JuliaBox](https://www.juliabox.com), or download and install Julia from [julialang.org/downloads](https://julialang.org/downloads/).

---

Please note a recent [update](https://juliacomputing.com/blog/2019/10/03/october-newsletter.html):

"Free Version of JuliaBox Is Ending - Pricing for Paid JuliaBox Starts at Just $7 per Month for Academic Users: In January, we notified the Julia community that Julia's growth was making free JuliaBox unsustainable, and that we would sunset the free version of JuliaBox. As a result, the free version of JuliaBox will end on Oct 31, 2019. We are grateful to all JuliaBox users - especially our paid users. We encourage you to do the following before Oct 31, 2019:

- If you want to continue using JuliaBox, please sign up for the paid version before October 31st, 2019. For academic users, the cost starts at just $7 per month.
    
- If you do not want to continue using JuliaBox, please download and save your code and datasets no later than October 31st, 2019. Your data and your code may no longer be available after October 31st, 2019."

---

### Local Julia installation

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

### JuliaBox

To use [JuliaBox](https://www.juliabox.com), create an account and log in.

The first step is to download the GitHub repository that contains the code. For this you need to copy your JuliaBox SSH keys to GitHub.

Step 1: on JuliaBox, click on `Settings` (under your login name in the top right). Copy your SSH key. It starts with `ssh-rsa`.

Step 2: go to your GitHub account, click on `Settings` (also under your name in the top right). Then select the tab `SSH and GPG keys`. Click on `New SSH key` and paste your JuliaBox key.

The next step is to download the repository to JuliaBox. 

Step 3: in JuliaBox, click on the button called `Git` (it's under `Dashboard`).

In the form, enter the address

    https://github.com/EricDarve/numerical_linear_algebra.git

in `Git Clone URL`. Then click the `+` button. This will create a new git repository in your account with all the files.

Step 4: from the dashboard, click on `Launch`. This will start a new Jupyter session. You should now see the new directory "numerical_linear_algebra". If you want to check that everything works, you can click on

    numerical_linear_algebra

then on 

    Demo.ipynb

inside the directory.
    
Wait for the kernel to be ready (check the top right corner of the window), then click on `Cell -> Run All` to update the plot. Read the notebook. You should see a plot at the end.

You are now ready to go!
