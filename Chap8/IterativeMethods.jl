module IterativeMethods

import Printf

function success_message(i)
    @Printf.printf "Converged after %d iterations\n" i
end

function failure_message(i)
    @Printf.printf "Failure to converge after %d iterations\n" i
end

function print_residual(i, maxiter, res, nb)
    @Printf.printf "Iteration %3d/%d, " i maxiter
    @Printf.printf "residual %1.2e; |r|/|b| %1.2e\n" res res / nb
end

include("cg.jl")
include("minres.jl")
include("cr.jl")
include("qmr.jl")
include("bicgstab.jl")
include("gmres.jl")

include("jacobi.jl")
include("gs.jl")
include("sor.jl")

export cg, minres, cg, qmr, bicgstab, gmres, jacobi, gs, sor

end
