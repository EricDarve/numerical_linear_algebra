using Printf
using LinearAlgebra
using Random

include("gees.jl")

function test()
    # Householder transforms

    for icase=1:20
        if icase <= 2
            if icase == 1
                n = 2
            else
                n = 128
            end
            x = zeros(n)
            x[1] = 1
        elseif icase <= 4
            if icase == 3
                n = 2
            else
                n = 128
            end
            x = zeros(n)
            x[1] = -1
        elseif icase <= 6
            if icase == 5
                n = 2
            else
                n = 128
            end
            x = zeros(n)
            x[1] = 1
            x[2] = eps(Float64)
        elseif icase <= 8
            if icase == 7
                n = 2
            else
                n = 128
            end
            x = zeros(n)
            x[1] = -1
            x[2] = eps(Float64)
        elseif icase == 9
            n = 3
            x = [244.60002521653865,4.109790267811492e-79,2.786261257333573e-81]
        else
            rng = MersenneTwister(2018)
            n = (icase-7)
            x = rand(rng, n)
        end

        beta, v = house(x)

        # if icase == 1 || icase == 2
        #     @test beta == 0
        # elseif icase == 3 || icase == 4
        #     @test beta == 2
        # end

        # Orthogonal transform tests
        I_n = diagm(0 => fill(1.0,n))
        B = copy(I_n)
        apply_left_householder!(B, 1, 1, beta, v)
        if norm(B'*B - UniformScaling(1.0)) >= (1<<3) * eps(Float64)
            @show v
            @show beta
            @show B
        end
        @test norm(B'*B - UniformScaling(1.0)) < (1<<3) * eps(Float64)

        B = copy(I_n)
        apply_right_householder!(B, n, 1, beta, v)
        @test norm(B'*B - UniformScaling(1.0)) < (1<<3) * eps(Float64)

        c, s = givens(x[1], x[2])

        B = copy(I_n)
        apply_left_givens!(B, 1, 1, c, s)
        @test norm(B'*B - UniformScaling(1.0)) < (1<<3) * eps(Float64)

        B = copy(I_n)
        apply_right_givens!(B, n, 1, c, s)
        @test norm(B'*B - UniformScaling(1.0)) < (1<<3) * eps(Float64)

        # Testing left Householder transform
        rng = MersenneTwister(2018)
        A = rand(rng, n, n)
        A[:,1] = x
        A0 = copy(A)
        beta, v = house(A[:,1])
        apply_left_householder!(A, 1, 1, beta, v)
        @test norm(A[2:n,1]) < (1<<4) * eps(Float64)

        # Testing right Householder transform
        A = copy(A0)'
        beta, v = house(A[1,:]')
        apply_right_householder!(A, n, 1, beta, v)
        @test norm(A[1,2:n]) < (1<<4) * eps(Float64)

        # Givens rotations
        A0 = copy(A)
        c, s = givens(A[1,1], A[2,1])
        apply_left_givens!(A, 1, 1, c, s)
        @test abs(A[2,1]) < 2*eps(Float64)

        # Testing right Givens
        A = copy(A0)'
        c, s = givens(A[1,1], A[1,2])
        apply_right_givens!(A, n, 1, c, s)
        @test abs(A[1,2]) < 2*eps(Float64)
    end

    # Test gehrd!(A)
    rng = MersenneTwister(2018)
    X = rand(rng, n, n)
    D0 = rand(rng,n); D0 = sort(D0,rev=true);
    A = X * diagm(0 => D0) / X
    F = eigen(A); D1 = F.values; V = F.vectors; D1 = sort(D1,rev=true);
    @test norm(D0 - D1) < 1e-11
    gehrd!(A)
    @test norm(tril(A,-2)) < 1e-14
    F = eigen(A); D2 = F.values; V = F.vectors; D2 = sort(D2,rev=true);
    @test norm(D1 - D2) < 1e-11

    # QR with Givens, sequential version
    rng = MersenneTwister(2018)
    n = 8
    X = rand(rng, n, n)
    D0 = rand(rng,n); D0 = sort(D0,rev=true);
    A = X * diagm(0 => D0) / X
    # Upper Hessenberg
    gehrd!(A)
    F = eigen(A); D1 = F.values; V = F.vectors; D1 = sort(D1,rev=true);
    @test norm(D0 - D1) < (1<<5) * eps(Float64)
    @test norm(tril(A,-2)) < (1<<3) * eps(Float64)
    givens_QR_iteration_s!(A)
    F = eigen(A); D1 = F.values; V = F.vectors; D1 = sort(D1,rev=true);
    @test norm(D0 - D1) < (1<<5) * eps(Float64)
    @test norm(tril(A,-2)) < (1<<3) * eps(Float64)

    # QR with bulge chasing
    rng = MersenneTwister(2018)
    n = 8
    X = rand(rng, n, n)
    D0 = rand(rng,n); D0 = sort(D0,rev=true);
    A = X * diagm(0 => D0) / X
    gehrd!(A)
    givens_QR_iteration!(A)
    F = eigen(A); D1 = F.values; V = F.vectors; D1 = sort(D1,rev=true);
    @test norm(D0 - D1) < (1<<5) * eps(Float64)
    @test norm(tril(A,-2)) < (1<<3) * eps(Float64)

    # QR with double shift, single step
    rng = MersenneTwister(2018)
    n = 8
    X = rand(rng, n, n)
    D0 = rand(rng,n); D0 = sort(D0,rev=true);
    A = X * diagm(0 => D0) / X
    gehrd!(A)
    @test norm(tril(A,-2)) < (1<<3) * eps(Float64)
    gees_single_step!(A,false)
    F = eigen(A); D1 = F.values; V = F.vectors; D1 = sort(D1,rev=true);
    @test norm(D0 - D1) < (1<<5) * eps(Float64)
end

test()
gees_testsuite()
