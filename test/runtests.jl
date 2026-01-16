using Test
using LinearAlgebra
using FFTW
using Revise
using FoSpFloquet


# --------------------------
# Setup
# --------------------------
geometry = (10,)
V = U1FockSpace(geometry, 4, 4)
basis = all_states_U1(V)
lattice = Lattice(geometry)

t = ManyBodyTensor_rnd(ComplexF64, V, 1, 1)
t2 = ManyBodyTensor_rnd(ComplexF64, V, 2, 2, 0.1)

op1 = nbody_Op(V, lattice, t)
op2 = nbody_Op(V, lattice, t2)

# --------------------------
# Tests
# --------------------------
@testset "TDFockOperator Tests" begin
    # -----------------------------
    # Create dummy operators
    # -----------------------------
    T = 1.0
    f1(t) = sin(2π*t/T)
    f2(t) = sin(2π*t/(0.5*T))

    periodic_op = PeriodicFockOperator(Dict(f1=>op1, f2=>op2), T)
    

    # -----------------------------
    # Test matrix representation
    # -----------------------------
    pmatrix = matrix_rep(periodic_op, basis)
    @test isa(pmatrix, PeriodicFockMatrix)
    @test all(isa.(values(pmatrix.terms), Matrix{ComplexF64}))

    # -----------------------------
    # Test Fourier transform
    # -----------------------------
    fourier_op = fourier_Op(periodic_op, 16)
    @test isa(fourier_op, FourierFockOperator)
    println(keys(fourier_op.terms))
    @test all(k ∈ -2:2 for k in keys(fourier_op.terms))
    # -----------------------------
    # Test Fourier of matrix representation
    # -----------------------------
    pmatrix = matrix_rep(periodic_op, basis)
    fourier_matrix_op = fourier_Op(pmatrix, 16)
    @test isa(fourier_matrix_op, FourierFockMatrix)
    @test all(isa.(values(fourier_matrix_op.terms), Matrix{ComplexF64}))
end