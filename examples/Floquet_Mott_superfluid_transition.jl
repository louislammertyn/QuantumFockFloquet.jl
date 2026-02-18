# This example reproduces partial results from the paper https://doi.org/10.1103/PhysRevLett.95.260404 
# It shows how periodic driving if the Bose-Hubbard Hamiltonian can interpolate between the Mott insulating phase and the superfluid phase 

using Revise
using FoSpFloquet
using LinearAlgebra

# Define space 
L = 5
N = 5
cutoff = N

geometry = (L,)
V = U1FockSpace(geometry, cutoff, N)
lattice = Lattice(geometry)
states = all_states_U1(V)

# define the Hamiltonian
J, U = 1., 3.
H_j, H_u = Bose_Hubbard_H(V, lattice, J, U)

## define the drive 
# params

V_drive = ZeroFockOperator()
for i in 1:L 
    V_drive +=  i * ni(V, i)
end



ω = 14.
K = 2.4 * ω
T = 2π/ω
f_t(t) = K * cos(ω * t )

# convert to relevant Operators
H_t = PeriodicFockOperator([H_j+ H_u, V_drive], [triv, f_t], T)
H_fourier = Fourier_op(H_t)
H_fourier_m = matrix_rep(H_fourier, states)

# Storage
U1s = []
U2s = []

# Tolerances
tols = [1/10^i for i in 3:4]

for tol in tols
    U1, _ = compute_Floquet(H_fourier_m, 0., 1; tol=tol)
    push!(U1s, U1)
    
    U2, _ = compute_Floquet(H_fourier_m, 0., 2; tol=tol)
    push!(U2s, U2)
end
U2s[end]' * U2s[end]
# Check convergence by differences
println("Convergence for U1:")
for i in 2:length(U1s)
    println("ΔU1 (tol $(tols[i])) = ", norm(U1s[i] - U1s[i-1]))
end

println("Convergence for U2:")
for i in 2:length(U2s)
    println("ΔU2 (tol $(tols[i])) = ", norm(U2s[i] - U2s[i-1]))
end

# Identity matrix
I1 = Matrix{ComplexF64}(I, size(U1s[1])...)
I2 = Matrix{ComplexF64}(I, size(U2s[1])...)

println("Convergence for U1 (unitary norm):")
for i in 2:length(U1s)
    err_unitary = norm(U1s[i]' * U1s[i-1] - I1, Inf)
    println("tol=$(tols[i]): err_unitary = $err_unitary")
end

println("Convergence for U2 (unitary norm):")
for i in 2:length(U2s)
    err_unitary = norm(U2s[i]' * U2s[i-1] - I2, Inf)
    println("tol=$(tols[i]): err_unitary = $err_unitary")
end
H_fourier_m(0.) - H_fourier_m(0.)'
U' 