# Apply a single midpoint-rule step
function midpoint_step!(U::Matrix{ComplexF64}, t::Float64, dt::Float64, ops_eigen::Vector, dU_cache::Matrix{ComplexF64})
    fill!(dU_cache, 0)
    for i in 1:size(dU_cache,1)
        dU_cache[i,i] = 1.0 + 0im
    end

    for (f, es, vs) in ops_eigen
        f_dt = f(t + dt/2) * dt
        exp_e = exp.(-1im * f_dt .* es)
        dU_cache .= vs * Diagonal(exp_e) * vs'
    end
    U .= dU_cache * U
    return U
end

# Full evolution using midpoint rule
function Ufl_midpoint(H::PeriodicFockMatrix, tol::Float64, ti::Float64=0.0)
    N = length(H.terms)
    first_matrix = first(values(H.terms))
    U = Matrix{ComplexF64}(I, size(first_matrix)...)
    dU_cache = similar(U)

    # Precompute eigendecompositions
    ops_eigen = Vector{Tuple{Function, Vector{Float64}, Matrix{ComplexF64}}}()
    for (f, h) in H.terms
        es, vs = eigen(h)
        push!(ops_eigen, (f, real.(es), vs))
    end

    # Timestep scheme
    te = ti + H.T
    dt = 0.01 * H.T
    t = ti

    u1 = similar(U)
    u2 = similar(U)
    while t < te
        midpoint_step!(U, t, dt, ops_eigen, dU_cache)
        t += dt
    end

    return U
end

function U_strang(H::FourierFockMatrix, t1::Float64, t2::Float64, tol::Float64=1e-8)

end

function Sambe_diagonalisation(H::FourierFockMatrix, M::Int64, basis::Vector{AbstractFockState})
    O_m = matrix_rep(O, basis)
    
    Of_m = fourier_Op(O_m)
    extended_size = (values(O.terms)[1] |> size |> collect ) .* (2*M+1)
    @assert prod(extended_size) > 2e4 "The matrix considered is very large, consider alternatives to full dense diagonalisation"
    K = zeros(ComplexF64, extended_size)
    for (ω, of_m) in Of_m
    end
end
