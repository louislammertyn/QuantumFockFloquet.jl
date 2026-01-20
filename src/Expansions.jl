function analytic_integral(t::Float64, dt::Float64, n::Int, ω::Float64)
    n == 0 && return dt
    return 1im * exp(-1im * n * ω* t) * (-1+exp(-1im*dt*n*ω))/(n*ω)
end

function analytic_double_integral(t::Float64, dt::Float64, n::Int, m::Int, ω::Float64)::ComplexF64
    # handle the simple special cases to avoid divide-by-zero
    if n == 0 && m == 0
        return 0.5 * dt^2
    elseif n == 0 && m != 0
        return -exp(-1im*ω*m*t)*(ω*dt*m*1im + exp(-1im*ω*dt*m) - 1)/(m^2*ω^2)
    elseif m == 0 && n != 0
        return exp(-1im*ω*n*t)*((n*ω*dt*1im + 1)*exp(-1im*n*ω*dt) - 1)/(n^2*ω^2)
    elseif n + m == 0 && n != 0 && m!=0
        # resonant case n + m = 0
        return (-n*ω*dt*1im - exp(-1im*n*ω*dt) + 1)/(n^2*ω^2)
    else
        # general case
        return exp(-1im * ω * t * (n + m)) *
               ((n + m) * exp(-1im * n * ω * dt) - n * exp(-1im * (n + m) * ω * dt) - m) /
               ((n + m) * n * m * ω^2)
    end
end
function all_commutators_Fourier(H::FourierFockMatrix)
    Commutators = Dict{Tuple{Int, Int}, Matrix{ComplexF64}}()
    tmp = similar(first(H.terms).op)
    tmp2 = similar(tmp)
    for h1 in H.terms, h2 in H.terms 
        mul!(tmp, h1.op, h2.op) 
        mul!(tmp2, h2.op, h1.op)
        tmp .-= tmp2
        for (i,n) in enumerate(h1.freqs), (j,m) in enumerate(h2.freqs )
            !(n<m) && continue
            C = get!(Commutators, (n,m), similar(tmp))   # create if not exists
            axpy!(h1.coeffs[i] * h2.coeffs[j], tmp, C)
        end
    end
    return Commutators
end


function FM_step_2nd!(
    H::FourierFockMatrix, frequencies::Vector{Int}, comms::Dict{Tuple{Int, Int}, Matrix{ComplexF64}}, U::Matrix{ComplexF64}, t::Float64, dt::Float64,
    cache::Tuple{
        Matrix{ComplexF64},
        Matrix{ComplexF64},
        Dict{Int, ComplexF64},
        Dict{Tuple{Int, Int}, ComplexF64}
    }
)
    dU, tmp, integrals_single, integrals_double = cache

    # integrals
    for n in frequencies
        integrals_single[n] = analytic_integral(t,dt,n,H.ω)
        for m in frequencies
            integrals_double[(n,m)] = analytic_double_integral(t,dt,n,m,H.ω)
        end
    end


    fill!(dU, zero(ComplexF64))

    # First order sum
    for h in H.terms 
        for (i,n) in enumerate(h.freqs)
            axpy!(integrals_single[n] * h.coeffs[i], h.op, dU)
        end
    end

    # second order sum 
    for ((n,m), M) in comms
        axpy!(integrals_double[(n,m)], M, dU)
        axpy!(-1. * integrals_double[(m,n)], M, dU)
    end
    # commutators
    mul!(tmp, dU, U)
    U .= tmp
end

# Convert PeriodicFockMatrix to FourierFockMatrix
function FM_U(H::PeriodicFockMatrix, t_i::Float64, t_e::Float64, tol::Float64=1e-9)
    H_fourier = Fourier_op(H)
    return FM_U(H_fourier, t_i, t_e, tol)
end

# Convert PeriodicFockMatrix to FourierFockMatrix
function FM_U(H::PeriodicFockMatrix, U_fl::Matrix, t_i::Float64, t_e::Float64, tol::Float64=1e-9)
    H_fourier = Fourier_op(H)
    return FM_U(H_fourier, U_fl, t_i, t_e, tol)
end

# Main function for FourierFockMatrix
function FM_U(H::FourierFockMatrix, t_i::Float64, t_e::Float64, tol::Float64=1e-9)
    @assert t_e > t_i "Make sure the times are in proper order."

    nT = div(t_e - t_i, H.T)          # number of full periods
    t_rem = rem(t_e - t_i, H.T)       # leftover micro-motion

    # case 1: multiple full periods, negligible micro-motion → use Floquet powers 
    if nT > 1 && t_rem < tol
        U_F = compute_Floquet(H, t_i, tol)
        return U_F^nT
    end

    # case 2: one or more periods with leftover micro-motion → First Floquet then leftover micromotion with FMTE
    if nT >= 1
        # Compute Floquet for first period
        U_F = compute_Floquet(H, t_i, tol)
        # Evolve leftover micro-motion
        t_start_rem = t_i + nT * H.T
        U_rem = t_rem > tol ? FMTE(H, t_start_rem, t_e, tol) : Matrix{ComplexF64}(I, size(first(H.terms))...)
        return U_rem * U_F^nT
    end

    # case 3: short interval (<1 period) → just FMTE
    return FMTE(H, t_i, t_e, tol)
end

function FM_U(H::FourierFockMatrix, U_fl::Matrix, t_i::Float64, t_e::Float64, tol::Float64=1e-9)
    @assert t_e > t_i "Make sure the times are in proper order."

    nT = div(t_e - t_i, H.T)          # number of full periods
    t_rem = rem(t_e - t_i, H.T)       # leftover micro-motion

    # case 1: multiple full periods, negligible micro-motion → use Floquet powers 
    if nT > 1 && t_rem < tol
        return U_fl^nT
    end

    # case 2: one or more periods with leftover micro-motion → First Floquet then leftover micromotion with FMTE
    if nT >= 1
        # Evolve leftover micro-motion
        t_start_rem = t_i + nT * H.T
        U_rem = t_rem > tol ? FMTE(H, t_start_rem, t_e, tol) : Matrix{ComplexF64}(I, size(first(H.terms))...)
        return U_rem * U_fl^nT
    end

    # case 3: short interval (<1 period) → just FMTE
    return FMTE(H, t_i, t_e, tol)
end

# Compute Floquet operator in initial time gauge of choice
function compute_Floquet(H::FourierFockMatrix, t_i::Float64, tol::Float64=1e-9)
    t_e = t_i + H.T
    return FMTE(H, t_i, t_e, tol)
end

# Fourier-Magnus time evolution (FMTE) over arbitrary interval 
function FMTE(H::FourierFockMatrix, t_i::Float64, t_e::Float64, tol::Float64=1e-9)
    U = Matrix{ComplexF64}(I, size(first(H.terms))...)  # initialize identity matrix

    # Precompute frequencies and commutators
    freqs = all_freqs(H)
    comms = all_commutators_Fourier(H)

    # Cache for Magnus steps
    dU = similar(U)
    tmp = similar(U)
    ints_single = Dict(k => zero(ComplexF64) for k in freqs)
    ints_double = Dict((n,m) => zero(ComplexF64) for n in freqs, m in freqs)
    cache = (dU, tmp, ints_single, ints_double)

    t = t_i
    dt = 0.01 * H.T  # initial timestep
    u1 = copy(U)
    u2 = copy(U)

    while t < t_e
        # 2nd-order Magnus step (full and two half-steps)
        FM_step_2nd!(H, freqs, comms, u1, t, dt, cache)
        FM_step_2nd!(H, freqs, comms, u2, t, dt/2, cache)
        FM_step_2nd!(H, freqs, comms, u2, t + dt/2, dt/2, cache)

        # Estimate error
        err = norm(u1 - u2) / norm(u1)

        if err < tol
            t += dt
            u1 .= u2
        else
            # Adaptive timestep with machine-precision safeguard
            dt *= 0.9 * (tol / err)^(1/2)
            dt = max(dt, 1e-12)
        end
    end

    return u2
end




function VanVleck(H::PeriodicFockOperator, order::Int64)
end