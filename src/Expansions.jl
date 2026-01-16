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
            n == m && continue
            C = get!(Commutators, (n,m), similar(tmp))   # create if not exists
            axpy!(h1.coeffs[i] * h2.coeffs[j], tmp, C)
        end
    end
    return Commutators
end


function FourierMagnus_step_2nd!(
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
    end
    # commutators
    mul!(tmp, dU, U)
    U .= tmp
end

function FourierMagnus_U(H::FourierFockMatrix, t_i::Float64, t_e::Float64, tol::Float64=1e-9)
    @assert t_e > t_i "Make sure the times are in proper order."

    nTs = div(t_e-t_i, H.T)
    micro_motion_t = rem(t_e-t_i, H.T)

    if nTs > 0 && micro_motion_t > tol
        U_Floquet = FourierMagnus_U(H_fourier, t_i, t_i + nTs * H.T, tol)
        t_i_ = t_i + nTs * H.T 
        U_ = FourierMagnus_U(H_fourier, t_i_, t_e, tol)
        return U_ * U_Floquet^nTs

    elseif nTs > 0 && micro_motion_t < tol 
        U_Floquet = FourierMagnus_U(H_fourier, t_i, t_i + nTs * H.T, tol)
        return U_Floquet^nTs

    else
    end

end

function FourierMagnus_U(H::PeriodicFockMatrix, t_i::Float64, t_e::Float64, tol::Float64=1e-9)
    @assert t_e > t_i "Make sure the times are in proper order."

    H = Fourier_op(H)
    nTs = div(t_e-t_i, H.T)
    micro_motion_t = rem(t_e-t_i, H.T)

    if nTs > 0 && micro_motion_t > tol
        U_Floquet = FourierMagnus_U(H_fourier, t_i, t_i + H.T, tol)
        t_i_ = t_i + nTs * H.T 
        U_ = FourierMagnus_U(H_fourier, t_i_, t_e, tol)
        return U_ * U_Floquet^nTs

    elseif nTs > 1 && micro_motion_t < tol 
        U_Floquet = FourierMagnus_U(H_fourier, t_i, t_i + H.T , tol)
        return U_Floquet^nTs

    else
        #initialise Magnus 

        U = similar(first(H.terms)) 
        for i in axes(U, 1)
            U[i,i] = one(ComplexF64)
        end

        dU = similar(U) 
        t = t_i 
        dt = 0.01 * H.T
        while t < t_e
        

    end

end
function VanVleck(H::PeriodicFockOperator, order::Int64)
end