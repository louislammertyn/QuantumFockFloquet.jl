abstract type AbstractTDFockOperator end

struct PeriodicFockOperator <: AbstractTDFockOperator
    terms::Dict{Function, MultipleFockOperator}
    T::Float64
end

struct FourierFockOperator <: AbstractTDFockOperator
    terms::Dict{Int64, MultipleFockOperator}
    ω::Float64
end

function fourier_Op(op::PeriodicFockOperator, N::Int=100)
    T = op.T
    tgrid = range(0, T, length=N+1)[1:end-1]

    f_terms = dict{Int64, AbstractFockOperator}()
    freqs = div(-N,2):div(N,2)
    for i in freqs
        f_terms[i] = ZeroFockOperator()
    end
    
    t_samples = range(0, T, 10)
    
    for (f, O) in op
        if all(t -> f(t)==1, t_samples)
            f_terms[0] += O
            continue
        end

        fourier_f = fftshift(fft(f.(tgrid))) / N  

        maxf = maximum(fourier_f)

        for (i, f) in enumerate(fourier_f)
            (f < maxf*.01) && (f_terms[freqs[i]] += O)
        end
    end
    
    return FourierFockOperator(f_terms, 2*π/T)

end


