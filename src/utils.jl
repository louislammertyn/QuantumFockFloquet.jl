triv(t::Float64) = one(ComplexF64)
function all_frequencies(H::Union{FourierFockOperator, FourierFockMatrix})::Vector{Int}
    freqs = Set{Int}()
    for term in H.terms
        for n in term.freqs
            push!(freqs, n)
        end
    end
    return sort!(collect(freqs))
end


