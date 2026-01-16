struct ConstantFunction
    value::ComplexF64  # the constant value
end
# Call syntax: f(t) = value
(c::ConstantFunction)(t) = c.value

Base.:==(c1::ConstantFunction, c2::ConstantFunction) = c1.value == c2.value
Base.show(io::IO, c::ConstantFunction) = print(io, "ConstantFunction(", c.value, ")")


function all_frequencies(H::Union{FourierFockOperator, FourierFockMatrix})
    freqs = Set{Int}()
    for term in H.terms
        for n in term.freqs
            push!(freqs, n)
        end
    end
    return sort!(collect(freqs))
end