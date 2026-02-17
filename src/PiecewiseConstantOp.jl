

struct PieceWiseOperator <: AbstractTDFockOperator
    ops::Vector{MultipleFockOperator}      # store operators
    intervals::Vector{Tuple{Float64,Float64}}   # store intervals (ti,te) corresponding to the relevant op
    T::Float64
end

struct PieceWiseMatrix <: AbstractTDFockOperator
    ops::Vector{Matrix{ComplexF64}}     # store operators
    intervals::Vector{Tuple{Float64,Float64}}    # store intervals (ti,te) corresponding to the relevant op
    T::Float64
end

function PieceWiseOperator(H::PeriodicFockOperator)
    bools = [false for i in 1:length(H.funcs)]
    ti = 0.
    ops = Vector{MultipleFockOperator}()
    intervals=Vector{Tuple{Float64,Float64}}()

    #find first op
    bools = [!iszero(f(0.)) for f in H.funcs]
    @assert count(bools) != 1 "Operator is not piecewise constant"
    
    i = findfirst(bools)
    push!(ops, H.ops[i])

    for t in LinRange(0,H.T, 1001)[1:end-1]
        bools = [!iszero(f(t)) for f in H.funcs]
        # check piecewiseconstant
        @assert count(bools) != 1 "Operator is not piecewise constant"

        new_i = findfirst(bools)

        # if i switches a new interval is entered with a different operator
        if new_i != i 
            te = t 
            push!(intervals, (ti,te))
            ti = te

            i = new_i
            push!(ops, H.ops[i]) 
        end  
    end

    push!(intervals, (ti, H.T))

    return PieceWiseOperator(ops, intervals, H.T)
end

function PieceWiseOperator(H::PeriodicFockMatrix)
    bools = [false for i in 1:length(H.funcs)]
    ti = 0.
    ops = Vector{Matrix{ComplexF64}}()
    intervals=Vector{Tuple{Float64,Float64}}()

    #find first op
    bools = [!iszero(f(0.)) for f in H.funcs]
    @assert count(bools) != 1 "Operator is not piecewise constant"
    
    i = findfirst(bools)
    push!(ops, H.ops[i])

    for t in LinRange(0,H.T, 1001)[1:end-1]
        bools = [!iszero(f(t)) for f in H.funcs]
        # check piecewiseconstant
        @assert count(bools) != 1 "Operator is not piecewise constant"

        new_i = findfirst(bools)

        # if i switches a new interval is entered with a different operator
        if new_i != i 
            te = t 
            push!(intervals, (ti,te))
            ti = te

            i = new_i
            push!(ops, H.ops[i]) 
        end  
    end

    push!(intervals, (ti, H.T))

    return PieceWiseOperator(ops, intervals, H.T)
end
