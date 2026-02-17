module FoSpFloquet

using Reexport
#using QuantumFockCore
@reexport using FoSpDynamics
using LinearAlgebra
using FFTW

include("./utils.jl")
include("./PeriodicFockOp.jl")
include("./Expansions.jl")


export AbstractTDFockOperator,
       PeriodicFockOperator,
       PeriodicFockMatrix,
       FourierTerm,
       FourierFockOperator,
       FourierFockMatrix,
       Fourier_op,
       matrix_rep
       
export all_frequencies


export analytic_integral,
       analytic_double_integral,
       all_commutators_Fourier,
       FM_step_2nd!,
       FM_U,
       compute_Floquet,
       FMTE, FMTE_1st, FMTE_2nd

export PieceWiseOperator, PieceWiseMatrix
end
