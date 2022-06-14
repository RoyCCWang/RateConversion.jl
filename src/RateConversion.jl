module RateConversion

using FFTW, LinearAlgebra

# Write your package code here.
include("Fourier/DFT.jl")
include("Fourier/DTFT.jl")
include("Fourier/FT.jl")

include("DSP/AM.jl")
include("DSP/cosine_transitions.jl")
include("DSP/DSP_helpers.jl")

end
