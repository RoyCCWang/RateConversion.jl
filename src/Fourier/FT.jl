
# Uses hquadrature.
# f is frequency in Hz.
# g: ℝ → ℝ, g(t) = 0 if t ∉ [a,b].
function evalFT(  g::Function,
                    a::T,
                    b::T,
                    f::T;
                    max_integral_evals::Int = 10000,
                    initial_divisions::Int = 1 )::Tuple{Complex{T},T} where T

    return HCubature.hquadrature(tt->g(tt)*exp(-im*2*pi*tt*f),
                                        a,
                                        b;
                                        norm = norm,
                                        atol = 0,
                                        maxevals = max_integral_evals,
                                        initdiv = initial_divisions)

end

function evalFTbatch(U, FT_g::Function, dummy_val::T) where T <: Real
    N = length(U)

    FT_evals = Vector{Complex{T}}(undef, N)
    NI_errors = Vector{T}(undef, N)

    for n = 1:length(U)
        FT_evals[n], NI_errors[n] = FT_g(U[n])
    end

    return FT_evals, NI_errors
end

#### observation model: single component

# first return is from equation 11.46 of Spin Dynamics.
# second return is from equation 5.4 of Spin Dynamics.
function evalFIDcomponent(t::T, α, β, λ, Ω)::Complex{T} where T <: Real

    if t < zero(T)
        return zero(T)
    end

    term1 = im*(Ω*t+β)
    term2 = -λ*t

    #return im*α*exp(term1+term2)*exp(im*β)
    return α*exp(term1+term2)
end
