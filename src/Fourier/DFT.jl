
function getDFTfreqrsp(h)
    N_samples = length(h)

    ω_set_fft = LinRange(0, 2*π-2*π/N_samples, N_samples)
    DFT_evals = fft(h)

    return ω_set_fft, DFT_evals
end

# function getfreqrsp(h::Vector{Float64}, resolution_multiple::Int = 20)
#     N_samples = length(h)
#
#     ω_set_fft = collect( LinRange(0,2*π-2*π/N_samples,N_samples))
#     DFT_evals = fft(h)
#
#     ω_set = collect( LinRange(0,2*π,N_samples*resolution_multiple) )
#     DTFT_evals = collect( computeDTFTviaformula(h,ω_set[i]) for i = 1:length(ω_set) )
#
#     return ω_set_fft, DFT_evals, ω_set, DTFT_evals
# end


function downsamplewithoutprefiltering( x::Vector{T},
                                        downsample_factor::Int )::Vector{T} where T
    #
    N = length(x)
    M = round(Int, ceil(N/downsample_factor))

    out = Vector{T}(undef, M)
    for m = 1:M

        out[m] = x[(m-1)*downsample_factor+1]
    end

    return out
end

"""
f is time-series in some unit of time; e.g., seconds.
x its corresponding set of time stamps.
u is in cycles per unit time; e.g., Hz.
"""
# cosine modulate a time-series by u ∈ ℝ, a frequency in cycles per unit time.
function modulatediscretesignal( f::Vector{T},
                                x,
                                u::T) where T <: Real
    N = length(x)
    @assert length(f) == N

    out = Vector{T}(undef, N)
    #out = Vector{Complex{T}}(undef, N)
    for n = 1:N
        out[n] = f[n]*cos( 2*π*u*x[n] )
        #out[n] = f[n]*exp( im*2*π*u -x[n] )
    end

    return out
end

# x is filter.
function modulatediscretesignalcomplex( f::Vector{Complex{T}},
                                x,
                                u::T)::Vector{Complex{T}} where T <: Real
    N = length(x)
    @assert length(f) == N

    out = Vector{Complex}(undef, N)
    #out = Vector{Complex{T}}(undef, N)
    for n = 1:N
        #out[n] = f[n]*cos( 2*π*u*x[n] )
        out[n] = f[n]*exp( im*2*π*u*x[n] )
    end

    return out
end

function performDFTanalysis(h::Vector{T}) where T

    N_samples = length(h)

    # relative frequency indices for fft().
    ω_set_fft = collect( LinRange(0, 2*π - 2*π/N_samples, N_samples))

    DFT_h = fft(h)
    mag_rsp_fft = abs.(DFT_h)
    phase_rsp_fft = angle.(DFT_h)

    return return ω_set_fft, DFT_h, mag_rsp_fft, phase_rsp_fft
end

function plotDFTmagnitudersp(   x_fft::Vector{Complex{T}},
                                fig_num::Int;
                                title_string::String = "DFT magnitude response") where T <: Real
    #
    N = length(x_fft)
    ω_set_fft = collect( LinRange(0, 2*π - 2*π/N, N))

    PyPlot.figure(fig_num)
    fig_num += 1

    PyPlot.plot( ω_set_fft, abs.(x_fft))
    PyPlot.title(title_string)

    return fig_num
end

function rad2freq(r,fs)
    fn = fs/2
    #r/π = x/fn, x is in Hz.
    return fn*r/π
end

function freq2rad(x,fs)
    fn = fs/2

    return x*π/fn
end
