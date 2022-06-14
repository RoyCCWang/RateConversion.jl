
function myfunc()
    BWs = []

end

function getupsampledsequence(h::Vector{T},
    R::Int) where T

    N = length(h)*R - (R-1) # last R-1 samples are zeros, can ignore.
    h2 = zeros(T, N)
    h2[1:R:end] = h

    return h2
end

# https://dsp.stackexchange.com/questions/37646/filter-order-rule-of-thumb
function getFIRFredHarrislength(fs, BW;
    attenuation_dB = 50.0,
    force_odd_flag = true)::Int

    N_h2 = round(Int, fs/(bp_b-bp_a) *(attenuation_dB/22))
    if iseven(N_h2) && force_odd_flag
        N_h2 += 1
    end

    return N_h2
end

##### lattice-related.

function gettimerange(N::Int, fs::T) where T
   Ts::T = 1/fs

   return zero(T):Ts:(N-1)*Ts
end

function getDFTfreqrange(N::Int, fs::T)::LinRange{T} where T
    a = zero(T)
    b = fs-fs/N

    return LinRange(a, b, N)
end

# this is probably not useful.
function getDFTfreqrangefromoffset(N::Int, fs::T, ν_begin::T)::Tuple{LinRange{T}, Int} where T

    U0 = getDFTfreqrange(N, fs)
    ind = findfirst(xx->(xx>ν_begin), U0)
    @assert typeof(ind) == Int

    if ind > 1 # need to wrap around.
        return LinRange(U0[ind], U0[ind-1] + fs, length(U0)), ind
    end

    return U0, ind
end

function getwraparoundDFTfreqs(N::Int, fs::T, ν_begin::T) where T

    U0 = getDFTfreqrange(N, fs)
    out, inds = wrapfreqrange(U0, ν_begin, fs)

    return U0, out, inds
end



# Assumes U0 is sorted in ascending order.
function wrapfreqrange(U0, ν_begin::T, fs::T) where T <: Real

    N = length(U0)
    ind = findfirst(xx->(xx>ν_begin), U0)

    # TODO handle these exceptions with more grace.
    @assert typeof(ind) == Int
    @assert ind <= N

    out = Vector{T}(undef, N)
    #out[1:ind] = U0[1:ind] .+ fs
    #out[ind+1:end] = U0[ind+1:end]

    M = N-ind
    out[1:M] = U0[ind+1:end]
    out[M+1:end] = U0[1:ind] .+ fs

    inds = collect(1:N)
    inds[1:M] = collect(ind+1:N)
    inds[M+1:end] = collect(1:ind)

    return out, inds
end

# centered around 0; i.e., output at index M/2 + 1 is the 0 second mark.
function gettimerangetunablefilter(N::Int, fs::T) where T
    @assert isodd(N)
    M = N - 1

    Ts::T = 1/fs

    A = fld(M,2)

    return (-A:1:A) * Ts
end

# discard the first and last length(h) samples.
# hs is the filtered output of s, with filter h.
# N_s := length(s), N_h = length(h).
# This way, hs and s are non-zero on the same support.
function matchfilteredoutputwithoriginal(   hs::Vector{T},
                                            N_s,
                                            N_h)::Vector{T} where T
    @assert isodd(N_h)

    M = round(Int, (N_h-1)/2) # for odd length(r)
    idx_range = M+1:M+1+N_s

    return hs[idx_range]
end

##### convolution-related.

# ignore entries of output that require out-of-bound h.
function fastlinearconv(a::Vector{T}, b)::Vector{T} where T
    m = length(a)
    n = length(b)

    if m < 1 || n < 1
        return copy(a)
    end

    out = zeros(T, m+n-1)
    #@inbounds for j = 1:m
    for j = 1:m
        for k = 1:n
            out[j+k-1] += a[j]*b[k]
        end
    end
    return out
end

# setup.
function setupconvDSindices(s, h, R::Int)::Vector{Tuple{Int,Int}}
    m = length(s)
    n = length(h)

    out_range = 1:R:(m+n-1)

    inds = Vector{Tuple{Int,Int}}(undef, length(s)*length(h))
    l = 0
    for j = 1:m
        for k = 1:n

            if (j+k-1) in out_range
                l += 1
                inds[l] = (j,k)
            end
        end
    end

    resize!(inds, l)

    return inds
end

# evaluate at selected positions. R is downsampling factor.
function convDS(s::Vector{T}, h, conv_rang, R::Int)::Vector{T} where T
    m = length(s)
    n = length(h)

    if m < 1 || n < 1
        return copy(a)
    end

    out = zeros(T, m+n-1)
    for ind in conv_rang
        j, k = ind

        out[j+k-1] += s[j]*h[k]
    end

    return out[1:R:end]
end

function applybandpassfilter(s_t::Vector{Complex{T}},
                             stop_band_lower::T,
                             stop_band_upper::T,
                             bp_δ::T,
                             fs::T,
                             N_h::Int) where T <: Real

    ### get bandpass filter.
    @assert isodd(N_h)

    h_func = xx->cosinetransitionbandpassimpulsefunc(xx,
                    stop_band_lower, stop_band_upper, bp_δ, fs)
    t_h = gettimerangetunablefilter(N_h, fs)
    h = h_func.(t_h)

    ### normalize for improved numerical conditioning.
    #   and set up DTFTs.
    t = gettimerange(length(s_t), fs)

    ### apply bandpass filter.
    hs = fastlinearconv(s_t, h)

    return hs, h
end
