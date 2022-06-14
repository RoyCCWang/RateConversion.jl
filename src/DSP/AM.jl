# prepare analytic amplitude
function prepareAs(s_t::Vector{Complex{T}},
                    fs::T,
                    ν_st::T,
                    ν_fin::T,
                    bp_δ::T,
                    N_h) where T

    ### get bandpass filter.
    ch_BW = ν_fin-ν_st
    ν0 = (ν_st+ν_fin)/2

    bp_a = ν0 - ch_BW/2 #- bp_δ
    bp_b = ν0 + ch_BW/2 #+ bp_δ

    h_func = xx->cosinetransitionbandpassimpulsefunc(xx, bp_a, bp_b, bp_δ, fs)
    t_h = gettimerangetunablefilter(N_h, fs)
    h = h_func.(t_h)

    ### normalize for improved numerical conditioning.
    #   and set up DTFTs.
    t = gettimerange(length(s_t), fs)

    r = h
    @assert isodd(length(r))

    ### apply bandpass filter.
    rs = fastlinearconv(s_t, r)

    M = round(Int, (length(r)-1)/2) # for odd length(r)
    N = length(s_t)
    t_rs_idx = -M:1:(M+N-1)
    Ts = one(T)/fs
    t_rs = Ts .* t_rs_idx

    ### get analytic signal.
    ANx_t_rs = DSP.hilbert(real.(rs))
    ANy_t_rs = DSP.hilbert(imag.(rs))

    # truncate range to match s_t
    idx_range = M+1:M+1+N
    Ax_t = abs.(ANx_t_rs[idx_range])
    Ay_t = abs.(ANy_t_rs[idx_range])
    t = t_rs[idx_range]

    return Ax_t, Ay_t, t, rs[idx_range]
end


function prepareAs2(s_t::Vector{Complex{T}},
                    bp_a::T,
                    bp_b::T,
                    bp_δ::T,
                    fs::T,
                    N_h::Int) where T

    ### get bandpass filter.

    h_func = xx->cosinetransitionbandpassimpulsefunc(xx, bp_a, bp_b, bp_δ, fs)
    t_h = gettimerangetunablefilter(N_h, fs)
    h = h_func.(t_h)

    ### normalize for improved numerical conditioning.
    #   and set up DTFTs.
    t = gettimerange(length(s_t), fs)

    r = h
    @assert isodd(length(r))

    ### apply bandpass filter.
    rs = fastlinearconv(s_t, r)

    M = round(Int, (length(r)-1)/2) # for odd length(r)
    N = length(s_t)
    t_rs_idx = -M:1:(M+N-1)
    Ts = one(T)/fs
    t_rs = Ts .* t_rs_idx

    ### get analytic signal.
    ANx_t_rs = DSP.hilbert(real.(rs))
    ANy_t_rs = DSP.hilbert(imag.(rs))

    # truncate range to match s_t
    idx_range = M+1:M+1+N
    Ax_t = abs.(ANx_t_rs[idx_range])
    Ay_t = abs.(ANy_t_rs[idx_range])
    t = t_rs[idx_range]

    return Ax_t, Ay_t, t, rs[idx_range], rs, ANx_t_rs, ANy_t_rs
end
