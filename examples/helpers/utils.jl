# N_h should be odd.
function filtermodulatesequence(s_t::Vector{Complex{T}},
    fs::T,
    bp_a::T, bp_b::T;
    bp_δ::T = (bp_b-bp_a)/10,
    N_h::Int = 2001) where T <: Real

    h_func = xx->RateConversion.cosinetransitionbandpassimpulsefunc(xx, bp_a, bp_b, bp_δ, fs)
    t_h = RateConversion.gettimerangetunablefilter(N_h, fs)
    h = h_func.(t_h)
    hs = RateConversion.fastlinearconv(s_t, h)

    M = round(Int, (length(h)-1)/2) # for odd length(r)
    N = length(s_t)
    t_hs_idx = -M:1:(M+N-1)
    Ts = one(T)/fs
    t_hs = Ts .* t_hs_idx

    modulation_shift = -(bp_a-bp_δ) # Hz.
    mhs0 = RateConversion.modulatediscretesignalcomplex(hs, t_hs, modulation_shift)

    return t_hs, mhs0, hs, h
end
