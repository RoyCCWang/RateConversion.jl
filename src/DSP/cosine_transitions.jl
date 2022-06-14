
function evalidealnotchimpulsefunc(x::T, a, b, fs)::Complex{T} where T <: Real

    z = im*2*π*x
    term2 = exp(z*a) - one(T) + exp(z*fs) - exp(z*b)

    if zero(T) - eps(T)*2 < x < zero(T) + eps(T)*2
        # singularity case: when x is 0.
        term1 = a*exp(z*a) + fs*exp(z*fs) - b*exp(z*b)

        return (im*2*π*term1 - term2)/(2*π*fs*im)
    end

    return term2/(z*fs)
end


"""
Returns the impulse response of a band-pass filter with cosine transition roll-offs.
The 0 Hz and fs Hz are also passed with a cosine transition roll-off.

The pass band is [a-δ,b+δ]. Units are in Hz.
"""
function cosinetransitionbandpassimpulsefunc(x::T, a, b, δ, fs)::Complex{T} where T <: Real

    term1 = evalDTFTintegralrising(x, a-δ, a+δ, a-δ, a+δ)
    term2 = evalDTFTintegralrect(x, a+δ, b-δ)
    term3 = evalDTFTintegralfalling(x, b-δ, b+δ, b-δ, b+δ)

    return (term1 + term2 + term3)/fs
end

"""
Returns the impulse response of a band-reject filter with cosine transition roll-offs.
The 0 Hz and fs Hz are also rejected with a cosine transition roll-off.

Rejection bands are [a-δ,b+δ] and [0,ϵ] and [fs-ϵ,fs]. All units in Hz.
"""
function cosinetransitionnotchimpulsefunc(x::T, a, b, δ, ϵ, fs)::Complex{T} where T <: Real

    term1 = evalDTFTintegralrising(x, zero(T), ϵ, zero(T), ϵ)
    term2 = evalDTFTintegralrect(x, ϵ, a-δ)
    term3 = evalDTFTintegralfalling(x, a-δ, a+δ, a-δ, a+δ)

    term4 = evalDTFTintegralrising(x, b-δ, b+δ, b-δ, b+δ)
    term5 = evalDTFTintegralrect(x, b+δ, fs-ϵ)
    term6 = evalDTFTintegralfalling(x, fs-ϵ, fs, fs-ϵ, fs)

    return (term1 + term2 + term3 + term4 + term5 + term6)/fs
end

"""
Evaluates the DTFT frequency response of an infinitely-long (non-truncated in time domain)
cosine transition roll-off band-reject filter.
i.e., the ideal DTFT response from which cosinetransitionnotchimpulsefunc() is based on.
"""
function evaltargetcosinenotchfunc(u::T, a, b, δ, ϵ, fs)::T where T <: Real

    if zero(T) <= u < ϵ
        return evalrising(u, zero(T), ϵ)

    elseif ϵ <= u < a-δ
        return one(T)

    elseif a-δ <= u < a+δ
        return evalfalling(u, a-δ, a+δ)

    elseif b-δ <= u < b+δ
        return evalrising(u, b-δ, b+δ)

    elseif b+δ <= u < fs-ϵ
        return one(T)

    elseif fs-ϵ <= u < fs
        return evalfalling(u, fs-ϵ, fs)
    end

    return zero(T)
end

function evalrising(u::T, stop_freq, pass_freq)::T where T <: Real
    return (cos( (pass_freq-u)/(pass_freq-stop_freq)*π )+1)/2
end

function evalfalling(u::T, pass_freq, stop_freq)::T where T <: Real
    return (cos( (pass_freq-u)/(stop_freq-pass_freq)*π )+1)/2
end

function evalDTFTintegralrising(x::T,
                stop_freq,
                pass_freq,
                integration_lower_limit,
                integration_upper_limit) where T <: Real
    # rename.
    a = integration_lower_limit
    b = integration_upper_limit
    p = pass_freq
    s = stop_freq

    # singularity locations.
    x_0_negative = -1/(2*(p-s))
    x_0_positive = 1/(2*(p-s))

    # common intermediate value.
    k2 = 2*π*x*im

    ### term with negative exponent.
    A_negative::Complex{T} = zero(Complex{T})

    if x_0_negative - eps(T)*2 < x < x_0_negative + eps(T)*2
        # singularity case: when x is x_0_negative.
        A_negative = exp(-p*π*im/(p-s))*(b-a)

    else
        k1 = -π*im
        B = p*k1/(p-s)
        c = -k1/(p-s)
        q = c+k2

        A_negative = exp(B)/q * ( exp(b*q) - exp(a*q))
    end

    ### term with positive exponent.
    A_positive::Complex{T} = zero(Complex{T})

    if x_0_positive - eps(T)*2 < x < x_0_positive + eps(T)*2
        # singularity case: when x is x_0_positive.
        A_positive = exp(p*π*im/(p-s))*(b-a)

    else
        k1 = π*im
        B = p*k1/(p-s)
        c = -k1/(p-s)
        q = c+k2

        A_positive = exp(B)/q * ( exp(b*q) - exp(a*q))
    end

    ### third term.
    term3 = evalDTFTintegralrect(x, a, b)

    return A_positive/4 + A_negative/4 + term3/2
end


function evalDTFTintegralfalling(x::T,
                pass_freq,
                stop_freq,
                integration_lower_limit,
                integration_upper_limit) where T <: Real
    # rename.
    a = integration_lower_limit
    b = integration_upper_limit
    p = pass_freq
    s = stop_freq

    # singularity locations.
    x_0_negative = -1/(2*(s-p))
    x_0_positive = 1/(2*(s-p))

    # common intermediate value.
    k2 = 2*π*x*im

    ### term with negative exponent.
    A_negative::Complex{T} = zero(Complex{T})

    if x_0_negative - eps(T)*2 < x < x_0_negative + eps(T)*2
        # singularity case: when x is x_0_negative.
        A_negative = exp(-p*π*im/(s-p))*(b-a)

    else
        k1 = -π*im
        B = p*k1/(s-p)
        c = -k1/(s-p)
        q = c+k2
        A_negative = exp(B)/q * ( exp(b*q) - exp(a*q))
    end


    ### term with positive exponent.
    A_positive::Complex{T} = zero(Complex{T})

    if x_0_positive - eps(T)*2 < x < x_0_positive + eps(T)*2
        # singularity case: when x is x_0_positive.
        A_positive = exp(p*π*im/(s-p))*(b-a)

    else
        k1 = π*im
        B = p*k1/(s-p)
        c = -k1/(s-p)
        q = c+k2
        A_positive = exp(B)/q * ( exp(b*q) - exp(a*q))
    end

    ### third term.
    term3 = evalDTFTintegralrect(x, a, b)

    return A_positive/4 + A_negative/4 + term3/2
end

function evalDTFTintegralrect(x::T, a, b)::Complex{T} where T <: Real

    q = 2*π*x*im

    if zero(T) - eps(T)*2 < x < zero(T) + eps(T)*2
        # singularity case: when x is 0.
        return exp(q*b)*b - exp(q*a)*a
    end

    return (exp(b*q)-exp(a*q))/q
end

# function evalDTFTintegralfalling(x::T,
#                 pass_freq,
#                 stop_freq,
#                 integration_lower_limit,
#                 integration_upper_limit) where T <: Real
#     # rename.
#     a = integration_lower_limit
#     b = integration_upper_limit
#     p = pass_freq
#     s = stop_freq
#
#     # return evalDTFTintegralrising(x, p, s,
#     #         integration_upper_limit,
#     #         integration_lower_limit)
#
#
#     # singularity locations.
#     x_0_negative = -1/(2*(s-p))
#     x_0_positive = 1/(2*(s-p))
#
#     # common intermediate value.
#     k2 = 2*π*x*im
#
#     ### term with negative exponent.
#     k1 = -π*im
#     B = p*k1/(s-p)
#     c = -k1/(s-p)
#     q = c+k2
#     A_negative = exp(B)/q * ( exp(b*q) - exp(a*q))
#
#     ### term with positive exponent.
#     k1 = π*im
#     B = p*k1/(s-p)
#     c = -k1/(s-p)
#     q = c+k2
#     A_positive = exp(B)/q * ( exp(b*q) - exp(a*q))
#
#     ### third term.
#     q = 2*π*x*im
#     term3 = (exp(b*q)-exp(a*q))/q
#
#     #return A_negative
#     #return A_positive
#     #return A_positive, A_negative, term3
#     return A_positive/4 + A_negative/4 + term3/2
# end
