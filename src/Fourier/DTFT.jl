


### legacy.
function computeDTFTviaformula(h::AbstractArray{Complex{T}}, ω::T)::Complex{T} where T <: Real
    N = length(h)
    return sum( h[n+1]*exp(-im*ω*n) for n = 0:N-1 )
end

function computeDTFTviaformula(h::AbstractArray{T}, ω::T)::Complex{T} where T <: Real
    N = length(h)
    return sum( h[n+1]*exp(-im*ω*n) for n = 0:N-1 )
end


### current.

# Eq'n 6.7 from Dubois' book. Assumes continuous signal is bandlimited.
function evalDTFTch6eq7( Fc::Function,
                        u::T,
                        Λ::LinRange{T2})::Complex{T} where {T <: Real, T2 <: Real}
    Ts = Λ[2]-Λ[1]
    fs = 1/Ts
    t_elapsed = Λ[end]-Λ[1]

    running_sum = zero(T)
    #for i = 0:length(Λ)-1
    for i = -length(Λ):length(Λ)
    #for i = -1000000:1000000
        r = fs*i

        running_sum += Fc(u-r)
    end

    return running_sum*fs
end

# case 1D.
function computeDTFTch3eq29(h, u::T, Λ)::Complex{T} where T <: Real

    # debug_array = zeros(Complex{T}, length(Λ))
    # debug_array2 = zeros(Complex{T}, length(Λ))
    # fs = 1/(Λ[2]-Λ[1])

    running_sum = zero(T)
    for i = 1:length(Λ)
        x = Λ[i]

        running_sum += h[i]*exp(-im*2*π*u*x)


        # # debug.
        # debug_array[i] = h[i]*exp(-im*2*π*u*x)
        #
        # #debug_array2[i] = h[i]*exp(-im*2*π*u*x +im*35*2*pi)
        # debug_array2[i] = h[i]*exp(-im*2*π*(u+3*fs)*x)
    end

    #println("debug: norm(debug_array -debug_array2 ) = ", norm(debug_array -debug_array2 ) )
    return running_sum
end

# case 1D.
function computeinvDTFTch3eq31(F::Function, x::T, Λ;
                integral_maxevals = 10000,
                integral_initial_divisions = 1)::Tuple{Complex{T},T} where T <: Real

    # debug_array = zeros(Complex{T}, length(Λ))
    # debug_array2 = zeros(Complex{T}, length(Λ))
    fs = 1/(Λ[2]-Λ[1])

    h = uu->F(uu)*exp(im*2*π*uu*x)
    NI_val, NI_err = HCubature.hquadrature(h,
                                        -fs/2,
                                        fs/2;
                                        norm = norm,
                                        atol = 0,
                                        maxevals = integral_maxevals,
                                        initdiv = integral_initial_divisions)
    #
    out = NI_val/fs

    return out, NI_err
end

# case 1D.
function linearconvch3eq26( f::Vector{T},
                            h,
                            out_range::StepRange{Int,Int})::Vector{T} where T
    #
    N = length(out_range)
    out = Vector{T}(undef, N)

    #
    M = fld(length(h)-1,2)
    h_range = -M:1:M

    #
    for i = 1:N

        x = out_range[i]

        running_sum = zero(T)
        for y in h_range

            factor1 = evalsymsequence(h, y)
            factor2 = evalcasualsequence(f, x-y)

            running_sum += factor1*factor2

        end
        out[i] = running_sum

    end

    return out
end

# evaluate a compactly supported function f on the unit integer lattice.
# The non-zero entries of f are stored in a.
function evalsequence(a::Vector{T}, zero_idx::Int, x::Int)::T where T
    idx = x+zero_idx
    if 1 <= idx <= length(a)
        return a[idx]
    end

    return zero(T)
end

# assumes first entry in the array a corresponds to f[x] = f[0].
function evalcasualsequence(a::Vector{T}, x::Int)::T where T
    return evalsequence(a, 1, x)
end

# assumes first entry in the array a corresponds to f[x] = f[-M/2],
#   M := length(a) -1, M is even.
function evalsymsequence(a::Vector{T}, x::Int)::T where T
    @assert isodd(length(a))

    return evalsequence(a, fld(length(a),2)+1, x)
end
