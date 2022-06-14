
using FFTW, LinearAlgebra

import PyPlot
import BSON

import NMRDataSetup

include("../src/RateConversion.jl")
import .RateConversion

import Random

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


Random.seed!(25)



### bandpass filter specification.

bp_a = 7620
bp_b = 7700

if target_name == "L-Arginine"
    bp_a = ppm2hzfunc(1.85)
    bp_b = ppm2hzfunc(1.97)
end

if target_name == "L-Serine"
    bp_a = ppm2hzfunc(3.9)
    bp_b = ppm2hzfunc(4.0)
end

bp_δ = (bp_b-bp_a)/10
@assert bp_δ > 0

N_h = 2001

# N_cost = 200
# U_cost = LinRange(bp_a, bp_b, N_cost)
# S_U_cost = DTFT_s.(U_cost)

### get bandpass filter.

h_func = xx->RateConversion.cosinetransitionbandpassimpulsefunc(xx, bp_a, bp_b, bp_δ, fs)
t_h = RateConversion.gettimerangetunablefilter(N_h, fs)
h = h_func.(t_h)
hs = RateConversion.fastlinearconv(s_t, h)

T = Float64
M = round(Int, (length(h)-1)/2) # for odd length(r)
N = length(s_t)
t_hs_idx = -M:1:(M+N-1)
Ts = one(T)/fs
t_hs = Ts .* t_hs_idx

DTFT_h = vv->RateConversion.computeDTFTch3eq29(h, vv, t_h)
DTFT_hs = vv->RateConversion.computeDTFTch3eq29(hs, vv, t_hs)

#t2 = t_hs[idx_range]

### end get and apply bandpass.

N_viz = 5000
#U = LinRange(bp_a-500, bp_b+500, N_viz)
#P = hz2ppmfunc.(U)

#P = LinRange(0.5, 2.0, N_viz)
#U = ppm2hzfunc.(P)

U = LinRange(0, fs, N_viz)
P = hz2ppmfunc.(U)

S_U = DTFT_s.(U)
HS_U = DTFT_hs.(U)


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, abs.(S_U), label = "S")
PyPlot.plot(P, abs.(HS_U), label = "filtered S")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("abs")
PyPlot.title("data vs. bandpassed")



PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(S_U), label = "S")
PyPlot.plot(P, real.(HS_U), label = "filtered S")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data vs. bandpassed")
