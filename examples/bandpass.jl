
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


s_t = q_t
t = t_test
DTFT_s = vv->RateConversion.computeDTFTch3eq29(s_t, vv, t)

### bandpass filter specification.

# # beta-glucose.
# bp_a = ppm2hzfunc(3.2)
# bp_b = ppm2hzfunc(3.26)

# serine.
# bp_a = ppm2hzfunc(3.915)
# bp_b = ppm2hzfunc(4.0)
bp_a = ppm2hzfunc(3.5)
bp_b = ppm2hzfunc(4.5)

# bp_a = 7620
# bp_b = 7700
#
# if target_name == "L-Arginine"
#     bp_a = ppm2hzfunc(1.85)
#     bp_b = ppm2hzfunc(1.97)
# end
#
# if target_name == "L-Serine"
#     bp_a = ppm2hzfunc(3.9)
#     bp_b = ppm2hzfunc(4.0)
# end

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

PyPlot.plot(t, real.(s_t), label = "s")
PyPlot.plot(t_hs, real.(hs), label = "filtered s")

PyPlot.legend()
PyPlot.xlabel("sec")
PyPlot.ylabel("real")
PyPlot.title("data vs. bandpassed")


### modulate.
# modulate.
modulation_shift = -(bp_a-bp_δ) #bp_δ*2 # Hz.
mhs0 = RateConversion.modulatediscretesignalcomplex(hs, t_hs, modulation_shift)

t_mhs = t_hs
DTFT_mhs0 = vv->RateConversion.computeDTFTch3eq29(mhs0, vv, t_mhs)

MHS0_U = DTFT_mhs0.(U)


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(U, abs.(S_U), label = "S")
PyPlot.plot(U, abs.(MHS0_U), label = "modulated filtered S")

PyPlot.legend()
PyPlot.xlabel("Hz")
PyPlot.ylabel("abs")
PyPlot.title("data vs. mb")

#@assert 1==2

mhs = mhs0 ./ maximum(abs.(mhs0))

# PyPlot.figure(fig_num)
# fig_num += 1
#
# PyPlot.plot(t_mhs, real.(mhs), label = "real")
# PyPlot.plot(t_mhs, imag.(mhs), label = "imag")
# PyPlot.plot(t_mhs, abs.(mhs), label = "abs")
#
# PyPlot.legend()
# PyPlot.xlabel("time (sec)")
# PyPlot.ylabel("")
# PyPlot.title("mhs")

ind = findfirst(xx->xx>0, t_hs)

t_range = (ind+800):(ind+1200)
#t_range = ind:length(t_hs)
t_mhs = t_hs[t_range]
mhs = mhs0[t_range]

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(t_mhs, real.(mhs), label = "mhs")
#PyPlot.plot(t, real.(s_t), label = "s")
PyPlot.plot(t_hs, real.(hs), label = "filtered s")

PyPlot.legend()
PyPlot.xlabel("time (sec)")
PyPlot.ylabel("")
PyPlot.title("real part of mhs")
