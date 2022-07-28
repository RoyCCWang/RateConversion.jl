
using FFTW, LinearAlgebra

import PyPlot
import BSON

import NMRDataSetup

include("../src/RateConversion.jl")
import .RateConversion

import Random

include("helpers/utils.jl")

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

t_hs, mhs0, hs, h = filtermodulatesequence(s_t, fs, bp_a, bp_b)

ind = findfirst(xx->xx>0, t_hs)

# #t_range = (ind+800):(ind+1200)
# t_range = ind:length(t_hs)
# t_mhs = t_hs[t_range]
# mhs = mhs0[t_range]
#
# PyPlot.figure(fig_num)
# fig_num += 1
#
# PyPlot.plot(t_mhs, real.(mhs), label = "mhs")
# #PyPlot.plot(t, real.(s_t), label = "s")
#
# PyPlot.legend()
# PyPlot.xlabel("time (sec)")
# PyPlot.ylabel("")
# PyPlot.title("real part of mhs")




PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(t, real.(s_t), label = "s")
PyPlot.plot(t_hs, real.(hs), label = "filtered s")

PyPlot.legend()
PyPlot.xlabel("sec")
PyPlot.ylabel("real")
PyPlot.title("s vs. bandpassed")



### rough envelop func.

# only when there is zero phase for all resonance comonents.
A0 = sum( sum(sum(As[n].αs[i]) for i = 1:length(As[n].αs)) for n = 1:length(As))

λ = λ0
A_func = tt->A0*exp(-λ*tt)

#t_hs[1001] #  this is 0.0, if N_h == 2001.

A_t = A_func.(t_hs)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(t_hs, real.(mhs0), label = "real mhs0")
PyPlot.plot(t_hs, imag.(mhs0), label = "imag mhs0")
PyPlot.plot(t_hs, abs.(mhs0), label = "abs mhs0")
PyPlot.plot(t_hs, A_t, label = "A(t)")

PyPlot.legend()
PyPlot.xlabel("time (sec)")
PyPlot.ylabel("")
PyPlot.title("envelope")

@assert 1==2

A_t = A_func.(t_hs)

PyPlot.figure(fig_num)
fig_num += 1

#PyPlot.plot(t_mhs, abs.(mhs), label = "abs mhs")
PyPlot.plot(t_hs, real.(hs), label = "real hs")
PyPlot.plot(t_hs, imag.(hs), label = "imag hs")
PyPlot.plot(t_hs, abs.(hs), label = "abs hs")
PyPlot.plot(t_hs, A_t, label = "A(t)")

PyPlot.legend()
PyPlot.xlabel("time (sec)")
PyPlot.ylabel("")
PyPlot.title("envelope")


@assert 1==2

A_t = A_func.(t)

PyPlot.figure(fig_num)
fig_num += 1

#PyPlot.plot(t_mhs, abs.(mhs), label = "abs mhs")
PyPlot.plot(t, real.(s_t), label = "real s")
PyPlot.plot(t, imag.(s_t), label = "imag s")
PyPlot.plot(t, abs.(s_t), label = "abs s")
PyPlot.plot(t, A_t, label = "A(t)")

PyPlot.legend()
PyPlot.xlabel("time (sec)")
PyPlot.ylabel("")
PyPlot.title("envelope")

# TODO is it possible to get A0 when there are phase?
# what about known phase?
