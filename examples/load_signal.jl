
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

### user inputs.


#save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments" # where to store the loaded experiment objects.

target_name = "L-Serine"

experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/NRC_4_amino_acid_mixture_Jan_2022/1" # where the Bruker experiment data folder is located.
project_name = "NRC-Jan2022-serine" # project folder name within `save_dir` for this experiment.

solvent_ppm_guess = 4.7 # where the solvent's peak frequency might be. In units of ppm.
solvent_window_ppm = 0.1 # This library will try to search for a resonance component +/- `solvent_window_ppm` about `solvent_ppm_guess`.

### end inputs.

## load.
println("load timing:")
@time s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW, α_0ppm, β_0ppm, λ_0ppm, Ω_0ppm,
    results_0ppm, dic, α_solvent, β_solvent, λ_solvent, Ω_solvent,
    results_solvent = NMRDataSetup.loadspectrum(experiment_full_path;
    solvent_ppm = solvent_ppm_guess,
    solvent_window_ppm = solvent_window_ppm)




##### fourier.

### DTFT
t = RateConversion.gettimerange(length(s_t), fs)
DTFT_s = vv->RateConversion.computeDTFTch3eq29(s_t, vv, t)

# DTFT vs. DFT. sanity check.
U_DFT = LinRange(0, fs-fs/length(s_t), length(s_t))
DFT_s = fft(s_t)
S_U = DTFT_s.(U_DFT)

discrepancy = DFT_s - S_U
relative_discrepancy = norm(discrepancy)/norm(DFT_s)
println("DFT, DTFT discrepancy: ", relative_discrepancy)
println()

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(U_DFT, abs.(S_U), label = "DTFT S", linewidth = "2")
PyPlot.plot(U_DFT, abs.(DFT_s), "--", label = "DFT S", linewidth = "2")

PyPlot.legend()
PyPlot.xlabel("Hz")
PyPlot.ylabel("abs")
PyPlot.title("DFT vs. DTFT")

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(U_DFT, real.(S_U), label = "DTFT S", linewidth = "2")
PyPlot.plot(U_DFT, real.(DFT_s), "--", label = "DFT S", linewidth = "2")

PyPlot.legend()
PyPlot.xlabel("Hz")
PyPlot.ylabel("real")
PyPlot.title("DFT vs. DTFT")



PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(U_DFT, imag.(S_U), label = "DTFT S", linewidth = "2")
PyPlot.plot(U_DFT, imag.(DFT_s), "--", label = "DFT S", linewidth = "2")

PyPlot.legend()
PyPlot.xlabel("Hz")
PyPlot.ylabel("imag")
PyPlot.title("DFT vs. DTFT")
