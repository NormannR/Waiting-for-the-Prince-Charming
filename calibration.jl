import Pkg;
Pkg.activate(".")
# %%
path = pwd() * "/functions"
push!(LOAD_PATH,path)
# %%
using Revise
using Calibration
using JLD2, FileIO
# %%
cal_y = calibrate_y()
# %%
cal_m = calibrate_m()
# %%
cal_c = calibrate_c()
# %%
cal_ssξ = calibrate_ssξ()
# %%
cal_u = calibrate_u()
# %%
JLD2.@save "output/calibrations.jld2" cal_y cal_c cal_m cal_ssξ cal_u
