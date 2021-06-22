# %%
import Pkg;
Pkg.activate(".")
# %%
path = pwd() * "/functions"
push!(LOAD_PATH,path)
# %%
using Revise
using JLD2, FileIO
using CompStat
using CustomTypes
# %%
JLD2.@load "output/calibrations.jld2"
# %%
using Plots
pyplot()
# %%
plot_cmpstat_bw(cal_y, "y")
# %%
plot_cmpstat_bw(cal_m, "m", x1 = [0.9, 1.3])
# %%
plot_cmpstat_bw(cal_c, "c")
# %%
plot_cmpstat_bw(cal_ssξ, "ssxi")
# %%
plot_cmpstat_bw(cal_u, "u", x0=[0.9, 6.95], x1=[1.33, 6.95])
# %%
plot_cmpstat_colors(cal_y, "y")
plot_cmpstat_colors(cal_m, "m", x1 = [0.9, 1.3])
plot_cmpstat_colors(cal_c, "c")
plot_cmpstat_colors(cal_ssξ, "ssxi")
plot_cmpstat_colors(cal_u, "u", x0=[0.9, 6.95], x1=[1.33, 6.95])
# %%
plot_transitions(cal_y.F*0.5, cal_y, "y")
# %%
