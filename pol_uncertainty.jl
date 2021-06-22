# %%
import Pkg;
Pkg.activate(".")
# %%
path = pwd() * "/functions"
push!(LOAD_PATH,path)
# %%
using Revise
# %%
using CustomTypes
using CompStat
using SMM
using JLD2, FileIO
# %%
using Plots
using LaTeXStrings
pyplot()
# %%
@load "output/calibrations.jld2" cal_y
# %%
cal_y.ϵ = 17/((2013-2005+1)*12)
# %%
cal = deepcopy(cal_y)
# %%
D = 1000
Δ = 0.025
# %%
shocks = Float64[]
states = Int64[]
simulate!(shocks, states, cal_y, D)
# %%
t = Float64[]
sim_np = Float64[]
sim_nf = Float64[]
sim_Wp = Float64[]
# %%
ss = RegUncSS()
ss.flag = true
fill_RegUncSS!(ss, Δ, cal_y, [log(2.2), log(2.2), 0.93, 0.93], false)
# %%
graph!(t, sim_np, sim_nf, sim_Wp, shocks, states, ss, cal_y, D, k=100)
# %%
p1 = plot(shocks,states,lab="",color=:black,line=(:solid, 1),linetype=:steppost)
yticks!([1.,2.],[L"F_1",L"F_2"])
ylabel!(L"F")
p2 = plot(t,sim_np,lab="",color=:black,line=(:solid, 1))
ylabel!(L"n^p")
p3 = plot(t,sim_nf,lab="",color=:black,line=(:solid, 1))
ylabel!(L"n^f")
p4 = plot(t,sim_Wp,lab="",color=:black,line=(:solid, 1))
ylabel!(L"W^p")
xlabel!("Time (months)")
# %%
p = plot(p1,p2,p3, layout=(3,1))
xgrid!(:off)
ygrid!(:off)
# %%
savefig("figures/risk_pol_sim.pdf")
