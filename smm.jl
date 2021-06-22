import Pkg;
Pkg.activate(".")
path = pwd() * "/functions"
push!(LOAD_PATH,path)
# %%
using JLD2, FileIO
using Plots
using LaTeXStrings
using Optim
using StatsBase
pyplot()
# %%
using SharedArrays
using Distributed
addprocs()
# %%
@everywhere begin
import Pkg;
Pkg.activate(".")
path = pwd() * "/functions"
push!(LOAD_PATH,path)
end
# %%
@everywhere begin
using Optim
using SMM
using Macros
using CompStat
using CustomTypes
end
# %%
@load "output/calibrations.jld2" cal_y
cal_y.ϵ = 17/((2013-2005+1)*12)
ss = RegUncSS()
dual = DualSS()
# %%
targets = [
	4/3, 			# F/wp
	0.83,			# jcf / jc
	(1-0.12)*(1-0.26),	# np
	0.12*(1-0.26), 		# nf
	0.7, 			# s / ξ
	1-(1-0.7)^(1/3)		# q
]'
# %%
N = 1000
D = 1000
init_param = deepcopy([ cal_y.s, cal_y.α, cal_y.m, cal_y.b, cal_y.γ, cal_y.F ])
shocks = Vector{Vector{Float64}}(undef, N)
states = Vector{Vector{Int64}}(undef, N)
for k in 1:N
	shocks[k] = Float64[]
	states[k] = Int64[]
	simulate!(shocks[k], states[k], cal_y, D)
end
# %%
n_Δ = 20
grid_Δ = range(0.,stop=0.025,length=n_Δ)
grid_param = SharedArray{Float64,3}((n_Δ, N, 6))
grid_np_sim = SharedArray{Float64,2}((n_Δ, N))
grid_nf_sim = SharedArray{Float64,2}((n_Δ, N))
grid_flag = SharedArray{Bool,2}((n_Δ, N))
grid_np = zeros((n_Δ, 3))
grid_nf = zeros((n_Δ, 3))
init_param = deepcopy([ cal_y.s, cal_y.α, cal_y.m, cal_y.b, cal_y.γ, cal_y.F ])
# %%
@everywhere begin
	cal_y = $cal_y
	targets = $targets
	shocks = $shocks
	states = $states
	ss = $ss
	dual = $dual
	init_param = $init_param
	D = $D
end
# %%
# Non distributed version
# for i in 1:n_Δ
# 	println("$i")
# 	@time begin
# 		for j in 1:N
# 			res = optimize(	x-> loss!(x, cal_y, ss, shocks[j], states[j], grid_Δ[i], D, targets), init_param, NelderMead())
# 			grid_param[i,j,:] .= res.minimizer
# 			grid_flag[i,j] = res.x_converged
# 			@load_vec_into_struct res.minimizer cal_y s α m b γ F  
# 			fill_DualSS!(dual, cal_y, [0.9, 1.3])
# 			grid_np_sim[i,j], grid_nf_sim[i,j] = (dual.np, dual.nf) 
# 		end
# 		grid_np[i,1] = mean(grid_np_sim[i,:])
# 		grid_nf[i,1] = mean(grid_nf_sim[i,:])
# 		grid_np[i,2:3] .= quantile(grid_np_sim[i,:],[0.025,0.975])
# 		grid_nf[i,2:3] .= quantile(grid_nf_sim[i,:],[0.025,0.975])
# 		init_param .= reshape(mean(grid_param[i,:,:], dims = 1),6)
# 	end
# end
# %%
for i in 1:n_Δ
	println("$i")
	@time begin
		@sync @distributed for j in 1:N
			res = optimize(	x-> loss!(x, cal_y, ss, shocks[j], states[j], grid_Δ[i], D, targets), init_param, NelderMead())
			grid_param[i,j,:] .= res.minimizer
			grid_flag[i,j] = Optim.converged(res)
			@load_vec_into_struct res.minimizer cal_y s α m b γ F  
			fill_DualSS!(dual, cal_y, [0.9, 1.3])
			grid_np_sim[i,j], grid_nf_sim[i,j] = (dual.np, dual.nf) 
		end
		grid_np[i,1] = mean(grid_np_sim[i,:])
		grid_nf[i,1] = mean(grid_nf_sim[i,:])
		grid_np[i,2:3] .= quantile(grid_np_sim[i,:],[0.025,0.975])
		grid_nf[i,2:3] .= quantile(grid_nf_sim[i,:],[0.025,0.975])
		init_param .= reshape(mean(grid_param[i,:,:], dims = 1),6)
		@everywhere init_param = $init_param
	end
end
# %%
JLD2.@save "output/smm.jld" grid_Δ grid_flag grid_param grid_np_sim grid_nf_sim grid_np grid_nf D
# %%
@load "output/smm.jld" grid_Δ grid_param grid_np_sim grid_nf_sim grid_np grid_nf D
# %%
p1 = plot(100*grid_Δ,100*(grid_np[:,1] .- targets[3])./targets[3],lab=L"$n^p$", color = :black, line=(:solid, 2))
plot!(p1, 100*grid_Δ,100*(grid_np[:,2:3] .- targets[3])./targets[3],lab="", color = :black, line=(:dash, 2))
ylabel!("Deviation w.r.t to steady-state (%)")
p2 = plot(100*grid_Δ,100*(grid_nf[:,1] .- targets[4])./targets[4],lab=L"$n^f$", color = :black, line=(:solid, 2))
plot!(p2, 100*grid_Δ,100*(grid_nf[:,2:3] .- targets[4])./targets[4],lab="", color = :black, line=(:dash, 2))
plot(p1,p2,layout=(1,2))
xgrid!(:off)
ygrid!(:off)
xlabel!(L"100 \mid \frac{F_i - F}{F} \mid")
# %%
savefig("figures/risk_pol_smm.pdf")
# %%
@load "output/smm.jld" grid_Δ grid_param grid_np_sim grid_nf_sim grid_np grid_nf D
# %%
p1 = plot(100*grid_Δ,100*(grid_np[:,1] .- targets[3])./targets[3],lab="Open-ended", color = :red, line=(:solid, 1))
plot!(p1, 100*grid_Δ,100*(grid_np[:,2:3] .- targets[3])./targets[3],lab="", color = :red, line=(:dash, 1))
ylabel!("Deviation w.r.t the stochastic steady-state (%)")
p2 = plot(100*grid_Δ,100*(grid_nf[:,1] .- targets[4])./targets[4],lab="Fixed-term", color = :blue, line=(:solid, 1))
plot!(p2, 100*grid_Δ,100*(grid_nf[:,2:3] .- targets[4])./targets[4],lab="", color = :blue, line=(:dash, 1))
plot(p1,p2,layout=(1,2))
xgrid!(:off)
ygrid!(:off)
xlabel!(L"100 \mid \frac{F_i - F}{F} \mid")
# %%
savefig("figures/risk_pol_smm_c.pdf")
