module SMM

using NLsolve
using Distributions
using QuadGK
using CompStat
using Macros
using CustomTypes

export graph!, simulate!, fill_RegUncSS!, average, test_eq, loss! 

@def reg_unc begin

    F_1 = (1-Δ)*F
    F_2 = (1+Δ)*F

    U_2 = ((1+ϵ/r)*b + η*γ*(θ_2 + ϵ*0.5*(θ_2+θ_1)/r)/(1-η))/(r+ϵ)
    U_1 = ((1+ϵ/r)*b + η*γ*(θ_1 + ϵ*0.5*(θ_2+θ_1)/r)/(1-η))/(r+ϵ)

    zf_2 = (r+δ+λ)*((r+ϵ)*U_2 - δ*0.5*ϵ*(U_2+U_1)/(r+δ))/(r+δ+ϵ) - λ*Ez/(r+δ)
    zf_1 = (r+δ+λ)*((r+ϵ)*U_1 - δ*0.5*ϵ*(U_2+U_1)/(r+δ))/(r+δ+ϵ) - λ*Ez/(r+δ)

    q_2 = m*θ_2^(-σ)
    q_1 = m*θ_1^(-σ)

    zc_1 = (r+ϵ)*U_1+λ*F_1-0.5*ϵ*(U_1+U_2)-λ*I0(zp_1)/(r+s+λ)
    zc_2 = (r+ϵ)*U_2+λ*F_2-0.5*ϵ*(U_1+U_2)-λ*I0(zp_1)/(r+s+λ)-λ*(I0(zp_2)-I0(zp_1))/(r+s+λ+ϵ/2) + 0.5*ϵ*zc_1/(r+s+λ+ϵ/2)
    zc_2 *= (r+s+λ+ϵ/2)/(r+s+λ+ϵ)

    zs_2 = (zc_2/(r+s+λ) - zf_2/(r+δ+λ))/(1/(r+s+λ) - 1/(r+δ+λ))
    zs_1 = (zc_1/(r+s+λ) - zf_1/(r+δ+λ))/(1/(r+s+λ) - 1/(r+δ+λ))

end

"""
    graph!(t, sim_np, sim_nf, sim_Wp, shocks, states, ss, param, D; k=10)

"""
function graph!(t, sim_np, sim_nf, sim_Wp, shocks, states, ss, param, D; k=10)
    @load_struct param η r s λ α m b γ σ F δ ϵ
    @spec_funs
   
    #Starting from one of the steady-states
    T = shocks[1]
    T_next = shocks[2]
    i0 = states[1]
    @load_index_cat ss i0 old zp zs zc npᶜ np⁰ np nf Wp Wpᶜ Wp⁰ μᵖ μᶠ ξ q

	n = size(shocks,1)
    push!(t,T)
    push!(sim_np, np_old)
    push!(sim_nf, nf_old)
	push!(sim_Wp, Wp_old)

    if n > 2
        for i in 2:n-1
            T = shocks[i]
            T_next = shocks[i+1]
            sh = states[i]

            # println((T, T_next, sh, zp_old, zs_old, zc_old, npᶜ_old, np⁰_old, nf_old, Wpᶜ_old, Wp⁰_old))
            
            @load_index ss sh θ zp zc zs ξ μᵖ μᶠ npᶜ np⁰ np nf u Wpᶜ Wp⁰ Wp F
            F_bar = 0.5*sum(ss.F)
            @dual_to_dual
            @dual_to_dual_wages

            #Updating old values
            zp_old = zp
            zs_old = zs
            zc_old = zc
            npᶜ_old = np_c(T_next-T)
            np⁰_old = np_0(T_next-T)
            nf_old = nᶠ(T_next-T)
            Wpᶜ_old = Wp_c(T_next-T)
            Wp⁰_old = Wp_0(T_next-T)

			time_interval = range(T,T_next,length=k)
			append!(t,time_interval)
			append!(sim_np,nᵖ.(time_interval .- T))
			append!(sim_nf,nᶠ.(time_interval .- T))
			append!(sim_Wp, Wᵖ.(time_interval .- T))

        end
    end
end

"""

    Eq!(R, x, Δ, param)

Fills vector `R` with the residuals of the system characterizing the model with regulatory uncertainty given variables `x`, firing cost wedge `Δ` and parameters `param`.
"""
function Eq!(R, x, Δ, param)

	@load_struct param η r s λ α m b γ σ F δ ϵ

    @spec_funs

    θ_2 = exp(x[1])
    θ_1 = exp(x[2])
    zp_2 = x[3]
    zp_1 = x[4]

    @reg_unc

    R[1] = zp_1-(r+ϵ)*U_1+(r+s+ϵ)*F_1+0.5*ϵ*(U_1+U_2-F_1-F_2)+λ*I0(zp_1)/(r+s+λ)+0.5*ϵ*(zp_1-zp_2)/(r+s+λ+ϵ/2)
    R[2] = zp_2-(r+ϵ)*U_2+(r+s+ϵ)*F_2+0.5*ϵ*(U_1+U_2-F_1-F_2)+λ*I0(zp_1)/(r+s+λ)+λ*(I0(zp_2)-I0(zp_1))/(r+s+λ+ϵ/2)
    R[3] = γ - q_2*(1-η)*( I0(max(zs_2,zc_2))/(r+s+λ) + (I0(zf_2)-I0(max(zs_2,zf_2)))/(r+δ+λ) )
    R[4] = γ - q_1*(1-η)*( I0(max(zs_1,zc_1))/(r+s+λ) + (I0(zf_1)-I0(max(zs_1,zf_1)))/(r+δ+λ) )

end

"""
    fill_RegUncSS!(ss, Δ, param, x0, show_trace)

Fills the `RegUncSS` structure `ss` with the proper steady-state values for both values of steady-state associated with `Δ`.
"""
function fill_RegUncSS!(ss, Δ, param, x0, show_trace)
    R = zeros(4)
    sol = nlsolve((R,x)->Eq!(R,x,Δ,param), x0, show_trace=show_trace)
    θ_2 = exp(sol.zero[1])
    θ_1 = exp(sol.zero[2])
    zp_2 = sol.zero[3]
    zp_1 = sol.zero[4]
    @load_struct param η r s λ α m b γ σ F δ ϵ

    @spec_funs

    @reg_unc

    ss.θ = [θ_1, θ_2]
    ss.zp = [zp_1, zp_2]
    ss.zc = [zc_1, zc_2]
    ss.zf = [zf_1, zf_2]
    ss.zs = [zs_1, zs_2]
    ss.F = [F_1, F_2]

    @load_struct ss θ zp zc zf zs F
    F_bar = 0.5*sum(F)

    p = @. m*θ^(1-σ)
    q = @. p/θ
    ξ = @. s+λ*G(zp)
    μᵖ = @. p*(1-G(max(zs,zc)))
    μᶠ = @. p*(G(max(zs,zf))-G(zf))
    npᶜ = @. λ*(1-G(zp))*μᵖ*δ/((s+λ)*(μᵖ*δ + ξ*δ + μᶠ*ξ))
    np⁰ = @. ξ*μᵖ*δ/((s+λ)*(μᵖ*δ + ξ*δ + μᶠ*ξ))
    np = @. npᶜ + np⁰
    nf = @. μᶠ*ξ/(μᵖ*δ + ξ*δ + μᶠ*ξ)
    u = @. 1-np-nf
    v = @. u*θ
    z_bar_p = @. (λ*np*I(zp)+p*u*I(max(zs,zc)))/((s+λ)*np)
    z_bar_f = @. (λ*nf*Ez+p*u*(I(zf)-I(max(zs,zf))))/((δ+λ)*nf)
    z_bar = @. (np*z_bar_p+nf*z_bar_f)/(np+nf)
    W = @. z_bar_p*np + z_bar_f*nf + b*u - γ*v - λ*G(zp)*np*F
    Eᶜ = @. η*(I(zp)/(1-G(zp)) + (r+s+ϵ)*F - ϵ*F_bar + γ*θ) + (1-η)*b
    E⁰ = @. η*(I(max(zs,zc))/(1-G(max(zs,zc))) - λ*F + γ*θ) + (1-η)*b
    Wpᶜ = @. λ*(1-G(zp))*Eᶜ*np/(s+λ)
    Wp⁰ = @. μᵖ*E⁰*u/(s+λ)
    Wp = @. Wpᶜ + Wp⁰

    ss.flag = ss.flag && converged(sol) && !(any(isnan.(Wp⁰))) && !(any(zs .< zc[2])) && !(any(np⁰ .< 0.01)) && !(any(npᶜ .< 0.01)) && !(any(nf .< 0.01))


    @save_struct ss zp zc zs zf npᶜ np⁰ np nf θ ξ μᵖ μᶠ z_bar_p z_bar_f z_bar W u Wpᶜ Wp⁰ Wp q

    return
end

macro test_regunc_eq(a,regunc,dual)
	expr = Expr(:block)
	for name in intersect(fieldnames(RegUncSS),fieldnames(DualSS))
		push!(expr.args,:(push!($a,$regunc.$name .- $dual.$name)))
	end
	esc(expr)
end

"""

	test_eq(param)

Verifies for parameters `param` that steady-state values of the dual model without uncertainty and of the regulatory uncertainty model with Δ = 0 and ϵ = 0 coincide.
"""
function test_eq(param)
	param.ϵ = 0.	
    Δ = 0.
	regunc = RegUncSS()
	init!(regunc, 2)
	fill_RegUncSS!(regunc, Δ, param, [log(2.2), log(2.2), 0.93, 0.93], false)
	dual = DualSS()
	fill_DualSS!(dual, param, [0.9, 1.3])
	t = []
	@test_regunc_eq t regunc dual
	@assert maximum([(t...)...]) < 1e-7 
    return
end

"""

    simulate!(shocks, states, D)

Simulates a path of shocks on firing costs of duration `D` and stores the resulting shocks and states in `shocks` and `states`.
"""
function simulate!(shocks, states, param, D)
    @load_struct param ϵ
    T = 0.
    c = 1
    u = rand()
    s = (u<0.5)*1+(u>=0.5)*2
    while (T < D)
        push!(shocks,T)
        T += -log(rand())/(ϵ/2)
        c += 1
        push!(states,s)
        s = 3-s
    end
    push!(states,states[end])
    push!(shocks,D)
end

"""

    average(shocks, states, ss, param, D)

Computes the empirical means of targets given a shock instants `shocks`, states `states`, steady-state values `ss`, parameters `param` and duration `D`.
"""
function average(shocks, states, ss, param, D)
    @load_struct param η r s λ α m b γ σ F δ ϵ
    @spec_funs
   
    #Starting from one of the steady-states
    T = shocks[1]
    T_next = shocks[2]
    i0 = states[1]
    @load_index_cat ss i0 old zp zs zc npᶜ np⁰ np nf Wpᶜ Wp⁰ Wp μᵖ μᶠ ξ q

    np_m = np_old*(T_next-T)
    nf_m = nf_old*(T_next-T)
	Fswp_m = F*(T_next-T)/(Wp_old/np_old)
	μᶠsμ_m = μᶠ_old*(T_next-T)/(μᵖ_old+μᶠ_old)
	ssξ_m = s*(T_next-T)/ξ_old
	q_m = q_old*(T_next-T)

    n = size(shocks,1)
    if n > 2
        for i in 2:n-1
            T = shocks[i]
            T_next = shocks[i+1]
            sh = states[i]

            # println((T, T_next, sh, zp_old, zs_old, zc_old, npᶜ_old, np⁰_old, nf_old, Wpᶜ_old, Wp⁰_old))
            
            @load_index ss sh θ zp zc zs ξ μᵖ μᶠ npᶜ np⁰ np nf u Wpᶜ Wp⁰ Wp F
            F_bar = 0.5*sum(ss.F)
            @dual_to_dual
            @dual_to_dual_wages
            
            #Means on [T,T_next]
            np_m += (T_next-T)*np + ρ₁ᵖ*(exp(α₁*(T_next-T))-1)/α₁ + ρ₂ᵖ*(exp(α₂*(T_next-T))-1)/α₂
            nf_m += (T_next-T)*nf + ρ₁ᶠ*(exp(α₁*(T_next-T))-1)/α₁ + ρ₂ᶠ*(exp(α₂*(T_next-T))-1)/α₂
            ssξ_m += (T_next-T)*s/ξ
            μᶠsμ_m += (T_next-T)*μᶠ/(μᶠ+μᵖ)
            q_m += (T_next-T)*m*θ^(-σ)
            
            Δ_wpm1,_ = quadgk(t->F/(Wᵖ(t)/nᵖ(t)), 0, T_next-T)
            Fswp_m += Δ_wpm1

            #Updating old values
            zp_old = zp
            zs_old = zs
            zc_old = zc
            npᶜ_old = np_c(T_next-T)
            np⁰_old = np_0(T_next-T)
            nf_old = nᶠ(T_next-T)
            Wpᶜ_old = Wp_c(T_next-T)
            Wp⁰_old = Wp_0(T_next-T)

        end
    end
    return [Fswp_m, μᶠsμ_m, np_m, nf_m, ssξ_m, q_m] ./ D

end

function moments!(ss, shocks, states, Δ, param, D, show_trace; x0=[log(2.2), log(2.2), 0.93, 0.93])
    mom = Inf*ones(6)
    fill_RegUncSS!(ss, Δ, param, x0, show_trace)
    if ss.flag
        mom .= average(shocks, states, ss, param, D)
    end

	return mom

end

function loss!(x, param, ss, shocks, states, Δ, D, targets; show_trace=false)
    @load_vec_into_struct x param s α m b γ F  
    ss.flag = !any(x .< 0.)
    err = Inf
	if ss.flag
        m = moments!(ss, shocks, states, Δ, param, D, show_trace)
		e = (m' - targets)./targets
		# err = (e*W*e')[1,1]
		err = maximum(abs.(e))
	end
    return err
end

# function EqDual!(R,x,param)

# 	@load_vec param s α m b γ F

# 	Ez = exp(α^2/2)
# 	G(x) = cdf(LogNormal(0,α),x)
# 	g(x) = pdf(LogNormal(0,α),x)
# 	I(x) = exp(α^2/2)*cdf(Normal(0,1),((α^2-log(max(x,0)))/α))
# 	I0(x) = I(x) - x*(1-G(x))

# 	zp = x[1]
# 	θ = x[2]

# 	zc = zp+(r+s+λ)*F
# 	q = m*θ^(-σ)
# 	rU = b + η*γ*θ/(1-η)
# 	zf = rU - λ*(Ez-rU)/(r+δ)
# 	zs = (zp/(r+s+λ) + F - zf/(r+δ+λ))/(1/(r+s+λ) - 1/(r+δ+λ))

# 	R[1] = zp+(r+s)*F+λ*I0(zp)/(r+s+λ)-rU
# 	R[2] = γ - q*(1-η)*( I0(max(zs,zc))/(r+s+λ) + (I0(zf)-I0(max(zs,zf)))/(r+δ+λ) )
# end

# @everywhere function compute_n_dual(sol)

# 	@load_vec param s α m b γ F

# 	Ez = exp(α^2/2)
# 	G(x) = cdf(LogNormal(0,α),x)
# 	g(x) = pdf(LogNormal(0,α),x)
# 	I(x) = exp(α^2/2)*cdf(Normal(0,1),((α^2-log(max(x,0)))/α))
# 	I0(x) = I(x) - x*(1-G(x))

# 	zp = sol[1]
# 	θ = sol[2]

# 	zc = zp+(r+s+λ)*F
# 	q = m*θ^(-σ)
# 	rU = b + η*γ*θ/(1-η)
# 	zf = rU - λ*(Ez-rU)/(r+δ)
# 	zs = (zp/(r+s+λ) + F - zf/(r+δ+λ))/(1/(r+s+λ) - 1/(r+δ+λ))
# 	ξ = s+λ*G(zp)
# 	p = m*θ^(1-σ)
# 	μᵖ = p*(1-G(max(zs,zc)))
# 	μᶠ = p*(G(max(zs,zf))-G(zf))
# 	np = μᵖ*δ/(μᵖ*δ + ξ*δ + μᶠ*ξ)
# 	nf = μᶠ*ξ/(μᵖ*δ + ξ*δ + μᶠ*ξ)

# 	return (np,nf)

# end


end