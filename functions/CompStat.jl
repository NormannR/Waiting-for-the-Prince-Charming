module CompStat

using Plots
pyplot()

using LaTeXStrings
using Roots
using NLsolve
using Distributions
using Macros
using CustomTypes

export plot_cmpstat_bw, plot_cmpstat_colors, plot_transitions, @spec_funs, @dual_to_dual, @dual_to_dual_wages, fill_DualSS!

function max_zeros(z)
    if isempty(z)
        return 0.
    else
        return maximum(z)
    end
end

@def spec_funs begin
    Ez = exp(α^2/2)
    G(x) = cdf(LogNormal(0,α),x)
    g(x) = pdf(LogNormal(0,α),x)
    I(x) = exp(α^2/2)*cdf(Normal(0,1),((α^2-log(max(x,0)))/α))
    I0(x) = I(x) - x*(1-G(x))
end

@def classic_to_classic begin
    α₀ = (G(max(zp,zp_old))-G(zp_old))/(1-G(zp_old))
    β₀ = (G(max(zc, zc_old))-G(zc_old))/(1-G(zc_old))
    c₀ᶜ = (1-α₀)*nᶜ_old
    c₀⁰ = (1-β₀)*n⁰_old
    ω = -(s+λ)

    ρ₀ᶜ = (c₀ᶜ*μ + c₀⁰*(-ξ - ω) + nᶜ*(-μ - ξ))/(μ + ξ + ω)
    ρ₁ᶜ = c₀ᶜ - nᶜ - ρ₀ᶜ
    
    ρ₀⁰ = -(c₀ᶜ*μ + c₀⁰*(-ξ - ω) + n⁰*μ + n⁰*ξ - μ)/(μ + ξ + ω)
    ρ₁⁰ = c₀⁰ - n⁰ - ρ₀⁰

    np_c(t)  = nᶜ + ρ₀ᶜ*exp(ω*t) + ρ₁ᶜ*exp(-(μ+ξ)*t) 
    np_0(t)  = n⁰ + ρ₀⁰*exp(ω*t) + ρ₁⁰*exp(-(μ+ξ)*t) 
    nᵖ(t) = np_c(t) + np_0(t)
end

@def classic_to_dual begin

    α₁ = -δ/2 - μᵖ/2 - μᶠ/2 - ξ/2 + sqrt(δ^2 - 2*δ*μᵖ + 2*δ*μᶠ - 2*δ*ξ + μᵖ^2 + 2*μᵖ*μᶠ + 2*μᵖ*ξ + μᶠ^2 - 2*μᶠ*ξ + ξ^2)/2
    α₂ = -δ/2 - μᵖ/2 - μᶠ/2 - ξ/2 - sqrt(δ^2 - 2*δ*μᵖ + 2*δ*μᶠ - 2*δ*ξ + μᵖ^2 + 2*μᵖ*μᶠ + 2*μᵖ*ξ + μᶠ^2 - 2*μᶠ*ξ + ξ^2)/2
    α₀ = (G(max(zp, zp_old))-G(zp_old))/(1-G(zp_old))
    β₀ = (G(max(zc, zc_old))-G(zc_old))/(1-G(zc_old))
    c₀ᶜ = (1-α₀)*nᶜ_old
    c₀⁰ = (1-β₀)*n⁰_old
    ω = -(s+λ)

    ρ₀ᶜ = (c₀ᶜ*(δ*μᵖ + μᵖ*ω) + c₀⁰*(-δ*ξ - δ*ω - μᶠ*ξ - μᶠ*ω - ξ*ω - ω^2) - npᶜ*α₁*α₂ + nf_old*(μᵖ*ξ + μᵖ*ω) - μᵖ*ξ - μᵖ*ω)/(α₁*α₂ - α₁*ω - α₂*ω + ω^2)
    ρ₁ᶜ = (c₀ᶜ*(α₁ + δ + μᵖ + μᶠ) + c₀⁰*(-ξ - ω) + npᶜ*α₂ + ρ₀ᶜ*(α₂ - ω))/(α₁ - α₂) 
    ρ₂ᶜ = c₀ᶜ - npᶜ - ρ₀ᶜ - ρ₁ᶜ 
    
    ρ₀⁰ = -(c₀ᶜ*(δ*μᵖ + μᵖ*ω) + c₀⁰*(-δ*ξ - δ*ω - μᶠ*ξ - μᶠ*ω - ξ*ω - ω^2) + np⁰*α₁*α₂ + nf_old*(μᵖ*ξ + μᵖ*ω) - δ*μᵖ - μᵖ*ξ - μᵖ*ω)/(α₁*α₂ - α₁*ω - α₂*ω + ω^2)
    ρ₁⁰ = -(c₀ᶜ*μᵖ + c₀⁰*(-α₁ - δ - μᶠ - ξ - ω) - np⁰*α₂ + nf_old*μᵖ - μᵖ + ρ₀⁰*(-α₂ + ω))/(α₁ - α₂)
    ρ₂⁰ = c₀⁰ - np⁰ - ρ₀⁰ - ρ₁⁰

    ρ₁ᶠ = -(c₀ᶜ*μᶠ + c₀⁰*μᶠ - nf*α₂ + nf_old*(-α₁ - μᵖ - ξ) - μᶠ)/(α₁ - α₂)
    ρ₂ᶠ = -nf + nf_old - ρ₁ᶠ

    ρ₁ᵖ = ρ₁ᶜ + ρ₁⁰
    ρ₂ᵖ = ρ₂ᶜ + ρ₂⁰

    np_c(t)  = npᶜ + ρ₀ᶜ*exp(ω*t) + ρ₁ᶜ*exp(α₁*t) +  ρ₂ᶜ*exp(α₂*t) 
    np_0(t)  = np⁰ + ρ₀⁰*exp(ω*t) + ρ₁⁰*exp(α₁*t) +  ρ₂⁰*exp(α₂*t) 
    nᵖ(t) = np_c(t) + np_0(t)
    nᶠ(t) = nf + ρ₁ᶠ*exp(α₁*t) + ρ₂ᶠ*exp(α₂*t)
    u(t) = 1 - nᵖ(t) - nᶠ(t)
end

@def dual_to_dual begin

    α₁ = -δ/2 - μᵖ/2 - μᶠ/2 - ξ/2 + sqrt(δ^2 - 2*δ*μᵖ + 2*δ*μᶠ - 2*δ*ξ + μᵖ^2 + 2*μᵖ*μᶠ + 2*μᵖ*ξ + μᶠ^2 - 2*μᶠ*ξ + ξ^2)/2
    α₂ = -δ/2 - μᵖ/2 - μᶠ/2 - ξ/2 - sqrt(δ^2 - 2*δ*μᵖ + 2*δ*μᶠ - 2*δ*ξ + μᵖ^2 + 2*μᵖ*μᶠ + 2*μᵖ*ξ + μᶠ^2 - 2*μᶠ*ξ + ξ^2)/2
    α₀ = (G(max(zp,zp_old))-G(zp_old))/(1-G(zp_old))
    β₀ = (G(max(zc, zc_old, zs_old))-G(max(zc_old, zs_old)))/(1-G(max(zc_old, zs_old)))
    c₀ᶜ = (1-α₀)*npᶜ_old
    c₀⁰ = (1-β₀)*np⁰_old
    ω = -(s+λ)

    ρ₀ᶜ = (c₀ᶜ*(δ*μᵖ + μᵖ*ω) + c₀⁰*(-δ*ξ - δ*ω - μᶠ*ξ - μᶠ*ω - ξ*ω - ω^2) - npᶜ*α₁*α₂ + nf_old*(μᵖ*ξ + μᵖ*ω) - μᵖ*ξ - μᵖ*ω)/(α₁*α₂ - α₁*ω - α₂*ω + ω^2)
    ρ₁ᶜ = (c₀ᶜ*(α₁ + δ + μᵖ + μᶠ) + c₀⁰*(-ξ - ω) + npᶜ*α₂ + ρ₀ᶜ*(α₂ - ω))/(α₁ - α₂) 
    ρ₂ᶜ = c₀ᶜ - npᶜ - ρ₀ᶜ - ρ₁ᶜ 
    
    ρ₀⁰ = -(c₀ᶜ*(δ*μᵖ + μᵖ*ω) + c₀⁰*(-δ*ξ - δ*ω - μᶠ*ξ - μᶠ*ω - ξ*ω - ω^2) + np⁰*α₁*α₂ + nf_old*(μᵖ*ξ + μᵖ*ω) - δ*μᵖ - μᵖ*ξ - μᵖ*ω)/(α₁*α₂ - α₁*ω - α₂*ω + ω^2)
    ρ₁⁰ = -(c₀ᶜ*μᵖ + c₀⁰*(-α₁ - δ - μᶠ - ξ - ω) - np⁰*α₂ + nf_old*μᵖ - μᵖ + ρ₀⁰*(-α₂ + ω))/(α₁ - α₂)
    ρ₂⁰ = c₀⁰ - np⁰ - ρ₀⁰ - ρ₁⁰

    ρ₁ᶠ = -(c₀ᶜ*μᶠ + c₀⁰*μᶠ - nf*α₂ + nf_old*(-α₁ - μᵖ - ξ) - μᶠ)/(α₁ - α₂)
    ρ₂ᶠ = -nf + nf_old - ρ₁ᶠ

    ρ₁ᵖ = ρ₁ᶜ + ρ₁⁰
    ρ₂ᵖ = ρ₂ᶜ + ρ₂⁰

    np_c(t)  = npᶜ + ρ₀ᶜ*exp(ω*t) + ρ₁ᶜ*exp(α₁*t) +  ρ₂ᶜ*exp(α₂*t) 
    np_0(t)  = np⁰ + ρ₀⁰*exp(ω*t) + ρ₁⁰*exp(α₁*t) +  ρ₂⁰*exp(α₂*t) 
    nᵖ(t) = np_c(t) + np_0(t)
    nᶠ(t) = nf + ρ₁ᶠ*exp(α₁*t) + ρ₂ᶠ*exp(α₂*t)
    u(t) = 1 - nᵖ(t) - nᶠ(t)
end

@def dual_to_dual_wages begin

    Eᶜ = η*(I(zp)/(1-G(zp)) + (r+s+ϵ)*F - ϵ*F_bar + γ*θ) + (1-η)*b
    E₀ᶜ = 0.
    if zp > zp_old
        E₀ᶜ = η*((I(zp_old) - I(zp))/(G(zp) - G(zp_old)) + (r+s+ϵ)*F - ϵ*F_bar + γ*θ) + (1-η)*b
    end
    E⁰ = η*(I(zs)/(1-G(zs)) - λ*F + γ*θ) + (1-η)*b
    E₀⁰ = 0.
    if zc > zs_old
        E₀⁰ = η*((I(zs_old) - I(zc))/(G(zc)-G(zs_old)) + γ*θ) + (1-η)*b
    end

    Wpᶜ = (Eᶜ*npᶜ*ξ + Eᶜ*npᶜ*ω + Eᶜ*np⁰*ξ + Eᶜ*np⁰*ω)/ω
    ν₀ᶜ = -Eᶜ*ξ*ρ₀ᶜ - Eᶜ*ξ*ρ₀⁰ - Eᶜ*ρ₀ᶜ*ω - Eᶜ*ρ₀⁰*ω
    ν₁ᶜ = -Eᶜ*npᶜ - Eᶜ*np⁰ - Eᶜ*ξ*(npᶜ + np⁰)/ω - Eᶜ*ρ₁ᶜ - Eᶜ*ρ₁⁰ - Eᶜ*ρ₂ᶜ - Eᶜ*ρ₂⁰ + Eᶜ*(α₁ + ξ)*(ρ₁ᶜ + ρ₁⁰)/(α₁ - ω) + Eᶜ*(α₂ + ξ)*(ρ₂ᶜ + ρ₂⁰)/(α₂ - ω) - E₀ᶜ*npᶜ_old*α₀ + Wpᶜ_old
    ν₂ᶜ = (-Eᶜ*ξ*ρ₁ᶜ - Eᶜ*ξ*ρ₁⁰ - Eᶜ*ρ₁ᶜ*ω - Eᶜ*ρ₁⁰*ω)/(α₁ - ω)
    ν₃ᶜ = (-Eᶜ*ξ*ρ₂ᶜ - Eᶜ*ξ*ρ₂⁰ - Eᶜ*ρ₂ᶜ*ω - Eᶜ*ρ₂⁰*ω)/(α₂ - ω)

    Wp⁰ = (E⁰*npᶜ*μᵖ + E⁰*nf*μᵖ + E⁰*np⁰*μᵖ - E⁰*μᵖ)/ω
    ν₀⁰ = -E⁰*μᵖ*ρ₀ᶜ - E⁰*μᵖ*ρ₀⁰
    ν₁⁰ = E⁰*μᵖ*(ρ₂ᶜ + ρ₂ᶠ + ρ₂⁰)/(α₂ - ω) + E⁰*μᵖ*(ρ₁ᶜ + ρ₁ᶠ + ρ₁⁰)/(α₁ - ω) - E⁰*μᵖ*(npᶜ + nf + np⁰ - 1)/ω - E₀⁰*np⁰_old*β₀ + Wp⁰_old
    ν₂⁰ = (-E⁰*μᵖ*ρ₁ᶜ - E⁰*μᵖ*ρ₁ᶠ - E⁰*μᵖ*ρ₁⁰)/(α₁ - ω)
    ν₃⁰ = (-E⁰*μᵖ*ρ₂ᶜ - E⁰*μᵖ*ρ₂ᶠ - E⁰*μᵖ*ρ₂⁰)/(α₂ - ω) 

    Wp_c(t) = Wpᶜ + ν₀ᶜ*t*exp(ω*t) + ν₁ᶜ*exp(ω*t) + ν₂ᶜ*exp(α₁*t) + ν₃ᶜ*exp(α₂*t)
    Wp_0(t) = Wp⁰ + ν₀⁰*t*exp(ω*t) + ν₁⁰*exp(ω*t) + ν₂⁰*exp(α₁*t) + ν₃⁰*exp(α₂*t)
    Wᵖ(t) = Wp_c(t) + Wp_0(t)

end

@def dual_to_classic begin
    α₀ = (G(max(zp,zp_old))-G(zp_old))/(1-G(zp_old))
    β₀ = (G(max(zc, zc_old))-G(zc_old))/(1-G(zc_old))
    c₀ᶜ = (1-α₀)*npᶜ_old
    c₀⁰ = (1-β₀)*np⁰_old
    ω = -(s+λ)

    ρ₀ᶜ = (c₀ᶜ*(δ*μ + μ*ω) + c₀⁰*(-δ*ξ - δ*ω - ξ*ω - ω^2) + nᶜ*(-δ*μ - δ*ξ) + nf_old*(μ*ξ + μ*ω) - μ*ξ - μ*ω)/(δ*μ + δ*ξ + δ*ω + μ*ω + ξ*ω + ω^2)
    ρ₁ᶜ = (c₀ᶜ*(δ - ξ) + c₀⁰*(-ξ - ω) - nᶜ*δ + ρ₀ᶜ*(-δ - ω))/(δ - μ - ξ)
    ρ₂ᶜ = c₀ᶜ - nᶜ - ρ₀ᶜ - ρ₁ᶜ

    ρ₀⁰ = -(c₀ᶜ*(δ*μ + μ*ω) + c₀⁰*(-δ*ξ - δ*ω - ξ*ω - ω^2) + n⁰*δ*μ + n⁰*δ*ξ + nf_old*(μ*ξ + μ*ω) - δ*μ - μ*ξ - μ*ω)/(δ*μ + δ*ξ + δ*ω + μ*ω + ξ*ω + ω^2)
    ρ₁⁰ = -(c₀ᶜ*μ + c₀⁰*(-δ + μ - ω) + n⁰*δ + nf_old*μ - μ + ρ₀⁰*(δ + ω))/(δ - μ - ξ)
    ρ₂⁰ = c₀⁰ - n⁰ - ρ₀⁰ - ρ₁⁰

    ρ₁ᵖ = ρ₁ᶜ + ρ₁⁰
    ρ₂ᵖ = ρ₂ᶜ + ρ₂⁰

    np_c(t)  = nᶜ + ρ₀ᶜ*exp(ω*t) + ρ₁ᶜ*exp(-(μ + ξ)*t) +  ρ₂ᶜ*exp(-δ*t) 
    np_0(t)  = n⁰ + ρ₀⁰*exp(ω*t) + ρ₁⁰*exp(-(μ + ξ)*t) +  ρ₂⁰*exp(-δ*t) 
    nᵖ(t) = np_c(t) + np_0(t)
    nᶠ(t) = nf_old*exp(-δ*t)
    u(t) = 1 - nᵖ(t) - nᶠ(t)

end

"""
    compute_dual!(ss, sol, param)

Fills the `DualSS` structure with the dual-equilibrium thresholds and employments given basic equilibrium values `sol` and parameters `param`
"""
function compute_dual!(ss, sol, param; index=nothing)

    Macros.@load_struct param η r s λ α m b γ σ F δ
    @spec_funs

    zp = sol[1]
    θ = sol[2]

    zc = zp+(r+s+λ)*F
    q = m*θ^(-σ)
    rU = b + η*γ*θ/(1-η)
    zf = rU - λ*(Ez-rU)/(r+δ)
    zs = (zp/(r+s+λ) + F - zf/(r+δ+λ))/(1/(r+s+λ) - 1/(r+δ+λ))
    ξ = s+λ*G(zp)
    p = m*θ^(1-σ)
    μᵖ = p*(1-G(max(zs,zc)))
    μᶠ = p*(G(max(zs,zf))-G(zf))
    npᶜ = λ*(1-G(zp))*μᵖ*δ/((s+λ)*(μᵖ*δ + ξ*δ + μᶠ*ξ))
    np⁰ = ξ*μᵖ*δ/((s+λ)*(μᵖ*δ + ξ*δ + μᶠ*ξ))
    np = npᶜ + np⁰
    nf = μᶠ*ξ/(μᵖ*δ + ξ*δ + μᶠ*ξ)
    u = 1-np-nf
    v = u*θ

    z_bar_p = (λ*np*I(zp)+p*u*I(max(zs,zc)))/((s+λ)*np)
    z_bar_f = (λ*nf*Ez+p*u*(I(zf)-I(max(zs,zf))))/((δ+λ)*nf)
    z_bar = (np*z_bar_p+nf*z_bar_f)/(np+nf)

    W = z_bar_p*np + z_bar_f*nf + b*u - γ*v - λ*G(zp)*np*F
    W /= r
    if isnothing(index)
        @save_struct ss zp zc zs zf npᶜ np⁰ np nf θ rU ξ μᵖ μᶠ z_bar_p z_bar_f z_bar W u 
    else
        @save_struct_index ss index zp zc zs zf npᶜ np⁰ np nf θ rU ξ μᵖ μᶠ z_bar_p z_bar_f z_bar W u 
    end

end

"""
    compute_classic!(ss, sol, param)

Fills the `ClassicSS` structure with the classic-equilibrium thresholds and employments given basic equilibrium values `sol` and parameters `param`
"""
function compute_classic!(ss, sol, param; index=nothing)
    Macros.@load_struct param η r s λ α m b γ σ F δ
    @spec_funs

    zp = sol[1]
    θ = sol[2]

    q = m*θ^(-σ)
    rU = b+η*γ*θ/(1-η)
    zc = zp+(r+s+λ)*F
    ξ = s+λ*G(zp)
    p = m*θ^(1-σ)
    μ = p*(1-G(zc))
    nᶜ = λ*(1-G(zp))*μ/((ξ+μ)*(s+λ))
    n⁰ = ξ*μ/((ξ+μ)*(s+λ))
    n = μ/(ξ+μ)
    u = 1-n
    v = u*θ

    z_bar = (λ*n*I(zp)+p*u*I(zc))/((s+λ)*n)
    W = z_bar*n + b*u - γ*v - λ*G(zp)*n*F
    W /= r
    if isnothing(index)
        @save_struct ss zp zc n nᶜ n⁰ θ μ z_bar W ξ rU u
    else
        @save_struct_index ss index zp zc n nᶜ n⁰ θ μ z_bar W ξ rU u
    end

end

"""
    EqClassic!(R,x,param)

Replaces `R` with the values of the fundamental equation block of the dual equilibrium given fundamental variables `x` and parameters `param`.
"""
function EqClassic!(R, x, param)
    Macros.@load_struct param η r s λ α m b γ σ F δ
    @spec_funs

    zp = x[1]
    θ = max(x[2],0.)

    q = m*θ^(-σ)
    rU = b + η*γ*θ/(1-η)
    zc = zp + (r+s+λ)*F

    R[1] = zp+(r+s)*F+λ*I0(zp)/(r+s+λ)-rU
    R[2] = γ - q*(1-η)*I0(zc)/(r+s+λ)

end

"""
    EqDual!(R,x,param)

Replaces `R` with the values of the fundamental equation block of the dual equilibrium given fundamental variables `x` and parameters `param`. Used to compute the steady state in the dual model.
"""
function EqDual!(R, x, param)
    Macros.@load_struct param η r s λ α m b γ σ F δ
    @spec_funs

    zp = x[1]
    θ = max(x[2],0.)

    zc = zp+(r+s+λ)*F
    q = m*θ^(-σ)
    rU = b + η*γ*θ/(1-η)
    zf = rU - λ*(Ez-rU)/(r+δ)
    zs = (zp/(r+s+λ) + F - zf/(r+δ+λ))/(1/(r+s+λ) - 1/(r+δ+λ))

    R[1] = zp+(r+s)*F+λ*I0(zp)/(r+s+λ)-rU
    R[2] = γ - q*(1-η)*( I0(max(zs,zc))/(r+s+λ) + (I0(zf)-I0(max(zs,zf)))/(r+δ+λ) )

end

"""
    solve_F1!(R, x, param)

Replaces `R` with the values of the fundamental equation block of the dual equilibrium under the constraint ``z^f = z^* = z^c``. Used to compute the dual steady state under the constraint ``z^f = z^* = z^c``. 
"""
function solve_F1!(R, x, param)
    Macros.@load_struct param η r s λ α m b γ σ δ
    @spec_funs
   
    F1 = x[1]
    θ = max(x[2],0.)

    q = m*θ^(-σ)
    rU = b + η*γ*θ/(1-η)
    zc = rU - λ*(Ez-rU)/(r+δ)
    zp = zc - (r+s+λ)*F1

    R[1] = zp+(r+s)*F1+λ*I0(zp)/(r+s+λ)-rU
    R[2] = γ - q*(1-η)*I0(zc)/(r+s+λ)
end

"""
    solve_F2!(R, x, param)

Replaces `R` with the values of the fundamental equation block of the dual equilibrium under the constraint ``z^p = 0``. Used to compute the dual steady state under the constraint ``z^p = 0``. 
"""
function solve_F2!(R, x, param)
    Macros.@load_struct param η r s λ α m b γ σ δ
    @spec_funs

    F2 = x[1]
    θ = max(x[2],0.)

    q = m*θ^(-σ)
    rU = b + η*γ*θ/(1-η)
    zf = rU - λ*(Ez-rU)/(r+δ)
    zs = ( F2 - zf/(r+δ+λ))/(1/(r+s+λ) - 1/(r+δ+λ))

    R[1] = (r+s)*F2+λ*Ez/(r+s+λ)-rU
    R[2] = γ - q*(1-η)*( I0(zs)/(r+s+λ) + (I0(zf)-I0(zs))/(r+δ+λ) )
end

"""
    fill_AdjDual!(adj_dual, F, ss_old, param, x0)

Fills the `index` index of the fields of the `comp_stat` structure with labor market values when firing costs equal `F`, `param` contains the other parameters, `ss` contains the initial steady state and `x0` is the initial guess for solving the system `EqDual` derives.
"""
function fill_AdjDual!(adj_dual, F, ss_old, param, x0; index=nothing)

    Macros.@load_struct param η r s λ α m b γ σ δ
    @spec_funs

    param.F = F
    ss = DualSS()
    fill_DualSS!(ss, param, x0)
    
    @load_struct_cat old ss_old npᶜ np⁰ np nf zp zc zs
    @load_struct ss zp zc zs zf npᶜ np⁰ np nf θ rU ξ μᵖ μᶠ z_bar_p z_bar_f z_bar W u  

    @dual_to_dual

    τ_p = max_zeros(fzeros(t->abs((nᵖ(t)-np)/np_old)-0.01,0,500))
    τ_f = max_zeros(fzeros(t->abs((nᶠ(t)-nf)/nf_old)-0.01,0,500))

    if isnothing(index)
        @save_struct adj_dual τ_p τ_f
    else
        @save_struct_index adj_dual index τ_p τ_f
    end

end

"""
    fill_AdjClassic!(adj_classic, F, ss_old, param, x0)

Fills the fields of the `comp_stat` structure with labor market values when firing costs equal `F`, `param` contains the other parameters, `ss` contains the initial steady state and `x0` is the initial guess for solving the system `EqDual` derives.
"""
function fill_AdjClassic!(adj_classic, F, ss_old, param, x0; index=nothing)

    Macros.@load_struct param η r s λ α m b γ σ δ
    @spec_funs

    param.F = F
    ss = ClassicSS()
    fill_ClassicSS!(ss, param, x0)
    
    @load_struct_cat old ss_old nᶜ n⁰ n zp zc
    @load_struct ss zp zc n nᶜ n⁰ θ μ z_bar W ξ rU
  
    @classic_to_classic

    τ = max_zeros(fzeros(t->abs((nᵖ(t)-n)/n_old)-0.01,0,500))

    if isnothing(index)
        @save_struct adj_classic τ
    else
        @save_struct_index adj_classic index τ
    end

end

"""
    fill_ClassicSS!(ss, param, x0)

Fills the `ss` structure with the steady-state values of the classic model when `param` contains parameters and `x0` is the initial value to solve the system `EqClassic` derives.
"""
function fill_ClassicSS!(ss, param, x0; index=nothing)
    R = zeros(2)
    sol = nlsolve((R,x)->EqClassic!(R,x,param), x0, show_trace=false)
    compute_classic!(ss, sol.zero, param; index=index)
end

"""
    fill_DualSS!(ss, param, x0)

Fills the `ss` structure with the steady-state values of the dual model when `param` contains parameters and `x0` is the initial value to solve the system `EqDual` derives.
"""
function fill_DualSS!(ss, param, x0; index=nothing)
    R = zeros(2)
    sol = nlsolve((R,x)->EqDual!(R,x,param), x0, show_trace=false)
    compute_dual!(ss, sol.zero, param; index=index)
end

"""
    plot_cmpstat_bw(param, f; x0 = ..., x1 = ...)

Plots in black and white various comparative statics using parameters `param` and saving the figures with the suffix `f`. `x0` and `x1` contain starting values to solve the involved systems of equations.
"""
function plot_cmpstat_bw(param, f; x0 = [0.9, 1.3], x1 = [1.33, 1.3])

    @load_struct param F
    classic = ClassicSS()
    fill_ClassicSS!(classic, param, x0)
    dual = DualSS()
    fill_DualSS!(dual, param, x0)

    sol = nlsolve((R,x)->solve_F1!(R,x,param),x1,show_trace=true)
    F1 = sol.zero[1]
    sol = nlsolve((R,x)->solve_F2!(R,x,param),x1,show_trace=true)
    F2 = sol.zero[1]

    N_F = 100

    dual_new = DualSS()
    classic_new = ClassicSS()
    adj_dual = AdjDual()
    adj_classic = AdjClassic()
    init!(dual_new, N_F)
    init!(classic_new, N_F)
    init!(adj_dual, N_F)
    init!(adj_classic, N_F)

    grid_F = range(0.8*F1, stop=1.1*F, length=N_F)
    for i in 1:N_F
        param.F = grid_F[i]
        fill_ClassicSS!(classic_new, param, x0; index=i)
        fill_DualSS!(dual_new, param, x0; index=i)
        fill_AdjClassic!(adj_classic, grid_F[i], classic, param, x0; index=i)
        fill_AdjDual!(adj_dual, grid_F[i], dual, param, x0; index=i)
    end

    cd("figures")

    # Transition probabilities
    @load_struct dual_new ξ μᵖ μᶠ np nf u z_bar_p z_bar_f z_bar W 
    @load_struct adj_dual τ_p τ_f
    @load_struct_cat c classic_new ξ μ n u z_bar W 
    @load_struct_cat c adj_classic τ 
    p1 = plot(grid_F, ξ,lab=L"$\xi$ - Dual",color=:black,line=(:solid, 1))
    plot!(grid_F, μᵖ,lab=L"$\mu^p$ - Dual",color=:black,line=(:dash, 2))
    plot!(grid_F, μᶠ,lab=L"$\mu^f$ - Dual",color=:black,line=(:dot, 1))
    plot!(grid_F, ξ_c,lab=L"$\xi$ - Classic",color=:black,line=(:dashdot, 1))
    plot!(grid_F, μ_c,lab=L"$\mu$ - Classic",color=:black,line=(:dash, 1), legend=:right)
    plot!([F1,F], seriestype=:vline, lab="",color=:black,line=(:dash, 1))
    xlabel!("Firing costs")
    ylabel!("Probabilities")
    xticks!([F1,F],[L"$\widehat{F}$",L"$F$"])
    xgrid!(:off)
    ygrid!(:off)
    savefig("transitions_"*f*".pdf")

    # Employment

    p2 = plot(grid_F, np,lab=L"$n^p$ - Dual",color=:black,line=(:dash, 2))
    plot!(grid_F, nf,lab=L"$n^f$ - Dual",color=:black,line=(:dot, 2))
    plot!(grid_F, u,lab=L"$u$ - Dual",color=:black,line=(:solid, 2))
    plot!(grid_F, n_c,lab=L"$n^p$ - Classic",color=:black,line=(:dash, 1))
    plot!(grid_F, u_c,lab=L"$u$ - Classic",color=:black,line=(:dashdot, 1), legend=:right)
    plot!([F1,F], seriestype=:vline, lab="",color=:black,line=(:dash, 1))
    xlabel!("Firing costs")
    ylabel!("Employment")
    xticks!([F1,F],[L"$\widehat{F}$",L"$F$"])
    xgrid!(:off)
    ygrid!(:off)
    savefig("employment_"*f*".pdf")

    #Probabilities and Employment
    plot(p1,p2,layout=(1,2))
    savefig("prob_emp_"*f*".pdf")

    # Productivity
    p4 = plot(grid_F, z_bar_p,lab=L"$\overline{z}^p$ - Dual",color=:black,line=(:dash, 2))
    plot!(grid_F, z_bar_f,lab=L"$\overline{z}^f$ - Dual",color=:black,line=(:dot, 2))
    plot!(grid_F, z_bar,lab=L"$\overline{z}$ - Dual",color=:black,line=(:solid, 2))
    plot!(grid_F, z_bar_c,lab=L"$\overline{z}$ - Classic",color=:black,line=(:dashdot, 1))
    plot!([F1,F], seriestype=:vline, lab="",color=:black,line=(:dash, 1))
    xlabel!("Firing costs")
    ylabel!("Productivity")
    xticks!([F1,F],[L"$\widehat{F}$",L"$F$"])
    xgrid!(:off)
    ygrid!(:off)
    savefig("productivity_"*f*".pdf")

    # Welfare
    p3 = plot(grid_F, W,lab=L"Dual",color=:black,line=(:solid, 2))
    plot!(grid_F, W_c,lab=L"Classic",color=:black,line=(:dash, 2))
    plot!([F1,F], seriestype=:vline, lab="",color=:black,line=(:dash, 1))
    xlabel!("Firing costs")
    ylabel!("Welfare")
    xticks!([F1,F],[L"$\widehat{F}$",L"$F$"])
    xgrid!(:off)
    ygrid!(:off)
    savefig("welfare_"*f*".pdf")

    #Welfare and Productivity
    plot(p3,p4,layout=(1,2))
    savefig("wel_prod_"*f*".pdf")

    # Adjustment
    plot(grid_F, τ_p,lab=L"$n^p$ - Dual",color=:black,line=(:solid, 2))
    plot!(grid_F, τ_f,lab=L"$n^f$ - Dual",color=:black,line=(:dot, 2))
    plot!(grid_F, τ_c,lab=L"$n$ - Classic",color=:black,line=(:dash, 2))
    plot!([F1,F], seriestype=:vline, lab="",color=:black,line=(:dash, 1))
    xlabel!("Firing costs")
    ylabel!(L"$\tau_{99\%}$")
    xticks!([F1,F],[L"$\widehat{F}$",L"$F$"])
    xgrid!(:off)
    ygrid!(:off)
    savefig("adjustment_"*f*".pdf")

    cd("..")

end

"""
    plot_cmpstat_colors(param, f; x0 = ..., x1 = ...)

Plots in colors various comparative statics using parameters `param` and saving the figures with the suffix `f`. `x0` and `x1` contain starting values to solve the involved systems of equations.
"""
function plot_cmpstat_colors(param, f; x0 = [0.9, 1.3], x1 = [1.33, 1.3])

    @load_struct param F
    classic = ClassicSS()
    fill_ClassicSS!(classic, param, x0)
    dual = DualSS()
    fill_DualSS!(dual, param, x0)

    sol = nlsolve((R,x)->solve_F1!(R,x,param),x1,show_trace=false)
    F1 = sol.zero[1]
    sol = nlsolve((R,x)->solve_F2!(R,x,param),x1,show_trace=false)
    F2 = sol.zero[1]

    N_F = 100

    dual_new = DualSS()
    classic_new = ClassicSS()
    adj_dual = AdjDual()
    adj_classic = AdjClassic()
    init!(dual_new, N_F)
    init!(classic_new, N_F)
    init!(adj_dual, N_F)
    init!(adj_classic, N_F)

    grid_F = range(0.8*F1, stop=1.1*F, length=N_F)
    for i in 1:N_F
        param.F = grid_F[i]
        fill_ClassicSS!(classic_new, param, x0; index=i)
        fill_DualSS!(dual_new, param, x0; index=i)
        fill_AdjClassic!(adj_classic, grid_F[i], classic, param, x0; index=i)
        fill_AdjDual!(adj_dual, grid_F[i], dual, param, x0; index=i)
    end

    cd("figures")

    # Transition probabilities
    @load_struct dual_new ξ μᵖ μᶠ np nf u z_bar_p z_bar_f z_bar W 
    @load_struct adj_dual τ_p τ_f
    @load_struct_cat c classic_new ξ μ n u z_bar W 
    @load_struct_cat c adj_classic τ 
    p1 = plot(grid_F, ξ,lab="Open-ended JD",color=:red, line=(:dash, 1))
    plot!(grid_F, μᵖ,lab="Open-ended JC",color=:red, line=(:solid, 1))
    plot!(grid_F, μᶠ,lab="Fixed-term JC",color=:blue, line=(:solid, 1))
    plot!([F1,F], seriestype=:vline, lab="",color=:black, line=(:dash, 1))
    xlabel!("Firing costs")
    ylabel!("Probabilities")
    xticks!([F1,F],[L"$\widehat{F}$",L"$F$"])
    xgrid!(:off)
    ygrid!(:off)

    savefig("transitions_color_"*f*".pdf")

    # Employment
    p2 = plot(grid_F, nf, lab="Fixed-term",color=:blue,line=(:solid, 1))
    plot!([F1,F], seriestype=:vline, lab="",color=:black,line=(:dash, 1))
    xlabel!("Firing costs")
    ylabel!("Employment")
    xticks!([F1,F],[L"$\widehat{F}$",L"$F$"])
    xgrid!(:off)
    ygrid!(:off)
    savefig("ft_emp_color_"*f*".pdf")
    p3 = plot(grid_F, np,lab="Open-ended",color=:red,line=(:solid, 1))
    plot!(grid_F, u,lab="Non-employment",color=:black,line=(:solid, 1))
    plot!([F1,F], seriestype=:vline, lab="",color=:black,line=(:dash, 1))
    xlabel!("Firing costs")
    ylabel!("Employment")
    xticks!([F1,F],[L"$\widehat{F}$",L"$F$"])
    xgrid!(:off)
    ygrid!(:off)
    savefig("oe_emp_color_"*f*".pdf")

    # Productivity
    p4 = plot(grid_F, z_bar_p,lab=L"$\overline{z}^p$ - Dual",color=:red,line=(:solid, 1))
    plot!(grid_F, z_bar_f,lab=L"$\overline{z}^f$ - Dual",color=:blue,line=(:solid, 1))
    plot!(grid_F, z_bar,lab=L"$\overline{z}$ - Dual",color=:black,line=(:solid, 1))
    plot!(grid_F, z_bar_c,lab=L"$\overline{z}$ - Classic",color=:black,line=(:dash, 1))
    plot!([F1,F], seriestype=:vline, lab="",color=:black,line=(:dash, 1))
    xlabel!("Firing costs")
    ylabel!("Productivity")
    xticks!([F1,F],[L"$\widehat{F}$",L"$F$"])
    xgrid!(:off)
    ygrid!(:off)
    savefig("productivity_color_"*f*".pdf")

    # Welfare
    p5 = plot(grid_F, W,lab=L"Dual",color=:black,line=(:solid, 2))
    plot!(grid_F, W_c,lab=L"Classic",color=:black,line=(:dash, 2))
    plot!([F1,F], seriestype=:vline, lab="",color=:black,line=(:dash, 1))
    xlabel!("Firing costs")
    ylabel!("Welfare")
    xticks!([F1,F],[L"$\widehat{F}$",L"$F$"])
    xgrid!(:off)
    ygrid!(:off)
    savefig("welfare_color_"*f*".pdf")

    # Adjustment
    p6 = plot(grid_F, τ_p,lab="Open-ended",color=:red,line=(:solid, 1))
    plot!(grid_F, τ_f,lab="Fixed-term",color=:blue,line=(:solid, 1))
    # plot!(grid_F, τ_c,lab=L"$n$ - Classic",color=:black,line=(:dash, 2))
    plot!([F1,F], seriestype=:vline, lab="",color=:black,line=(:dash, 1))
    xlabel!("Firing costs")
    ylabel!(L"$\tau_{99\%}$")
    xticks!([F1,F],[L"$\widehat{F}$",L"$F$"])
    xgrid!(:off)
    ygrid!(:off)
    savefig("adj_colors_"*f*".pdf")

    cd("..")
end

"""
    dual_to_classic(F_ref, param, f; x0 = ... )

Plots the transition between a dual and a classic equilibrium when firing costs unexpectedly move from their value stored in `param` to `F_ref`. The resulting figures will have the suffix `f`. 
"""
function dual_to_classic(F_ref, param, f; x0 = [0.9, 1.3])

    dual = DualSS()
    fill_DualSS!(dual, param, x0)
    param_ref = deepcopy(param)
    param_ref.F = F_ref
    classic = ClassicSS()
    fill_ClassicSS!(classic, param_ref, x0)

    @load_struct_cat old dual npᶜ np⁰ np nf zp zc zs u 
    @load_struct classic nᶜ n⁰ ξ μ zp zc 
    Macros.@load_struct param_ref δ α s λ
    @spec_funs

    @dual_to_classic
    cd("figures")

    t = range(0.,stop=24,length=1000)
    p1 = plot(t,vcat(np_old, nᵖ.(t[2:end])),lab=L"$n^p$",color=:black,line=(:solid, 1))
    plot!(t,vcat(nf_old, nᶠ.(t[2:end])),lab=L"$n^f$",color=:black,line=(:dot, 2))
    plot!(t,vcat(u_old, u.(t[2:end])),lab=L"$u$",color=:black,line=(:dash, 2), legend=:right)
    xlabel!("Time (months)")
    ylabel!("Employment")
    xgrid!(:off)
    ygrid!(:off)

    savefig("dual_to_classic_"*f*".pdf")

    p2 = plot(t,vcat(np_old, nᵖ.(t[2:end])),lab="Open-ended",color=:red,line=(:solid, 1))
    plot!(t,vcat(nf_old, nᶠ.(t[2:end])),lab="Fixed-term",color=:blue,line=(:dot, 1))
    plot!(t,vcat(u_old, u.(t[2:end])),lab="Non-employment",color=:black,line=(:dash, 1), legend=:right)
    xlabel!("Time (months)")
    ylabel!("Employment")
    xgrid!(:off)
    ygrid!(:off)

    savefig("dual_to_classic_color_"*f*".pdf")

    cd("..")
    
    return p1


end


"""
    classic_to_dual(F_ref, param, f; x0 = ... )

Plots the transition between a classic and a dual equilibrium when firing costs unexpectedly move from their value stored in `param` to `F_ref`. The resulting figures will have the suffix `f`. 
"""
function classic_to_dual(param, f; x0 = [0.9, 1.3])

    classic = ClassicSS()
    fill_ClassicSS!(classic, param, x0)
    dual = DualSS()
    fill_DualSS!(dual, param, x0)

    @load_struct_cat old classic nᶜ n⁰ n zp zc
    npᶜ_old = nᶜ_old
    np⁰_old = n⁰_old
    np_old = n_old
    nf_old = 0.
    u_old = 1-np_old-nf_old
    @load_struct dual npᶜ np⁰ nf zc zs zf zp np ξ μᶠ μᵖ
    Macros.@load_struct param δ λ s α
    @spec_funs

    @classic_to_dual

    cd("figures")

    t = range(0,stop=60,length=1000)
    p1 = plot(t,vcat(np_old,nᵖ.(t[2:end])),lab=L"$n^p$",color=:black,line=(:solid, 1))
    plot!(t,vcat(nf_old, nᶠ.(t[2:end])),lab=L"$n^f$",color=:black,line=(:dot, 2))
    plot!(t,vcat(u_old,u.(t[2:end])),lab=L"$u$",color=:black,line=(:dash, 2), legend=:right)
    xlabel!("Time (months)")
    xgrid!(:off)
    ygrid!(:off)

    savefig("classic_to_dual_"*f*".pdf")

    cd("..")

    return p1
end

function plot_transitions(F_ref, param, f; x0 = [0.9, 1.3])
    p1 = dual_to_classic(F_ref, param, f; x0 = x0)
    p2 = classic_to_dual(param,f; x0 = x0)
    plot(p1,p2,layout=(1,2))
    cd("figures")
    savefig("dyn_reforms_"*f*".pdf")
    cd("..")
end

end