module Calibration

using Roots
using NLsolve
using Distributions
using Macros
using CustomTypes

export calibrate_y, calibrate_m, calibrate_ssξ, calibrate_u, calibrate_c

@def cal_common begin
    μ_z = 0.
    Ez(α) = exp(μ_z + α^2/2)
    G(x,α) = cdf(LogNormal(μ_z,α),x)
    I(x,α) = exp(μ_z+α^2/2)*cdf(Normal(0,1),((μ_z+α^2-log(max(x,0)))/α))
    I0(x,α) = I(x,α) - x*(1-G(x,α))

    zp(α) = fzero( zp -> (ξ-s)/λ - G(zp,α) , 0, Ez(α)+3*α )
    rU(α) = zp(α)+(r+s)*F+λ*I0(zp(α),α)/(r+s+λ)
    zf(α) = rU(α) - λ*(Ez(α)-rU(α))/(r+δ)
    zs(α) = (zp(α)/(r+s+λ) + F - zf(α)/(r+δ+λ))/(1/(r+s+λ) - 1/(r+δ+λ))

    share_tc_jc(α) = (G(zs(α),α)-G(zf(α),α))/(1-G(zf(α),α))

    α = fzero(α->share_tc_jc(α) - μᶠsμ, 0.01, 1.)

    Ez_cal = Ez(α)
    G_cal(x) = G(x,α)
    I_cal(x) = I(x,α)
    I0_cal(x) = I0(x,α)

    zp_cal = zp(α)
    rU_cal = rU(α)
    zf_cal = zf(α)
    zs_cal = zs(α)

    p = μᵖ/(1-G_cal(zs_cal))
    γ = q*(1-η)*( I0_cal(zs_cal)/(r+s+λ) + (I0_cal(zf_cal)-I0_cal(zs_cal))/(r+δ+λ) )
    θ = p/q
    b = rU_cal - η*γ*θ/(1-η)
    m = q*θ^σ

    wf_bar = η*( λ*Ez_cal + δ*(I_cal(zf_cal) - I_cal(zs_cal))/(G_cal(zs_cal)-G_cal(zf_cal)) )/(δ+λ) + (1-η)*rU_cal

    E_wp = η*(I_cal(zp_cal)/(1-G_cal(zp_cal)) + (r+s)*F) + (1-η)*rU_cal
    E_wp_0 = η*(I_cal(zs_cal)/(1-G_cal(zs_cal)) - λ*F) + (1-η)*rU_cal
    wp_bar = (λ*(1-G_cal(zp_cal))*E_wp + ξ*E_wp_0)/(s+λ)


end

"""
    calibrate_y()

Returns a dictionary linking parameters and their values when the average time interval between two productivity shocks is one year
"""
function calibrate_y() 
    η = 0.6 #Worker's bargaining power
    σ = 0.6 #Elasticity of the matching function w.r.t vacancies
    r = 0.05/12 #Monthly interest rate
    δ = 1/1.5 #Probability of temporary jobs' destruction
    rho_F = 4/3 #Firing-cost-to-average-permanent-contract-wage ratio
    λ = 1/12 #Probability of occurence of productivity shocks

    u = 0.26 #Unemployment
    nfsn = 0.12 #Temporary-to-total-employment ratio
    q = 1-(1-0.7)^(1/3) #Firms' worker-meeting probability
    ssξ = 0.7 #Exogenous to total permanent splits
    μᶠsμ = 0.83 #Share of temporary jobs in job creation

    n = 1-u
    nf = nfsn*n
    np = n-nf

    μᶠ = δ*nf/u
    μ = μᶠ/μᶠsμ
    μᵖ = μ-μᶠ

    ξ = μᵖ*u/np
    s = ssξ*ξ

    n = 1-u
    nf = nfsn*n
    np = n-nf

    μᶠ = δ*nf/u
    μ = μᶠ/μᶠsμ
    μᵖ = μ-μᶠ

    ξ = μᵖ*u/np
    s = ssξ*ξ
    
    function Eq(F)
        @cal_common
        return (F/wp_bar)-rho_F
    end

    F = fzero(Eq,1.2,2.)

    @cal_common

    cal = Parameters()
    Macros.@save_struct cal η r s λ α m b γ σ F δ
    return cal

end

"""

    calibrate_m()

Returns a dictionary linking parameters and their values when the average time interval between two productivity shocks is one month.
"""
function calibrate_m() 

    η = 0.6 #Worker's bargaining power
    σ = 0.6 #Elasticity of the matching function w.r.t vacancies
    r = 0.05/12 #Monthly interest rate
    δ = 1/1.5 #Probability of temporary jobs' destruction
    rho_F = 4/3 #Firing-cost-to-average-permanent-contract-wage ratio
    λ = 1 #Probability of occurence of productivity shocks

    u = 0.26 #Unemployment
    nfsn = 0.12 #Temporary-to-total-employment ratio
    q = 1-(1-0.7)^(1/3) #Firms' worker-meeting probability
    ssξ = 0.7 #Exogenous to total permanent splits
    μᶠsμ = 0.83 #Share of temporary jobs in job creation

    n = 1-u
    nf = nfsn*n
    np = n-nf

    μᶠ = δ*nf/u
    μ = μᶠ/μᶠsμ
    μᵖ = μ-μᶠ

    ξ = μᵖ*u/np
    s = ssξ*ξ

    n = 1-u
    nf = nfsn*n
    np = n-nf

    μᶠ = δ*nf/u
    μ = μᶠ/μᶠsμ
    μᵖ = μ-μᶠ

    ξ = μᵖ*u/np
    s = ssξ*ξ
    
    function Eq(F)
        @cal_common
        return (F/wp_bar)-rho_F
    end

    F = fzero(Eq,1.2,2.)

    @cal_common
 
    cal = Parameters()
    Macros.@save_struct cal η r s λ α m b γ σ F δ
    return cal
   
end

"""
    calibrate_c()

Returns a dictionary linking parameters and their values when the average time interval between two productivity shocks is three years.
"""
function calibrate_c() 

    η = 0.6 #Worker's bargaining power
    σ = 0.6 #Elasticity of the matching function w.r.t vacancies
    r = 0.05/12 #Monthly interest rate
    δ = 1/1.5 #Probability of temporary jobs' destruction
    rho_F = 4/3 #Firing-cost-to-average-permanent-contract-wage ratio
    λ = 1/36 #Probability of occurence of productivity shocks

    u = 0.26 #Unemployment
    nfsn = 0.12 #Temporary-to-total-employment ratio
    q = 1-(1-0.7)^(1/3) #Firms' worker-meeting probability
    ssξ = 0.7 #Exogenous to total permanent splits
    μᶠsμ = 0.83 #Share of temporary jobs in job creation

    n = 1-u
    nf = nfsn*n
    np = n-nf

    μᶠ = δ*nf/u
    μ = μᶠ/μᶠsμ
    μᵖ = μ-μᶠ

    ξ = μᵖ*u/np
    s = ssξ*ξ

    n = 1-u
    nf = nfsn*n
    np = n-nf

    μᶠ = δ*nf/u
    μ = μᶠ/μᶠsμ
    μᵖ = μ-μᶠ

    ξ = μᵖ*u/np
    s = ssξ*ξ
    
    function Eq(F)
        @cal_common
        return (F/wp_bar)-rho_F
    end

    F = fzero(Eq,1.2,2.)

    @cal_common
    
    cal = Parameters()
    Macros.@save_struct cal η r s λ α m b γ σ F δ
    return cal

end

"""
    calibrate_ssξ()

Returns a dictionary linking parameters and their values when the average time interval between two productivity shocks is one year and s = 0.
"""
function calibrate_ssξ() 

    η = 0.6 #Worker's bargaining power
    σ = 0.6 #Elasticity of the matching function w.r.t vacancies
    r = 0.05/12 #Monthly interest rate
    δ = 1/1.5 #Probability of temporary jobs' destruction
    rho_F = 4/3 #Firing-cost-to-average-permanent-contract-wage ratio
    λ = 1/12 #Probability of occurence of productivity shocks

    u = 0.26 #Unemployment
    nfsn = 0.12 #Temporary-to-total-employment ratio
    q = 1-(1-0.7)^(1/3) #Firms' worker-meeting probability
    ssξ = 0. #Exogenous to total permanent splits
    μᶠsμ = 0.83 #Share of temporary jobs in job creation

    n = 1-u
    nf = nfsn*n
    np = n-nf

    μᶠ = δ*nf/u
    μ = μᶠ/μᶠsμ
    μᵖ = μ-μᶠ

    ξ = μᵖ*u/np
    s = ssξ*ξ

    n = 1-u
    nf = nfsn*n
    np = n-nf

    μᶠ = δ*nf/u
    μ = μᶠ/μᶠsμ
    μᵖ = μ-μᶠ

    ξ = μᵖ*u/np
    s = ssξ*ξ
    
    function Eq(F)
        @cal_common
        return (F/wp_bar)-rho_F
    end

    F = fzero(Eq,1.2,2.)

    @cal_common
          
    cal = Parameters()
    Macros.@save_struct cal η r s λ α m b γ σ F δ
    return cal

end

"""
    calibrate_u()

Returns a dictionary linking parameters and their values when the average time interval between two productivity shocks is one year and the job seekers' pool does not include inactive workers.
"""
function calibrate_u() 
    η = 0.6 #Worker's bargaining power
    σ = 0.6 #Elasticity of the matching function w.r.t vacancies
    r = 0.05/12 #Monthly interest rate
    δ = 1/1.5 #Probability of temporary jobs' destruction
    rho_F = 4/3 #Firing-cost-to-average-permanent-contract-wage ratio
    λ = 1/12 #Probability of occurence of productivity shocks

    u = 0.1 #Unemployment
    nfsn = 0.12 #Temporary-to-total-employment ratio
    q = 1-(1-0.7)^(1/3) #Firms' worker-meeting probability
    ssξ = 0.7 #Exogenous to total permanent splits
    μᶠsμ = 0.83 #Share of temporary jobs in job creation

    n = 1-u
    nf = nfsn*n
    np = n-nf

    μᶠ = δ*nf/u
    μ = μᶠ/μᶠsμ
    μᵖ = μ-μᶠ

    ξ = μᵖ*u/np
    s = ssξ*ξ

    n = 1-u
    nf = nfsn*n
    np = n-nf

    μᶠ = δ*nf/u
    μ = μᶠ/μᶠsμ
    μᵖ = μ-μᶠ

    ξ = μᵖ*u/np
    s = ssξ*ξ
    
    function Eq(F)
        @cal_common
        return (F/wp_bar)-rho_F
    end

    F = fzero(Eq,1.2,2.)

    @cal_common
        
    cal = Parameters()
    Macros.@save_struct cal η r s λ α m b γ σ F δ
    return cal

end
end