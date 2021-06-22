module CustomTypes

export Parameters, ClassicSS, DualSS, AdjDual, AdjClassic, RegUncSS

mutable struct Parameters
    η
    r
    s
    λ
    α
    m
    σ
    γ
    F
    δ
    b
    ϵ
    Parameters() = new()
end

mutable struct ClassicSS
    zp
    zc
    n
    nᶜ
    n⁰
    u
    θ
    μ
    z_bar
    W
    ξ
    rU
    ClassicSS() = new()
end

mutable struct DualSS
    zp
    zc
    zs
    zf
    npᶜ
    np⁰
    np
    nf
    θ
    rU
    ξ
    μᵖ
    μᶠ
    z_bar_p
    z_bar_f
    z_bar
    W
    u
    DualSS() = new()
end

mutable struct AdjDual
    τ_p
    τ_f
    AdjDual() = new()
end

mutable struct AdjClassic
    τ
    AdjClassic() = new()
end

mutable struct RegUncSS
    F
    zp
    zc
    zs
    zf
    npᶜ
    np⁰
    np
    nf
    u
    q
    θ
    ξ
    μᵖ
    μᶠ
    z_bar_p
    z_bar_f
    z_bar
    W
    Wpᶜ
    Wp⁰
    Wp
    flag
    RegUncSS() = new()
end

end