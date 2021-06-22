import Pkg;
Pkg.activate(".")
# %%
path = pwd() * "/functions"
push!(LOAD_PATH,path)
# %%
using PyCall
sym = pyimport("sympy")
# %%
p, ξ, δ, ω, μᵖ, μᶠ, c₀ᶜ , c₀⁰ , n₀ᶠ, Nᶜ , N⁰ , Nᶠ, α₁, α₂, nᶜ, ρ₀ᶜ, ρ₁ᶜ, ρ₂ᶜ, n⁰, ρ₀⁰, ρ₁⁰, ρ₂⁰, nᶠ, ρ₁ᶠ, ρ₂ᶠ, E₀ᶜ, Eᶜ, E₀⁰, E⁰, α₀, β₀, Wᶜ, W⁰, n₀ᶜ, n₀⁰, ν₀ᶜ, ν₁ᶜ, ν₂ᶜ, ν₀⁰, ν₁⁰, ν₂⁰, W₀ᶜ, W₀⁰, W∞ᶜ, W∞⁰ = sym.symbols("p, ξ, δ, ω , μᵖ , μᶠ, c₀ᶜ , c₀⁰ , n₀ᶠ, Nᶜ , N⁰ , Nᶠ, α₁, α₂, nᶜ, ρ₀ᶜ, ρ₁ᶜ, ρ₂ᶜ, n⁰, ρ₀⁰, ρ₁⁰, ρ₂⁰, nᶠ, ρ₁ᶠ, ρ₂ᶠ, E₀ᶜ, Eᶜ, E₀⁰, E⁰, α₀, β₀, Wᶜ, W⁰, n₀ᶜ, n₀⁰, ν₀ᶜ, ν₁ᶜ, ν₂ᶜ, ν₀⁰, ν₁⁰, ν₂⁰, W₀ᶜ, W₀⁰, W∞ᶜ, W∞⁰")
# %%
# ====================================================== #
# ================= Dual to dual ======================= #
# ====================================================== #
system = [(p+ξ)*Nᶜ + (ω+ξ)*N⁰ - c₀ᶜ, μᵖ*Nᶜ + (p + μᵖ - ω)*N⁰ + μᵖ*Nᶠ - ( μᵖ/p + c₀⁰ ), μᶠ*Nᶜ + μᶠ*N⁰ + (p+δ+μᶠ)*Nᶠ - (μᶠ/p + n₀ᶠ) ]
sol = sym.linsolve(system, [Nᶜ, N⁰, Nᶠ])
# %%
sol = sym.factor(sol, p)
# %%
Nᶜ = sol.args[1][1]
N⁰ = sol.args[1][2]
Nᶠ = sol.args[1][3]
# %%
top_Nᶜ, bot_Nᶜ = sym.fraction(Nᶜ)
top_Nᶜ =  sym.Poly(top_Nᶜ,p)
for c in top_Nᶜ.coeffs()
    println(sym.factor(c, [c₀ᶜ, c₀⁰, n₀ᶠ]))
end
# %%
top_N⁰, bot_N⁰ = sym.fraction(N⁰)
top_N⁰ =  sym.Poly(top_N⁰,p)
for c in top_N⁰.coeffs()
    println(sym.factor(c, [c₀ᶜ, c₀⁰, n₀ᶠ]))
end
# %%
top_Nᶠ, bot_Nᶠ = sym.fraction(Nᶠ)
top_Nᶠ =  sym.Poly(top_Nᶠ,p)
for c in top_Nᶠ.coeffs()
    println(sym.factor(c, [c₀ᶜ, c₀⁰, n₀ᶠ]))
end
# Laplace transforms seem correct !
# %%
# ======================= Nᶜ ========================= #
bot_Nᶜ = sym.Poly(bot_Nᶜ, p)
roots = sym.roots(bot_Nᶜ)
# %%
expr = sym.cancel(top_Nᶜ*p/(p*(p-ω)*(p-α₁)*(p-α₂)))
expr = expr.subs(p , 0)
# nᶜ : Checks out!
# %%
expr = sym.factor(sym.cancel(top_Nᶜ/(p*(p-ω)*(p-α₁)*(p-α₂)) - nᶜ/p),p)
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(0)
expr = sym.cancel(top*(p-ω)/bot)
expr = sym.factor(expr.subs(p , ω),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ])
# ρ₀ᶜ : Checks out!
# %%
expr = sym.cancel(top/bot - ρ₀ᶜ/(p-ω))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(ω)
expr = sym.cancel(top*(p-α₁)/bot)
expr = sym.factor(expr.subs(p , α₁),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ, ρ₀ᶜ])
# ρ₁ᶜ : Checks out !
# %%
expr = sym.cancel(top/bot - ρ₁ᶜ/(p-α₁))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(α₁)
expr = sym.cancel(top*(p-α₂)/bot)
expr = sym.factor(expr.subs(p , α₂),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ, ρ₀ᶜ, ρ₁ᶜ])
# %%
# ======================= N⁰ ========================= #
bot_N⁰ = sym.Poly(bot_N⁰, p)
roots = sym.roots(bot_N⁰)
# %%
expr = sym.cancel(top_N⁰*p/(p*(p-ω)*(p-α₁)*(p-α₂)))
expr = expr.subs(p , 0)
# n⁰ : Checks out!
# %%
expr = sym.factor(sym.cancel(top_N⁰/(p*(p-ω)*(p-α₁)*(p-α₂)) - n⁰/p),p)
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(0)
expr = sym.cancel(top*(p-ω)/bot)
expr = sym.factor(expr.subs(p , ω),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ])
# ρ₀⁰ : Checks out!
# %%
expr = sym.cancel(top/bot - ρ₀⁰/(p-ω))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(ω)
expr = sym.cancel(top*(p-α₁)/bot)
expr = sym.factor(expr.subs(p , α₁),[c₀ᶜ, c₀⁰, n₀ᶠ, n⁰, ρ₀⁰])
# ρ₁⁰ : Checks out !
# %%
expr = sym.cancel(top/bot - ρ₁⁰/(p-α₁))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(α₁)
expr = sym.cancel(top*(p-α₂)/bot)
expr = sym.factor(expr.subs(p , α₂),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ, ρ₀ᶜ, ρ₁ᶜ])
# ρ₂⁰ : Checks out !
# %%
# ======================= Nᶠ ========================= #
bot_Nᶠ = sym.Poly(bot_Nᶠ, p)
roots = sym.roots(bot_Nᶠ)
# %%
expr = sym.cancel(top_Nᶠ*p/(p*(p-α₁)*(p-α₂)))
expr = expr.subs(p , 0)
# nᶠ : Checks out!
# %%
expr = sym.factor(sym.cancel(top_Nᶠ/(p*(p-α₁)*(p-α₂)) - nᶠ/p),p)
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(0)
expr = sym.cancel(top*(p-α₁)/bot)
expr = sym.factor(expr.subs(p, α₁),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ])
# ρ₁ᶠ : Checks out!
# %%
expr = sym.cancel(top/bot - ρ₁ᶠ/(p-α₁))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(α₁)
expr = sym.cancel(top*(p-α₂)/bot)
expr = sym.factor(expr.subs(p , α₂),[c₀ᶜ, c₀⁰, n₀ᶠ, n⁰, ρ₀⁰])
# ρ₂ᶠ : Checks out!
# %%
# ====================================================== #
# ================= Dual to classic ==================== #
# ====================================================== #
system = [(p+ξ)*Nᶜ + (ω+ξ)*N⁰ - c₀ᶜ, μᵖ*Nᶜ + (p + μᵖ - ω)*N⁰ + μᵖ*Nᶠ - ( μᵖ/p + c₀⁰ ), (p+δ)*Nᶠ - n₀ᶠ]
sol = sym.linsolve(system, [Nᶜ, N⁰, Nᶠ])
# %%
sol = sym.factor(sol, p)
# %%
Nᶜ = sol.args[1][1]
N⁰ = sol.args[1][2]
Nᶠ = sol.args[1][3]
# %%
top_Nᶜ, bot_Nᶜ = sym.fraction(Nᶜ)
top_Nᶜ =  sym.Poly(top_Nᶜ,p)
for c in top_Nᶜ.coeffs()
    println(sym.factor(c, [c₀ᶜ, c₀⁰, n₀ᶠ]))
end
# %%
top_N⁰, bot_N⁰ = sym.fraction(N⁰)
top_N⁰ =  sym.Poly(top_N⁰,p)
for c in top_N⁰.coeffs()
    println(sym.factor(c, [c₀ᶜ, c₀⁰, n₀ᶠ]))
end
# %%
top_Nᶠ, bot_Nᶠ = sym.fraction(Nᶠ)
top_Nᶠ =  sym.Poly(top_Nᶠ,p)
for c in top_Nᶠ.coeffs()
    println(sym.factor(c, [c₀ᶜ, c₀⁰, n₀ᶠ]))
end
# Laplace transforms seem correct !
# %%
# ======================= Nᶜ ========================= #
bot_Nᶜ = sym.Poly(bot_Nᶜ, p)
roots = sym.roots(bot_Nᶜ)
# %%
expr = sym.cancel(top_Nᶜ*p/(p*(p-ω)*(p + μᵖ + ξ)*(p + δ)))
expr = expr.subs(p , 0)
# nᶜ : Checks out!
# %%
expr = sym.factor(sym.cancel(top_Nᶜ/(p*(p-ω)*(p + μᵖ + ξ)*(p + δ)) - nᶜ/p),p)
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(0)
expr = sym.cancel(top*(p-ω)/bot)
expr = sym.factor(expr.subs(p , ω),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ])
# ρ₀ᶜ : Checks out!
# %%
expr = sym.cancel(top/bot - ρ₀ᶜ/(p-ω))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(ω)
expr = sym.cancel(top*(p + μᵖ + ξ)/bot)
expr = sym.factor(expr.subs(p , -(μᵖ + ξ)),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ, ρ₀ᶜ])
# ρ₁ᶜ : Checks out !
# %%
expr = sym.cancel(top/bot - ρ₁ᶜ/(p + μᵖ + ξ))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(-(μᵖ + ξ))
expr = sym.cancel(top*(p + δ)/bot)
expr = sym.factor(expr.subs(p , -δ),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ, ρ₀ᶜ, ρ₁ᶜ])
# %%
# ======================= N⁰ ========================= #
bot_N⁰ = sym.Poly(bot_N⁰, p)
roots = sym.roots(bot_N⁰)
# %%
expr = sym.cancel(top_N⁰*p/(p*(p-ω)*(p + μᵖ + ξ)*(p + δ)))
expr = expr.subs(p , 0)
# n⁰ : Checks out!
# %%
expr = sym.factor(sym.cancel(top_N⁰/(p*(p-ω)*(p + μᵖ + ξ)*(p + δ)) - n⁰/p),p)
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(0)
expr = sym.cancel(top*(p-ω)/bot)
expr = sym.factor(expr.subs(p , ω),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ])
# ρ₀⁰ : Checks out!
# %%
expr = sym.cancel(top/bot - ρ₀⁰/(p-ω))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(ω)
expr = sym.cancel(top*(p + μᵖ + ξ)/bot)
expr = sym.factor(expr.subs(p , -(μᵖ + ξ)),[c₀ᶜ, c₀⁰, n₀ᶠ, n⁰, ρ₀⁰])
# ρ₁⁰ : Checks out !
# %%
expr = sym.cancel(top/bot - ρ₁⁰/(p + μᵖ + ξ))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(-(μᵖ + ξ))
expr = sym.cancel(top*(p+δ)/bot)
expr = sym.factor(expr.subs(p , -δ),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ, ρ₀ᶜ, ρ₁ᶜ])
# ρ₂⁰ : Checks out !
# %%
# ====================================================== #
# ================= Classic to classic ================= #
# ====================================================== #
system = [(p+ξ)*Nᶜ + (ω+ξ)*N⁰ - c₀ᶜ, μᵖ*Nᶜ + (p + μᵖ - ω)*N⁰ - ( μᵖ/p + c₀⁰ )]
sol = sym.linsolve(system, [Nᶜ, N⁰])
# %%
sol = sym.factor(sol, p)
# %%
Nᶜ = sol.args[1][1]
N⁰ = sol.args[1][2]
# %%
top_Nᶜ, bot_Nᶜ = sym.fraction(Nᶜ)
top_Nᶜ =  sym.Poly(top_Nᶜ,p)
for c in top_Nᶜ.coeffs()
    println(sym.factor(c, [c₀ᶜ, c₀⁰, n₀ᶠ]))
end
# %%
top_N⁰, bot_N⁰ = sym.fraction(N⁰)
top_N⁰ =  sym.Poly(top_N⁰,p)
for c in top_N⁰.coeffs()
    println(sym.factor(c, [c₀ᶜ, c₀⁰, n₀ᶠ]))
end
# Laplace transforms seem correct !
# %%
# ======================= Nᶜ ========================= #
bot_Nᶜ = sym.Poly(bot_Nᶜ, p)
roots = sym.roots(bot_Nᶜ)
# %%
expr = sym.cancel(top_Nᶜ*p/(p*(p-ω)*(p + μᵖ + ξ)))
expr = expr.subs(p , 0)
# nᶜ : Checks out!
# %%
expr = sym.factor(sym.cancel(top_Nᶜ/(p*(p-ω)*(p + μᵖ + ξ)) - nᶜ/p),p)
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(0)
expr = sym.cancel(top*(p-ω)/bot)
expr = sym.factor(expr.subs(p , ω),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ])
# ρ₀ᶜ : Checks out!
# %%
expr = sym.cancel(top/bot - ρ₀ᶜ/(p-ω))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(ω)
expr = sym.cancel(top*(p + μᵖ + ξ)/bot)
expr = sym.factor(expr.subs(p , -(μᵖ + ξ)),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ, ρ₀ᶜ])
# ρ₁ᶜ : Checks out !
# %%
# ======================= N⁰ ========================= #
bot_N⁰ = sym.Poly(bot_N⁰, p)
roots = sym.roots(bot_N⁰)
# %%
expr = sym.cancel(top_N⁰*p/(p*(p-ω)*(p + μᵖ + ξ)))
expr = expr.subs(p , 0)
# n⁰ : Checks out!
# %%
expr = sym.factor(sym.cancel(top_N⁰/(p*(p-ω)*(p + μᵖ + ξ)) - n⁰/p),p)
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(0)
expr = sym.cancel(top*(p-ω)/bot)
expr = sym.factor(expr.subs(p , ω),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ])
# ρ₀⁰ : Checks out!
# %%
expr = sym.cancel(top/bot - ρ₀⁰/(p-ω))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(ω)
expr = sym.cancel(top*(p + μᵖ + ξ)/bot)
expr = sym.factor(expr.subs(p , -(μᵖ + ξ)),[c₀ᶜ, c₀⁰, n₀ᶠ, n⁰, ρ₀⁰])
# ρ₁⁰ : Checks out !
# %%
expr = sym.cancel(top/bot - ρ₁⁰/(p + μᵖ + ξ))
top,bot = sym.fraction(expr)
top = sym.Poly(top, p)
top = top - top(-(μᵖ + ξ))
expr = sym.cancel(top*(p+δ)/bot)
expr = sym.factor(expr.subs(p , -δ),[c₀ᶜ, c₀⁰, n₀ᶠ, nᶜ, ρ₀ᶜ, ρ₁ᶜ])
# ρ₂⁰ : Checks out !
# %%
# ====================================================== #
# ======================= Wages ======================== #
# ====================================================== #
system = [(p-ω)*Wᶜ - W₀ᶜ + (ω+ξ)*Eᶜ*(nᶜ/p + ρ₀ᶜ/(p-ω) + ρ₁ᶜ/(p-α₁) + ρ₂ᶜ/(p-α₂) + n⁰/p + ρ₀⁰/(p-ω) + ρ₁⁰/(p-α₁) + ρ₂⁰/(p-α₂)) + α₀*E₀ᶜ*n₀ᶜ, (p-ω)*W⁰ - W₀⁰ - μᵖ*E⁰*(1/p - (nᶜ/p + ρ₀ᶜ/(p-ω) + ρ₁ᶜ/(p-α₁) + ρ₂ᶜ/(p-α₂)) - (n⁰/p + ρ₀⁰/(p-ω) + ρ₁⁰/(p-α₁) + ρ₂⁰/(p-α₂)) - (nᶠ/p + ρ₁ᶠ/(p-α₁) + ρ₂ᶠ/(p-α₂))) + β₀*E₀⁰*n₀⁰]
sol = sym.linsolve(system, [Wᶜ, W⁰])
# %%
Wᶜ = sol.args[1][1]
W⁰ = sol.args[1][2]
# %%
top_Wᶜ, bot_Wᶜ = sym.fraction(Wᶜ)
top_Wᶜ =  sym.Poly(top_Wᶜ,p)
# %%
top_W⁰, bot_W⁰ = sym.fraction(W⁰)
top_W⁰ =  sym.Poly(top_W⁰,p)
# Laplace transforms seem correct !
# %%
# ======================= Wᶜ ========================= #
bot_Wᶜ = sym.Poly(bot_Wᶜ, p)
roots = sym.roots(bot_Wᶜ)
# %%
expr = sym.cancel(top_Wᶜ/((p-ω)^2*(p-α₁)*(p-α₂)))
W∞ᶜ = sym.cancel(expr.subs(p , 0))
# %%
expr = sym.cancel(top_Wᶜ/(p*(p-α₁)*(p-α₂)))
ν₀ᶜ = sym.cancel(expr.subs(p , ω))
# %%
expr = sym.cancel(top_Wᶜ/(p*(p-ω)*(p-α₁)*(p-α₂)) - ν₀ᶜ/(p-ω))
expr = expr.subs(p, ω)
sym.apart(expr, ω)
# %%
expr = sym.cancel(top_Wᶜ/(p*(p-ω)^2*(p-α₂)))
ν₁ᶜ = sym.cancel(expr.subs(p , α₁))
# ν₁ᶜ : Checks out !
# %%
expr = sym.cancel(top_Wᶜ/(p*(p-ω)^2*(p-α₁)))
ν₂ᶜ = sym.cancel(expr.subs(p , α₂))
# ν₂ᶜ : Checks out !
# %%
# ======================= W⁰ ========================= #
bot_W⁰ = sym.Poly(bot_W⁰, p)
roots = sym.roots(bot_W⁰)
# %%
expr = sym.cancel(top_W⁰/((p-ω)^2*(p-α₁)*(p-α₂)))
sym.cancel(expr.subs(p , 0))
# %%
expr = sym.cancel(top_W⁰/(p*(p-α₁)*(p-α₂)))
ν₀⁰ = sym.cancel(expr.subs(p , ω))
# %%
expr = sym.cancel(top_W⁰/(p*(p-ω)*(p-α₁)*(p-α₂)) - ν₀⁰/(p-ω))
expr = expr.subs(p, ω)
sym.apart(expr, ω)
# %%
expr = sym.cancel(top_W⁰/(p*(p-ω)^2*(p-α₂)))
ν₁⁰ = sym.cancel(expr.subs(p , α₁))
# ν₁ᶜ : Checks out !
# %%
expr = sym.cancel(top_W⁰/(p*(p-ω)^2*(p-α₁)))
ν₂⁰ = sym.cancel(expr.subs(p , α₂))
#