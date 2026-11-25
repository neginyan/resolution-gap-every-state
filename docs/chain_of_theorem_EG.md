# The chain of Theorem EG, item by item

This is the working proof of Theorem EG (Theorem 1 of the paper) as one chain: every lemma in the order of use, with unified notation, its units, its hypotheses, and the joints between items.
It is the English version of the chain document of the lab, with the seven corrections of the external review applied (marked **[review fix n]**) and with the range of the paper (J₀ = 1.9×10⁻⁸, σ ≤ 1.5×10⁻¹⁴).
It adds no new claim. The paper gives the same chain in its own numbering; the map is at the end.

- "Uses" lists the items used. "Units" says whether the formulas are in raw coordinates, dimensionless with σ written explicitly, or for the state X/σ with unit noise (σ = 1).
- Every numerical constant was recomputed with outward-rounded interval arithmetic (`src/eh_constants.py`, `src/ei_constants.py`; item 30).
- Labels: [proof], [computer-assisted proof], [derivation], [numerical], [conjecture].

---

## 0. Setting and notation

### 0.1 State and noise
- **State:** ρ is a probability measure on R²; X ~ ρ. **Assumption (M1): E|X| < ∞.** No density is needed.
- **Observation:** Y = X + σZ with Z ~ N(0, I₂) independent of X, σ > 0. **Z := (Y − X)/σ** (this sign throughout).
- **Blurred density:** ϱ := ρ * N(0, σ²I). E is the expectation over Y ~ ϱ.
- **Coordinates:** (v, w) are the unstable and stable directions at the upper equilibrium of the pendulum; **z := (v + w)/√2, p := (v − w)/√2** (the (z, p) of paper E).
- **Assumption (H2):** Var(X_z) ≤ π² (the observer looks at the neighbourhood of the equilibrium, one period in z).
- ℓ := ln(1/σ), **α := 1/(2(π² + σ²))** (the constant of the Cramér–Rao floor).

### 0.2 Scores and Tweedie
- Raw: s := ∇ln ϱ, H := −∇²ln ϱ, m := E[X | Y] = Y + σ²s, J_vv := E[s_v²], J_ww := E[s_w²].
- Dimensionless: **s̃ := σs = −E[Z | Y], H̃ := σ²H = I − Cov(Z | Y)**.
- **A := σ²J_ww = E[s̃_w²], B := σ²J_vv = E[s̃_v²]**, J := A + B (paper E's a is B and its d is A + B).
- φ := 1 − cos m_z (at the posterior mean), **φ₀ := 1 − cos Y_z** (at the observed point).
- θ := E[s̃_v²φ]/B, **M := A + θB/2**, **O := {M < 112σ²}** (an open set of states).
- A_w := E[s̃_w²(1 − φ/4)], N_φ := ½E[s̃_v²φ] = θB/2; A_w⁰, N_φ⁰ with φ₀ in place of φ.
- Branches: (a) := 2A_w/(σ(A + B)), (a′) := N_φ/(σ(A + B)).
- W_zz := E[Z_z² − 1 | Y], U := E[s̃_p² sin²Y_z].
- C̃_w := E[s̃_v² cos m_z (1 − H̃_ww − H̃_vw)]/B, Γ := Cov(Z_w, Z_v + Z_w | Y).

### 0.3 Gap and remainders
- Pendulum: ż = p, ṗ = sin z, κ = 1. Critical rate R := E[s·E[f(X) | Y]]/(σ²E|s|²), **gap g := (1 − R)/σ**.
- f₀(x) := x − sin x; 𝒟(y) := E[f₀(X_z) | Y = y]; N := E[s_p𝒟] (a dimensionless number).
- Paper F's remainder: **B_F** (called B in papers E and F), R₃, T₂, R_L. D_F := 2√(2A_wN_φ), ε := −(B_F + R₃)/D_F, r_F := −(B_F + R₃)/(2A_w + N_φ).
- Remainders around the observed point: E₂ := (σ/2)E[s̃_p sin Y_z W_zz], R′, r⁰ := (E₂ + R′)/(2A_w⁰ + N_φ⁰).

### 0.4 The dimensionless heat flow
- Divide y by σ and write ρ̃ for the law of X/σ. 𝒩_t(u) := (2πt)^{−1/2}e^{−u²/(2t)} (principal square root for Re t > 0).
- **ϱ_t(v, w) := ∫𝒩_t(v − x_v)𝒩₁(w − x_w)dρ̃(x)** (v-noise of variance t, w-noise 1), q_t := ∂_wϱ_t; ϱ₁ is ϱ rescaled by σ.
- **F(t) := ∫q_t²/ϱ_t**, **F(1) = A**, **L := ln(1/A)**.
- D := E[(∂_v s̃_w)²], Q₂ := E[(E[Z_vZ_w | Y])²], SV := E[s̃_w² Var(Z_v | Y)], S2 := E[s̃_v²s̃_w²], **K := E[(∂_v s̃_w + ½s̃_vs̃_w)²]**, **k_eff² := K/A**.
- **ψ(L) := 289L^{5/2}/(ln 2L)^{3/2} + 4L + 3.57**, **L̄ := 2ℓ + ln(1/α)**, **k̄ := √ψ(L̄)**.
- Numbers of item 18: δ := 1/(4L), T := 1 + δ, η := ln(2L)/(4L), **τ := ln(2L)/(8L²)**, r₁ := δ(1 − η), r₂ := δ + τ, r₃ := T/2, **ϑ := ln(r₂/r₁)/ln(r₃/r₁)**, **ω := (τ + δ)/(2(1 − τ))**, c_T := √(T/(1 − τ)), K_b² := 4L, W_T(y) := E[e^{−ωZ_v²} | Y_T = y].
- Isotropic flow: ϱ̃^iso_t := ρ̃ * N(0, tI₂), I_e(t) := E[(∂_e ln ϱ̃^iso_t)²], L_e := ln(1/I_e(1)).
- h_O(x) := 562·x·(ln(1/x))³/(ln(2 ln(1/x)))², **J₀ := 1.9×10⁻⁸**, Ξ := ‖W_zz‖₂.
- Region constants: c_r := √0.1, K_r := 400, β′.

### 0.5 Names that collided in the working record (old name → name here)

| old name | meaning 1 | meaning 2 (and more) | here |
|---|---|---|---|
| B | part of the remainder of papers E, F | the information σ²J_vv | B_F and B |
| φ | 1 − cos m_z | 1 − cos Y_z (EG-D), the Gaussian kernel φ_z (EG-1) | φ, φ₀, 𝒩_z |
| θ | E[s̃_v²φ]/B | exponent of the three circles | θ and ϑ |
| α | Cramér–Rao floor 1/(2(π² + σ²)) | exponent of the Gaussian ratio (EG-1) | α and ω |
| s | score | the return fraction (EG-1), the split fraction (Q-sm) | s, τ, u |
| δ | 1/(4L) (EG-1) | Taylor's δ = x_z − m_z of paper F, δ_w, δ_c | δ, ζ |
| K | ED's K | region constant (EE 3200, EF/EG 400), K² = 4L (EG-1), a tail cut | K, K_r, K_b |
| c | √0.1 | √(T/(1 − s)) (EG-1), c₁–c₄ | c_r, c_T, c₁–c₄ |
| g | gap | x − sin x (EG) | g and f₀ |
| D | ED's E[(∂_v s̃_w)²] | paper E's D(y), paper F's D, paper E's d | D, 𝒟(y), D_F, A + B |
| p | density | coordinate (v − w)/√2 | ϱ and p |
| h | h(x) of EG-O | h(y′) = E[Z_w | Y′] of Q-sm, h_t of the 1-D core | h_O, h_{1−u} |
| G | Cov(Z_w, Z_v + Z_w | Y) of EE-C | G(t) of EG-iso and EF-S | Γ, I_e(t) |
| L | ln(1/A) | ln(1/G) (EG-iso), window width of paper F | L, L_e |
| W | W(y) of EG-1 | W_zz, paper F's W₂ | W_T, W_zz |
| r | radii r₁, r₂, r₃ | the ratio r (EF), ED's Hölder exponent r | r₁–r₃, r_F, r⁰, r_H |
| ρ | state | Taylor remainder of paper F, ρ₃ of EG-D, ρ* | ρ, ρ_L, ρ_D |
| ψ | ψ(L) | ψ = √2 sin m_z (AT2′), the test vector of EG-1c | ψ(L), ψ_AT, χ |
| M | A + θB/2 | moments M_b = E|Z|^{2b} (ED, EE), M(r) of the three circles | M, M_b, 𝔐(r) |
| a | paper E's a = σ²J_vv | branch (a), a = F₁(1) of the 1-D core | B, (a), a₁ |
| sign of Z | Z = (Y − X)/σ | (X − Y)/σ in a first draft of EF §2.0 and parts of paper F's code | fixed to (Y − X)/σ (s̃ = −E[Z | Y]; even quantities unchanged) |

---

## 1. Integration by parts [proof]
- **Uses:** none. **Units:** raw (any σ).
- **Claim.** Let ρ be any probability measure.
  - (i) Let G(y) = Ψ(y)·∏_j E[q_j(Z) | Y = y] with polynomials q_j and a C¹ function Ψ such that Ψ and ∇Ψ are bounded by constants times such products, or a finite sum of such G. Then ϱG ∈ W^{1,1}(R²) and **∫∂_i(ϱG) = 0**.
  - (ii) Under (M1), for every affine ℓ(y): ℓϱ ∈ W^{1,1} and **E[s_iℓ(Y)] = −∂_iℓ**.
- **Proof.**
  1. ϱ(y)E[q(Z) | y] = ∫q((y − x)/σ)𝒩_{σ²}(y − x)dρ(x). The integrand is smooth in y and its derivatives are polynomials times a Gaussian, so differentiation under the integral gives ∂_{y_i}[ϱE[q | y]] = (ϱ/σ)E[∂_iq − Z_iq | y]. **[review fix 1]** Hence ∂_i(ϱG) is a finite sum of terms ϱΨ̃∏_jE[q̃_j | y] with |Ψ̃| bounded by products of the same kind (this is where ∇Ψ enters). *We do not claim that the class is closed under differentiation;* we only use that ∂_i(ϱG) is bounded by a finite sum of the same form and is therefore in L¹ by step 2.
  2. Z ~ N(0, I) does not depend on ρ. By Hölder and Jensen, ∫ϱ|∏_{j≤n}E[q_j | y]| ≤ ∏_j‖q_j(Z)‖_{L^n} < ∞. So ϱG ∈ L¹ and, by step 1, ∂_i(ϱG) ∈ L¹: ϱG ∈ W^{1,1}.
  3. For f ∈ W^{1,1} and a smooth cut-off χ_R(y) = χ(y/R), ∫χ_R∂_if = −∫f∂_iχ_R → 0, and the left side tends to ∫∂_if by dominated convergence.
  4. (ii): |∂_iϱ(y)| ≤ ∫(|y_i − x_i|/σ²)𝒩_{σ²}(y − x)dρ(x), so ∫|ℓ||∂_iϱ| ≤ ∫dρ(x)∫(|ℓ(x)| + |ℓ(u)|)(|u_i|/σ²)𝒩_{σ²}(u)du < ∞ by (M1); likewise ∫|ℓ|ϱ < ∞. Step 3 for ℓϱ gives ∫ϱ∂_iℓ + ∫ℓ∂_iϱ = 0. ∎
- **Joint.** Every integration by parts below is of form (i) or (ii); only (ii) uses (M1).

## 2. Tweedie's formulas [proof]
- **Uses:** [1]. **Units:** raw and dimensionless.
- **Claim.** s̃ = −E[Z | Y], H̃ = I − Cov(Z | Y), m = Y + σ²s; σ²∂_i∂_jϱ/ϱ = E[Z_iZ_j | Y] − δ_ij, in particular **W_zz = σ²∂_z²ϱ/ϱ** and, at σ = 1, **∂_v∂_wϱ/ϱ = E[Z_vZ_w | Y]**; the posterior cumulants of X of order n ≥ 3 are σ^{2n}∂ⁿln ϱ; Cov(X_p, F(X) | y) = σ²∂_pE[F(X) | y] for F bounded or of linear growth.
- **Proof.** ∂_{y_i}𝒩_{σ²}(y − x) = −(Z_i/σ)𝒩, ∂_i∂_j𝒩 = ((Z_iZ_j − δ_ij)/σ²)𝒩; divide by ϱ. The cumulants are derivatives of ln ϱ(y) = −|y|²/(2σ²) + ln∫e^{y·x/σ²}(…)dρ. For the last formula integrate x_p𝒩_{σ²}(y − x) = (y_p + σ²∂_{y_p})𝒩_{σ²}(y − x) against F and divide by ϱ. ∎
- **Joint.** s̃, H̃, W_zz are posterior moments of Z, so they are of the form of [1](i).

## 3. Paper E's gap formula and three forms of N [proof]
- **Uses:** [1], [2]. **Units:** σg and N are dimensionless numbers; s is raw inside formulas. **Assumption:** (M1).
- **Claim.** R is well defined, **A + B > 0**, and
  > **σ(A + B)·g = 2A + N**, N = E[s_p𝒟] = −E[∂_p𝒟] = −σ⁻²E[Cov(X_p, f₀(X_z) | Y)] = **E[s_p(Y)(f₀(X_z) − f₀(Y_z))]**,
  > and **E|Cov(X_p, f₀(X_z) | Y)| ≤ 2σ²**.
- **Proof.**
  0. **[review fix 2]** A + B = σ²E|s|² > 0: s ≡ 0 would make ϱ constant, which is impossible for an integrable density. So the denominator of R is positive.
  1. E[f(X) | y] = (m_p, m_z − 𝒟) in the (z, p) components, so s·E[f | y] = s_zm_p + s_pm_z − s_p𝒟.
  2. Integrability: m = y + σ²s; E|s_zy_p| < ∞ by the estimate in the proof of [1](ii), E|s_zs_p| ≤ E|s|² < ∞; 𝒟 − f₀(y_z) = E[f₀(X_z) − f₀(Y_z) | y] and |f₀′| ≤ 2, so |𝒟 − f₀(y_z)| ≤ 2σE[|Z_z| | y].
  3. By [1](ii), E[s_zy_p] = E[s_py_z] = 0. So E[s_zm_p + s_pm_z] = 2σ²E[s_zs_p] = σ²(J_vv − J_ww). With κ = 1 and σ²E|s|² = A + B this gives σg = (2A + N)/(A + B).
  4. **[review fix 2]** Split 𝒟 = Ψ + f₀(y_z) with Ψ(y) := E[f₀(X_z) − f₀(Y_z) | Y = y]. Then **E[s_p f₀(Y_z)] = E[s_pY_z] − E[s_p sin Y_z] = 0 − 0**, by [1](ii) (y_z does not depend on p) and by [1](i) (sin y_z is bounded, smooth and independent of p). Since |Ψ| ≤ 2σE[|Z_z| | y], ϱΨ ∈ W^{1,1} by the argument of [1], so N = E[s_pΨ] = −E[∂_pΨ] = −E[∂_p𝒟], and ∂_p𝒟 = σ⁻²Cov(X_p, f₀(X_z) | Y) by [2]. This gives the covariance form and the new form.
  5. **[review fix 3]** Integrability of the covariance. f₀ is 2-Lipschitz, so Var(f₀(X_z) | Y) ≤ E[(f₀(X_z) − f₀(Y_z))² | Y] ≤ 4σ²E[Z_z² | Y]; Var(X_p | Y) = σ²Var(Z_p | Y) ≤ σ²E[Z_p² | Y]. The conditional Cauchy–Schwarz inequality and then Cauchy–Schwarz in Y give E|Cov(X_p, f₀(X_z) | Y)| ≤ 2σ²(E Z_z² E Z_p²)^{1/2} = 2σ². So only the first moment of X is ever needed. ∎

## 4. Paper E's branch identity and two simple lower bounds [proof]
- **Uses:** [1], [2], [3]. **Units:** dimensionless (σ explicit). **Assumption:** (M1).
- **Claim.** (a) **σ(A + B)·g = 2A_w + N_φ + (B_F + R₃)** with R₃ := N − σ²E[H_pzφ] and B_F := ½σ²E[s_v∂_vφ − s_w∂_wφ] (raw H, s). (b) A_w ≥ A/2, so (a) ≥ A/(σ(A + B)) and M ≤ 2A_w + N_φ. (c) D_F = 2√(2A_wN_φ) = 2σ(A + B)√((a)(a′)), g = (a) + (a′) − 2ε√((a)(a′)) = (1 − r_F)((a) + (a′)), |r_F| ≤ |ε|.
- **Proof.** (a) H_pz = ½(H_vv − H_ww); by [1](i), E[H_vvφ] = E[s_v²φ] + E[s_v∂_vφ] (φ is bounded and of the admissible form), and likewise for w; add to 2A of [3]: 2A − ½E[s̃_w²φ] = 2A_w. (b) 0 ≤ φ ≤ 2, so 1 − φ/4 ∈ [½, 1]. (c) From the definitions; D_F ≤ 2A_w + N_φ (AM–GM) gives |r_F| ≤ |ε|. ∎
- **Joint.** Paper F's Lemma BB1 (g ≥ (1 − |ε|)θ/(2σ)) is not used; region 2 of [23] replaces it.

## 5. Paper F's Lemma G2 (the three parts of the remainder) [proof]
- **Uses:** [1], [2], [4]. **Units:** raw H, s. **Assumption:** none (the "sufficient decay" of papers E, F is removed by [1]).
- **Claim.** B_F = (σ⁴/2)E[cos m_z H_pz] − (σ⁴/2)E[sin m_z(s_vH_vz − s_wH_wz)], R₃ = T₂ + R_L with T₂ = −(σ⁴/2)E[sin m_z ∂_pzz ln ϱ] and **|R_L| ≤ (√960/6)σ²**.
- **Proof.** 1. ∂_vφ = sin m_z ∂_vm_z, ∂_vm_z = 1/√2 − σ²H_vz, likewise for w; the first part is (σ²/2)E[s_p sin m_z] = −(σ²/2)E[∂_p sin m_z] = (σ⁴/2)E[cos m_z H_pz] by [1](i). 2. With ζ := X_z − m_z, f₀(X_z) = f₀(m_z) + φζ + ½ sin m_z ζ² + ρ_L, |ρ_L| ≤ |ζ|³/6; by [2], Cov(X_p, X_z | y) = −σ⁴H_pz and Cov(X_p, ζ² | y) = σ⁶∂_pzz ln ϱ; the first order is σ²E[H_pzφ], the second is T₂. 3. R_L = −σ⁻²E Cov(X_p, ρ_L | Y); Cauchy–Schwarz twice gives |R_L| ≤ (σ⁻²/6)√(E Var(X_p | Y))√(Eζ⁶), E Var(X_p | Y) ≤ σ², Eζ⁶ ≤ 2⁶·15σ⁶ = 960σ⁶. ∎

## 6. Lemma EF-1 (the numerator) [proof]
- **Uses:** [1], [2], [5]. **Units:** dimensionless A, B, θ, explicit σ.
- **Claim.** **|B_F + R₃| ≤ c₁σ² + c₂σ√(θB) + c₃σ√A** with c₁ = √3/2 + √3 + √960/6 ≤ 7.77, c₂ = (2√2 + √(3/2))/√2 ≤ 2.87, c₃ = (2√2 + √(3/2))/2 ≤ 2.03.
- **Proof.** Integrate T₂ once in p ([1](i)): T₂ = −(σ⁴/2)E[sin m_z s_pH_zz] + (σ⁶/2)E[cos m_z H_pzH_zz]. In dimensionless form, B_F + R₃ = (σ²/2)E[cos m_z H̃_pz] + (σ²/2)E[cos m_z H̃_pzH̃_zz] + R_L − (σ/2)E[sin m_z s̃_vX₊] + (σ/2)E[sin m_z s̃_wX₋], X± := √2H̃_zz ± H̃_pz/√2. Tweedie: ‖H̃_zz‖₂ ≤ 2, ‖H̃_pz‖₂ ≤ √3, ‖X±‖₂ ≤ 2√2 + √(3/2). With sin²m_z ≤ 2φ, Cauchy–Schwarz term by term. ∎
- **Numerical.** On 692 recorded engine states (200 inside O) the smallest margin is +0.678 (joint J2).

## 7. Lemma EF-2 (the threshold) [proof]
- **Uses:** [4], [6]. **Claim.** If M ≥ 112σ², then **|B_F + R₃| ≤ ½M ≤ ½(2A_w + N_φ)** and **g ≥ ½((a) + (a′))**.
- **Proof.** 2.03√A + 2.87√(θB) ≤ c₄√M with c₄ = √(2.03² + 2·2.87²) ≤ 4.5382. With x = √M/σ, |B_F + R₃| ≤ σ²(7.77 + c₄x); ½x² − c₄x − 7.77 increases for x ≥ c₄ and is 0.2028… > 0 at x = √112 (interval arithmetic). Insert into [4]. ∎

## 8. Lemma EE-3 (the Cramér–Rao floor) [proof]
- **Uses:** [1]. **Assumption:** (H2).
- **Claim.** **A ≥ ασ²**, hence **L = ln(1/A) ≤ L̄**. Likewise σ²J_zz ≥ 2ασ².
- **Proof.** [1](ii): E[s_w(Y_z − EY_z)] = −1/√2; Cauchy–Schwarz gives ½ ≤ J_ww Var(Y_z), and Var(Y_z) = Var(X_z) + σ² ≤ π² + σ². For J_zz, E[s_z(Y_z − EY_z)] = −1. ∎

## 9. ED's Lemma 0 (dimensionless form of paper F's k_eff) [proof]
- **Claim.** Paper F's k_eff² = σ²∫(∂_v∂_w√ϱ)²/∫(∂_w√ϱ)² equals **K/A** after rescaling by σ.
- **Proof.** ∂_w√ϱ = ½√ϱ s_w, ∂_v∂_w√ϱ = ½√ϱ(∂_vs_w + ½s_vs_w); y = σy′ removes σ. ∎
- **Numerical.** On 692 states paper F's code and the sums of ED agree to a relative 2×10⁻⁷ (J6).

## 10. ED's Lemma 2 (the half lemma) [proof]
- **Claim.** **S2 + SV = E[s̃_w²Z_v²] ≤ 8AL + 7.136A**; in particular SV ≤ 8AL + 7.136A.
- **Proof.** Cut at Z_v² = 8L: the inside is ≤ 8LA; outside, Hölder and E s̃_w⁴ ≤ 3 give ≤ √3(E[Z_v⁴e^{(Z_v² − 8L)/4}])^{1/2} = √3(3·2^{5/2})^{1/2}A = 3·2^{5/4}A ≤ 7.136A. ∎

## 11. ED's Lemma 3 (the identity) [proof]
- **Claim.** **K = D + ½A − ½SV − ¼S2.**
- **Proof.** K = D + E[s̃_vs̃_w∂_vs̃_w] + ¼S2 and E[s̃_vs̃_w∂_vs̃_w] = −½E[s̃_w²∂_v²ϱ/ϱ] = −½(S2 + SV) + ½A by [1](i), [2]. ∎

## 12. Lemma Q-id [proof]
- **Claim.** **D = Q₂ + SV − A**, hence **K = Q₂ + ½SV − ½A − ¼S2 ≤ Q₂ + ½SV**.
- **Proof.** ∂_v∂_wϱ/ϱ = ∂_vs̃_w + s̃_vs̃_w = E[Z_vZ_w | Y] ([2]); Q₂ = D + 2E[s̃_vs̃_w∂_vs̃_w] + S2. ∎

## 13. Lemma Q-sm (smoothing) [proof]
- **Claim.** For every u ∈ (0, 1), **Q₂ ≤ F(1 − u)/u**.
- **Proof.** Z_v = √(1 − u)Z₁ + √uZ₂, Y′ := (X_v + √(1 − u)Z₁, Y_w) has density ϱ_{1−u}, Y = Y′ + (√uZ₂, 0). E[Z₂Z_w | Y] = √uE[Z_vZ_w | Y] = E[Z₂h_{1−u}(Y′) | Y] with h_{1−u} = −∂_w ln ϱ_{1−u}; Jensen and Z₂ ⫫ Y′. ∎

## 14. Lemma EG-K (kernel integrals) [proof]
- **Claim.** For T > 0, |z − T| = r < T: ∫|𝒩_z|²/𝒩_T = (1 − (r/T)²)^{−1/2} and ∫u²|𝒩_z|²/𝒩_T = (1 − (r/T)²)^{−1/2}·T|z|²/(T² − r²).
- **Proof.** |𝒩_z(u)|² = (2π|z|)⁻¹e^{−u²Re z/|z|²}; for z = T + re^{iθ₀}, 2T Re z − |z|² = T² − r² for every angle; Gaussian integrals of variance T|z|²/(T² − r²). ∎

## 15. Lemma EG-1a (outer circle) [proof]
- **Claim.** For |z − T| = r < T, **‖q_z‖²_{H_T} ≤ (1 − (r/T)²)^{−1/2}**, H_T := L²(dy/ϱ_T).
- **Proof.** q_z(y) = −∫𝒩_z(v − x_v)(w − x_w)𝒩₁(w − x_w)dρ̃; Cauchy–Schwarz in ρ̃ with the weight √(𝒩_T𝒩₁) gives |q_z|² ≤ ϱ_T∫(|𝒩_z|²/𝒩_T)(w − x_w)²𝒩₁dρ̃; divide by ϱ_T, integrate (Tonelli), ∫w′²𝒩₁ = 1 and [14]. ∎

## 16. Lemma EG-1b (inner circle) [proof]
- **Claim.** For T = 1 + δ and |z − T| = r < δ, **‖q_z‖²_{H_T} ≤ A(1 − (r/δ)²)^{−1/2}**.
- **Proof.** Re(z − 1) > 0; the semigroup law 𝒩_{a+b} = 𝒩_a * 𝒩_b holds for Re a, Re b > 0; q_z = 𝒩_{z−1} *_v q₁, ϱ_T = 𝒩_δ *_v ϱ₁; Cauchy–Schwarz in v′ with the weight 𝒩_δ(v − v′)ϱ₁(v′, w); [14] for (δ, z − 1); ∫q₁²/ϱ₁ = A. ∎

## 17. Lemma EG-1c (holomorphy and the three circles) [proof]
- **Uses:** [14], [15].
- **Claim.** 𝔇 := {|z − T| < T}. (i) For every χ ∈ H_T, f_χ(z) := ∫χ̄q_z/ϱ_T is holomorphic on 𝔇. (ii) For 0 < r₁ < r₂ < r₃ < T and |z₂ − T| = r₂, ‖q_{z₂}‖ ≤ (sup_{r₁}‖q‖)^{1−ϑ}(sup_{r₃}‖q‖)^ϑ with ϑ = ln(r₂/r₁)/ln(r₃/r₁).
- **Proof.** **[review fix 4]** (i) Let 𝒦 ⊂ 𝔇 be compact. On 𝒦, r = |z − T| ≤ r_𝒦 < T and |z| lies between two positive constants, so by the formula in the proof of [14], **|𝒩_z(u)|²/𝒩_T(u) ≤ C_𝒦e^{−c_𝒦u²}** with c_𝒦 = min_𝒦 (T² − r²)/(2T|z|²) > 0. Let **G_𝒦(y) := C_𝒦∫e^{−c_𝒦(v−x_v)²}(w − x_w)²𝒩₁(w − x_w)dρ̃(x)**; it does not depend on z, and ∫G_𝒦 dy = C_𝒦√(π/c_𝒦) < ∞. The Cauchy–Schwarz step of [15] gives |q_z|² ≤ ϱ_TG_𝒦 for every z ∈ 𝒦, hence **|χ̄q_z/ϱ_T| ≤ (|χ|/√ϱ_T)√G_𝒦 ∈ L¹** and **sup_{z∈𝒦}‖q_z‖²_{H_T} ≤ ‖G_𝒦‖_{L¹}** (the local bound of the H_T norm). For each y, z ↦ q_z(y) is holomorphic on Re z > 0 (the integrand is holomorphic in z and bounded on compact sets uniformly in x), so ∮_{∂Δ}q_z(y)dz = 0 for every triangle Δ ⊂ 𝔇 by **Cauchy's theorem**. By the dominating function and Fubini, ∮_{∂Δ}f_χ dz = ∫(χ̄/ϱ_T)∮_{∂Δ}q_z dz dy = 0, and f_χ is continuous by dominated convergence; Morera's theorem gives (i). (ii) Hadamard's three-circle theorem for f_χ with ‖χ‖ = 1, |f_χ(z)| ≤ ‖q_z‖, and the supremum over χ. ∎

## 18. Theorem EG-1 (backward stability and the bound on k_eff) [proof]
- **Uses:** [9], [10], [12], [13], [15], [16], [17]. **Units:** σ = 1. **Assumption:** none (every probability measure).
- **Claim.** If A ≤ e⁻⁴ (L ≥ 4): (i) ‖q_{1−τ}‖²_{H_T} ≤ (4e/√3)√(L/ln 2L)A; (ii) **F(1 − τ) ≤ 36.1√(L/ln 2L)A**; (iii) **Q₂ ≤ 289AL^{5/2}/(ln 2L)^{3/2}**; (iv) **k_eff² ≤ ψ(L)**; (v) if A > e⁻⁴, k_eff² ≤ 74.2 ≤ ψ(4). So **every state has k_eff² ≤ ψ(max(L, 4))**.
- **Proof.**
  1. *Numbers.* For L ≥ 4: η ≤ 0.12997, τ ≤ 0.016246, δ ≤ 1/16 (all decrease in L); r₁ < r₂ < r₃ < T; 1 − τ lies on |z − T| = r₂; Re z ≥ T/2 on |z − T| ≤ r₃.
  2. *Circles.* Outer ([15], r = T/2): ≤ 2/√3. Inner ([16], r = r₁): ≤ A/√(2η − η²) ≤ A/√η.
  3. *The exponent.* r₃/r₁ > 2L, ln(r₂/r₁) ≤ τ/δ + 2η = ln(2L)/L (for η ≤ 0.13, −ln(1 − η) ≤ 2η), so ϑ < 1/L.
  4. *Three circles* ([17], squared): ‖q_{1−τ}‖² ≤ (A/√η)^{1−ϑ}(2/√3)^ϑ ≤ e·A·2√(L/ln 2L)·(2/√3). This is (i).
  5. *Back to the weight 1/ϱ_{1−τ}.* 𝒩_{1−τ}/𝒩_T = c_Te^{−ωζ²} (u = √Tζ), so exactly ϱ_{1−τ} = c_Tϱ_TW_T and q_{1−τ} = −c_Tϱ_TE[Z_we^{−ωZ_v²} | Y_T].
  6. *ωK_b² < 1:* ω < 1/(4L) ⟺ ln(2L)(1 + 1/(2L)) < 2L, true for L ≥ 4.
  7. *Valleys.* E_bad := {W_T < 1/(2e)}; there P(|Z_v| > K_b | y) > ½, so P_T(E_bad) ≤ 2A². The contribution of E_bad to F(1 − τ) is ≤ √(6T/(1 − τ))A ≤ 2.55A; outside, ϱ_{1−τ} ≥ c_Tϱ_T/(2e) and the contribution is ≤ 34.13√(L/ln 2L)A. Since √(L/ln 2L) increases and is 1.38693… at L = 4, the sum is ≤ 36.1√(L/ln 2L)A. This is (ii).
  8. (iii) is [13] with u = τ, 1/τ = 8L²/ln 2L, 36.1·8 ≤ 289. (iv) [12], [10] and [9]: k_eff² ≤ Q₂/A + 4L + 3.568. (v) Q₂ ≤ 1, ½SV/A < 19.57 for L < 4, e⁴ + 19.57 < 74.17 ≤ ψ(4) = 3103.7; ψ increases on [4, ∞). ∎
- **On the choice of δ. [review fix 7]** δ = 1/(4L) is **one sufficient choice**, not a forced one. The return step (step 7) needs only ωK_b² < 1; for instance δ = 0.4/L gives ωK_b² = (0.8 + ln(2L)/(4L))/(1 − τ) = 0.9453 < 1 at L = 4 (interval arithmetic, `ei_constants.py` key d04), with other constants. We have not optimized the choice. What the argument does need is that δ ln(1/A) stays bounded: a v-noise of order one would fill the valleys of ϱ₁ and the weight 1/ϱ_T could not be turned back into 1/ϱ_{1−τ}. That is why the theorem gives τ ~ ln L/L² and the power L^{5/2}. (The working record of stage EG called δ = 1/(4L) "the best within this toolset"; that statement is withdrawn.)
- **Numerical.** Every step, the points of the circles included, on 1441 two-dimensional states (889 with L ≥ 4), the conclusions (ii)–(iv) on 2449 two-dimensional states, and 2921 one-variable states: no step fails; smallest margins ln(right/left) 0.702 (inner circle) and 1.498 (three circles); ϑL ≤ 0.75; F(1 − τ)/A ≤ 1.06.

## 19. The Heisenberg step of AT2′ (paper F's Theorem AT2 with ψ_AT = √2 sin m_z) [proof]
- **Uses:** [1], [2], [9]. **Units:** raw → dimensionless (the joint where a power of σ was once lost).
- **Claim.** **σ²C̃_w²B² ≤ 16θBk_eff²A.**
- **Proof.** f := s_v√ϱ = 2∂_v√ϱ; ∫f²∂_wψ_AT = C̃_wJ_vv; by [1](i) ∫f²∂_wψ_AT = −2∫ψ_ATf∂_wf; Cauchy–Schwarz; ψ_AT² ≤ 4φ; ∫(∂_wf)² = k_eff²J_ww/σ². In dimensionless form (C̃_wB/σ²)² ≤ 16(θB/σ²)(k_eff²A/σ⁴). ∎
- **Joint.** Without the σ⁻² of the last step the inequality is false; the check J18 finds it broken by a factor 9×10⁷¹.

## 20. EE's identity and Lemma EE-C (improved form) [proof]
- **Claim.** E[s̃_v²H̃_ww] = A − SV ≤ A, and **C̃_w ≥ 1 − (A + 1.291√A)/B − 8.4258√(θ/B)**.
- **Proof.** Γ = 1 − H̃_ww − H̃_vw, cos m_z = 1 − φ: C̃_wB = E[s̃_v²Γ] − E[s̃_v²φΓ]. E[s̃_v²H̃_vw] = E[s̃_v³s̃_w]/3, |…| ≤ √(15A)/3, so E[s̃_v²Γ] ≥ B − A − 1.291√A. Γ² ≤ E[Z_w² | Y]E[(Z_v + Z_w)² | Y] gives EΓ⁴ ≤ (105·1680)^{1/2} = 420, E[s̃_v²Γ²] ≤ √(3·420), and |E[s̃_v²φΓ]| ≤ (2√1260)^{1/2}√(θB) ≤ 8.4258√(θB). ∎

## 21. Lemma EF-3 (the denominator) [proof]
- **Claim.** D_F ≥ √(2AθB) and **D_F ≥ σ|C̃_w|B/(2√2k_eff)**.

## 22. Lemma EF-O (the v-information inside O) [proof]
- **Uses:** [20], [21]. **Assumption:** σ ≤ 10⁻³.
- **Claim.** If M < 112σ², then **B ≤ σ max(280, 1268k_eff)**.
- **Proof.** A < 112σ², θB < 224σ²; [20] gives C̃_w > 1 − 140σ/B. If B ≥ 280σ, C̃_w ≥ ½ and [19] gives ½σB ≤ 4k_eff√(224σ²·112σ²), i.e. B ≤ 1267.2k_effσ. ∎

## 23. Theorem EF (I) (four regions outside O) [proof]
- **Uses:** [4], [7], [8], [18], [20], [21]. **Assumptions:** (M1), (H2), σ ≤ 10⁻⁴, M ≥ 112σ².
- **Claim.** g ≥ min{1, c_r/(2√ℓ), 0.14(1 − 2σ), (1 − 2σ)(1 − 0.284)/(4√2k̄)}.
- **Proof.** g ≥ ½((a) + (a′)) and (a) + (a′) ≥ M/(σ(A + B)). (1) A/(A + B) ≥ 2σ. (2) θ ≥ 2σc_r/√ℓ. (3) A/(A + B) < 2σ, θ < 2σc_r/√ℓ, B < K_rσ: M/(σ(A + B)) > 0.28(1 − 2σ). (4) otherwise D_F/(σ(A + B)) ≥ (1 − 2σ)|C̃_w|/(2√2k_eff), C̃_w ≥ 1 − β′, β′ ≤ 0.284, and k_eff ≤ k̄ by [18] and [8]. ∎
- **Numerical.** 481 states outside O with σ ≤ 10⁻⁴ and (H2): g exceeds the bound by a factor ≥ 1809 (J22).

## 24. Lemma EG-D (expansion around the observed point) [proof]
- **Uses:** [1], [2], [3], [4](b). **Assumption:** (M1).
- **Claim.** **σ(A + B)g = 2A_w⁰ + N_φ⁰ + E₂ + R′**, |R′| ≤ (√15/6)σ²√(E s̃_p²) ≤ 0.6455σ²√J, |E₂| ≤ (σ/2)√U·Ξ, U ≤ 4N_φ⁰ + A, A_w⁰ ≥ A/2.
- **Proof.** [3]: N = σ⁻¹E[s̃_p(Y)(f₀(X_z) − f₀(Y_z))], X_z = Y_z − σZ_z. Taylor: f₀(Y_z − σZ_z) − f₀(Y_z) = −σZ_zφ₀ + (σ²/2)Z_z² sin Y_z + ρ_D, |ρ_D| ≤ σ³|Z_z|³/6. First order: E[s̃_ps̃_zφ₀] = ½E[(s̃_v² − s̃_w²)φ₀]. Second order: E[s̃_p sin Y_z] = 0 ([1](i)), so E[Z_z² | Y] may be replaced by W_zz. Third order: Cauchy–Schwarz with E Z⁶ = 15. ∎
- **Numerical.** Without quadrature error on three point masses (σ = 1 to 0.1) and on two Gaussian blobs by plain quadrature (σ = 10⁻², 10⁻⁴): |R′| < 3% of its bound. One of the two reviews of the manuscript (OpenAI ChatGPT) checked it on another three-point state (points (0, 0), (1.3, 0.45), (−0.7, 1.2) in (z, p), weights 0.2, 0.5, 0.3, σ = 0.55): |R′| = 0.00988 ≤ 0.2286; our recomputation (`src/ei_point.py`) gives R′ = 0.00988 and 2A + N = 1.2496652.

## 25. Lemma EG-iso (isotropic three circles) [proof]
- **Uses:** [14], [17]. **Units:** σ = 1, isotropic noise.
- **Claim.** For a unit vector e with I_e(1) ≤ e⁻⁴ and L_e = ln(1/I_e(1)), with δ, T, η, τ, r₁–r₃, ω of [18] built from L_e: **I_e(1 − τ) ≤ 70.2(L_e/ln 2L_e)I_e(1)**.
- **Proof.** The proof of [18] with q_z := ∂_eϱ̃^iso_z = −∫((y − x)·e/z)𝒩^{(2)}_z(y − x)dρ̃. Outer circle: Cauchy–Schwarz and both formulas of [14] give ‖q_z‖² ≤ T³/(T² − r²)², = 16/(9T) ≤ 16/9 at r = T/2. Inner circle: a two-dimensional kernel integral is the square of a one-dimensional one, ≤ I_e(1)/(1 − (r/δ)²) ≤ I_e(1)/η. **[review fix 4]** Holomorphy as in [17], with the dominating function **G_𝒦(y) = C_𝒦∫e^{−c_𝒦|y−x|²}|(y − x)·e|²dρ̃(x)**, integrable by the second formula of [14]. Three circles: ≤ e(16/9)^{1/4}(4L_e/ln 2L_e)I_e(1) ≤ 12.56(L_e/ln 2L_e)I_e(1). Return: ϱ^iso_{1−τ} = (T/(1 − τ))ϱ^iso_TW, q²/ϱ ≤ (T²/(1 − τ)³)ϱ_TE[Z_e² | Y_T], and in two dimensions P(|Z| > K_b) = e^{−K_b²/2} = I_e(1)² exactly. E_bad contributes ≤ √6T²(1 − τ)^{−3}I_e(1) ≤ 2.91I_e(1), the rest ≤ 2e·12.56(L_e/ln 2L_e)I_e(1), and L_e/ln 2L_e ≥ 4/ln 8 = 1.92359…; the total coefficient is ≤ 69.80 ≤ 70.2. ∎

## 26. Lemma EG-sm (diagonal smoothing) [proof]
- **Claim.** For every u ∈ (0, 1), E[(∂_e²ϱ̃^iso_1/ϱ̃^iso_1)²] ≤ I_e(1 − u)/u; if I_e(1) ≤ e⁻⁴ this is **≤ h_O(I_e(1))**.
- **Proof.** Y = Y′ + √uZ₂ with Y′ ~ ϱ^iso_{1−u}; ∂_e²ϱ₁/ϱ₁ = −u^{−1/2}E[Z_{2,e}∂_e ln ϱ_{1−u}(Y′) | Y]; Jensen and independence; u = τ and 70.2·8 ≤ 562. ∎
- **Joint.** In units σ = 1, W_zz = ∂_z²ϱ̃₁/ϱ̃₁ ([2]), so with e = e_z: Ξ² ≤ h_O(I_z(1)), I_z(1) = E[s̃_z²] ≤ J.

## 27. Theorem EG-O (small information) [proof]
- **Uses:** [8], [24], [26]. **Assumptions:** (M1), (H2), J = A + B ≤ J₀ = 1.9×10⁻⁸.
- **Claim.** |E₂ + R′| ≤ ¼(2A_w⁰ + N_φ⁰), hence **g ≥ ¾A/(σ(A + B))**.
- **Proof.**
  1. **[review fix 5]** *The size of σ.* [8] gives A ≥ ασ², so σ² ≤ J/α = 2J(π² + σ²), i.e. **σ² ≤ 2π²J/(1 − 2J)** ≤ 3.7505×10⁻⁷, and α ≥ 0.0506605.
  2. **[review fix 5]** *The second derivative.* In units σ = 1, W_zz = ∂_z²ϱ̃₁/ϱ̃₁ with e = e_z, and I_z(1) = E[s̃_z²] ≤ J. By (H2) and [8], **I_z(1) ≥ 2ασ² > 0**, so L_z = ln(1/I_z(1)) is finite. h_O increases on (0, e⁻⁴] (d ln h_O/d ln x = 1 − 3/L + 2/(L ln 2L) > 0 for L ≥ 4), so [26] gives Ξ² ≤ h_O(J₀) ≤ 4.706×10⁻³ and Ξ ≤ 0.06860.
  3. |E₂| ≤ σΞ√N_φ⁰ + (σΞ/2)√A ≤ ¼N_φ⁰ + (Ξ²/α + Ξ/(2√α))A, using σ² ≤ A/α; |R′| ≤ 0.6455σ²√J ≤ (0.6455√J₀/α)A.
  4. The bracket is 0.09288 + 0.15238 + 0.00176 ≤ 0.24702 ≤ ¼ (interval arithmetic). As A ≤ 2A_w⁰, |E₂ + R′| ≤ ¼(N_φ⁰ + A) ≤ ¼(2A_w⁰ + N_φ⁰); insert into [24]. ∎
- **Numerical.** On the 36 recorded states with J ≤ 2×10⁻⁹ and (H2), |r⁰| ≤ 1.2×10⁻⁶; on all 200 engine states inside O, ≤ 2.9×10⁻⁵ (J29, J30).

## 28. Theorem EG [proof]
- **Uses:** [8], [18], [22], [23], [27]. **Assumptions:** (M1), (H2).
- **Claim.** (a) For 0 < σ ≤ 10⁻⁴ and M ≥ 112σ²: g ≥ min{1, c_r/(2√ℓ), 0.14(1 − 2σ), (1 − 2σ)(1 − 0.284)/(4√2k̄)}. (b) For 0 < σ ≤ 1.5×10⁻¹⁴, for every state: g ≥ ¾α/(ασ + 1268k̄), the smallest of the five terms. (c) For 0 < σ ≤ 10⁻¹⁵: **g ≥ 2.9×10⁻⁵/k̄ ≥ 1.7×10⁻⁶(ln 2L̄)^{3/4}/L̄^{5/4}**.
- **Proof.**
  1. (a) is [23].
  2. (b) inside O: [22] and k_eff ≤ k̄ give B ≤ 1268σk̄ (k̄ ≥ 1), and A < 112σ², so J ≤ 112σ² + 1268σk̄. **[review fix 6]** The right side increases with σ: α depends on σ, **dL̄/dℓ = 2 − 2σ²/(π² + σ²) ∈ (0, 2)**, and **d ln(σk̄)/dℓ = −1 + ½(ψ′(L̄)/ψ(L̄))(dL̄/dℓ) < −1 + ψ′(L̄)/ψ(L̄) < 0**, because ψ′/ψ < 2.5/L ≤ 1 for L ≥ 4. At σ = 1.5×10⁻¹⁴, L̄ = 66.644, k̄ ≤ 984.16 and J ≤ 1.8719×10⁻⁸ ≤ J₀. [27] gives g ≥ ¾A/(σ(A + B)) ≥ ¾α/(ασ + 1268k̄).
  3. The smallest term: for k̄ ≥ 984 the bound of (b) is below 3.1×10⁻⁵/k̄, the fourth term of (a) is above 0.12/k̄, and the other three do not depend on k̄.
  4. (c): at σ = 10⁻¹⁵, L̄ = 72.060, k̄ ≤ 1072.3, ¾α/(ασ + 1268k̄) ≥ 2.996×10⁻⁵/k̄. For L̄ ≥ 72, (4L̄ + 3.57)(ln 2L̄)^{3/2}/L̄^{5/2} ≤ 0.0735 decreases, k̄ ≤ 17.01L̄^{5/4}/(ln 2L̄)^{3/4}, and 0.75α/(1268·17.01) ≥ 1.76×10⁻⁶. ∎
- **Numerical.** J22 (part (a) on 481 states) and J29, J30 (the conclusion inside O on 36 states). The recorded states inside O all have σ ≥ 10⁻⁶, so the step "J ≤ J₀ inside O" is checked by interval arithmetic only.

---

## 29. The chain (which item is used where)

```
[1] integration by parts -+- [2] Tweedie -+- [3] gap formula - [4] branch identity -+- [5] F's G2 - [6] EF-1 - [7] EF-2 --+
                          |               |                                        |                                     |
                          |               +- [9] ED lemma 0 ---------------------- +-------------+                       |
                          |               +- [10] half lemma --+                   |             |                       |
                          |               +- [11] ED lemma 3 - [12] Q-id --+       |             |                       |
                          |               +- [13] Q-sm --------------------+       |             |                       |
                          |  [14]->[15],[16]->[17] three circles ----------+- [18] Theorem EG-1 (k_eff^2 <= psi) -+      |
                          +- [8] Cramer-Rao floor (H2) --------------------------------------------------------+       |
                          +- [19] Heisenberg step (sigma^2) - [21] EF-3 -+- [22] EF-O -----------------------+          |
                          +- [20] EE-C improved -----------------------+                                     |          |
                          |                                         [23] Theorem EF (I): outside O <---------+----------+
                          +- [24] EG-D -+  [14],[17]->[25] EG-iso -> [26] EG-sm -+
                                        +----------------------------------------+- [27] Theorem EG-O: J <= 1.9e-8
                                                                   [28] Theorem EG = [23] (outside O) + [22]+[27] (inside O)
```

## 30. The numerical constants (interval arithmetic)
- 60-digit decimal arithmetic with outward rounding; √, exp, ln correctly rounded and widened by one unit in the last place on each side; π as an interval of its decimal expansion (`src/eh_iv.py`).
- **All 101 claims hold** (72 in `src/eh_constants.py`, 29 in `src/ei_constants.py`; Table 4 of the paper, generated from the records).
- **Three rounded values printed in the working record were rounded in the wrong direction for the way they were used.** The constants derived from them (36.1, 70.2) hold when recomputed. (1) √(4/ln 8) = 1.38693…, printed as 1.387 and used as a lower bound (item 18); (2) 4/ln 8 = 1.92359…, printed as 1.924 (item 25); (3) c₁ = 9.71673… of the original EE-C, printed as 9.7167 (not used by the chain).

## 31. The joints on recorded states
34 joints, each inserting the conclusion of one item into the hypothesis of the next on the recorded states (692 engine states of paper F's code, 2449 two-dimensional and 2921 one-variable states with closed-form fields): `src/eh_joints.py`, record `results_EH/eh_joints.json`, Table 5 of the paper. All hold; J18, the window inequality without σ², fails as it must.

---

## Map: items of this document → the paper

| item | paper | item | paper |
|---|---|---|---|
| [1] | Lemma A1 | [16] | Lemma 6 (ii) |
| [2] | Lemma A2 | [17] | Lemma 7 |
| [3] | Proposition A3 | [18] | Theorem 8, Remark 9 |
| [4] | Proposition A4 | [19] | Lemma 12 |
| [5] | Lemma B1 | [20] | Lemma 13 |
| [6] | Lemma 10 | [21] | Lemma 12 |
| [7] | Lemma 11 | [22] | Lemma 14 |
| [8] | Lemma A5 | [23] | Proposition 15 |
| [9]–[12] | Lemma 3 (proof in App. B) | [24] | Lemma 16 |
| [13] | Lemma 4 | [25] | Lemma 17 |
| [14] | Lemma 5 | [26] | Lemma 18 |
| [15] | Lemma 6 (i) | [27] | Theorem 19 |
| | | [28] | Theorem 1 |
