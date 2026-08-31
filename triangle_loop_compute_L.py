# -*- coding: utf-8 -*-
"""
triangle_loop_compute_L.py
=====================================================================
Compute the azimuthally-averaged Triangle angular factors

              \overline{\mathcal{E}}^\lambda_{\sigma1\sigma2}

for the physical Triangle DIAGRAM (one hard + two soft propagators),
for all external helicities lambda in {L,+,-} on the single hard line k1.

This follows exactly the same pipeline as box_loop_compute.py, except
that a Triangle factor has only THREE contraction factors (one hard
propagator), whereas the Box has FOUR (two hard propagators).

---------------------------------------------------------------------
Geometry (following the .tex Triangle section conventions)
    k_s || z-axis          :  k_s = k_s (0,0,1)
    k_1 = k_1 (sinθ1, 0, cosθ1)              phi_1 = 0   (hard line)
    q   = q ( sinθ cosφ,   sinθ sinφ,   cosθ )           (soft line)
    K   = q + k_s                            (K is φ-independent)

Notation:  z    = cosθ ,  u = sinθ
           s1  = sinθ1 ,  c1 = cosθ1
           Cphi = cosφ ,  Sphi = sinφ        (φ = azimuth of q)

---------------------------------------------------------------------
Polarization vectors
    soft on q :
        e^{(L)}_q = ( u Cphi, u Sphi, z )
        e^{(σ)}_q = 1/√2 ( z Cphi - σ i Sphi ,
                           z Sphi + σ i Cphi , -u )
    soft on K :
        e^{(L)}_K = 1/K ( q u Cphi, q u Sphi, q z + k_s )
        e^{(σ)}_K = 1/(√2 K) ( (q z+k_s) Cphi - σ i K Sphi,
                               (q z+k_s) Sphi + σ i K Cphi, -q u )
    hard k1 :
        e^{(L)}_{k1} = ( s1, 0, c1 )              (real -> conj no-op)
        e^{(σ)}_{k1} = 1/√2 ( c1, σ i, -s1 )

---------------------------------------------------------------------
Angular factor -- THREE factor product:
    E^lambda_{σ1σ2} = A^lambda_{σ1} * B_{σ1σ2} * C^lambda_{σ2}
        A = e_{k1}^{(lambda)*} · e_K^{(σ1)}
        B = e_K^{(σ1)*}  · e_q^{(σ2)}
        C = e_q^{(σ2)*}  · e_{k1}^{(lambda)}

  (cf. the .tex: E^lambda_{σ1σ2}
       = A_{lambda σ1} B_{σ1 σ2} C_{σ2 lambda})

IMPORTANT (same convention as box_loop_compute.py):  EACH factor is kept
as a raw polynomial in Cphi, Sphi -- NOT reduced, NOT simplified, NOT
azimuthally averaged.  The phi-average is applied separately below.

---------------------------------------------------------------------
All three hard-line helicities and all nine soft-helicity combinations
(σ1,σ2) in {L,+,-} x {L,+,-} are computed (27 factors total).
The final weighted sum applies, to each soft line, one SYMBOLIC beta
factor for a longitudinal (L) line and none for transverse (+/-) lines:
    β_L = beta,   β_+ = β_- = 1,
with NO eta sign weights (helicity expansion all-plus, same as the Box
program).  The result is
    S(lambda) = sum_{σ1,σ2} β_{σ1} β_{σ2}
                \mathcal{I}^lambda_{σ1σ2;d},
where d is the Triangle nonlocal-branch index used in the .tex.
=====================================================================
"""

import sympy as sp

# ----------------------------------------------------------------------
# symbols
# ----------------------------------------------------------------------
q    = sp.Symbol('q',   positive=True)
K    = sp.Symbol('K',   positive=True)      # |q + k_s|  kept as a real symbol
ks   = sp.Symbol('k_s', positive=True)
z    = sp.Symbol('z',   real=True)          # cos(theta)
u    = sp.Symbol('u',   real=True)          # sin(theta)
Cphi = sp.Symbol('Cphi', real=True)         # cos(phi)   (azimuth of q)
Sphi = sp.Symbol('Sphi', real=True)         # sin(phi)

s1   = sp.Symbol('s1', real=True)           # sin(theta1)
c1   = sp.Symbol('c1', real=True)           # cos(theta1)

I    = sp.I

beta = sp.Symbol('beta', commutative=True)   # soft longitudinal-mode factor

# ----------------------------------------------------------------------
# polarization vectors  (polynomials in Cphi, Sphi)
# ----------------------------------------------------------------------
def eL_q():
    """longitudinal polarization of the soft line q"""
    return sp.Matrix([u*Cphi, u*Sphi, z])

def es_q(sig):
    """transverse polarization of the soft line q ; sig = +1 or -1"""
    return sp.Matrix([
        (z*Cphi - sig*I*Sphi)/sp.sqrt(2),
        (z*Sphi + sig*I*Cphi)/sp.sqrt(2),
        -u/sp.sqrt(2)])

def eL_K():
    """longitudinal polarization of the soft line K = q + k_s"""
    return sp.Matrix([q*u*Cphi, q*u*Sphi, q*z + ks])/K

def es_K(sig):
    """transverse polarization of the soft line K ; sig = +1 or -1"""
    return sp.Matrix([
        ((q*z + ks)*Cphi - sig*I*K*Sphi)/(sp.sqrt(2)*K),
        ((q*z + ks)*Sphi + sig*I*K*Cphi)/(sp.sqrt(2)*K),
        -q*u/(sp.sqrt(2)*K)])

def k1_L():
    """longitudinal polarization of the (single) hard external line k1"""
    return sp.Matrix([s1, 0, c1])

def k1_s(sig):
    """transverse polarization of the hard line k1; sig = +1 or -1"""
    return sp.Matrix([c1, sig*I, -s1])/sp.sqrt(2)

def conj(m):
    return sp.conjugate(m)

# ----------------------------------------------------------------------
# build E^lambda_{σ1σ2} = A * B * C
# ----------------------------------------------------------------------
def build_E(lam, sig1, sig2):
    """
    lam  = helicity of the hard k1 line (in {'L', +1, -1})
    sig1 = soft helicity on the K-line  (in {'L', +1, -1})
    sig2 = soft helicity on the q-line   (in {'L', +1, -1})
    """
    k1 = k1_L() if lam == 'L' else k1_s(lam)
    eK_s1 = eL_K() if sig1 == 'L' else es_K(sig1)
    eq_s2 = eL_q() if sig2 == 'L' else es_q(sig2)

    # A = e_{k1}^{(lambda)*} · e_K^{(σ1)}
    A = conj(k1).dot(eK_s1)
    # B = e_K^{(σ1)*} · e_q^{(σ2)}
    B = (conj(eK_s1)).dot(eq_s2)
    # C = e_q^{(σ2)*} · e_{k1}^{(lambda)}
    C = conj(eq_s2).dot(k1)

    # three-factor product -> polynomial in Cphi,Sphi
    return sp.expand(A*B*C)

helicity = ['L', 1, -1]                       # order: L, +, -

def hname(x):
    """helicity name: 'L' -> 'L', +1 -> '+', -1 -> '-'."""
    return {'L': 'L', 1: '+', -1: '-'}[x]

# ----------------------------------------------------------------------
# compute all 3 hard helicities x 9 soft-helicity factors
# ----------------------------------------------------------------------
print("Computing all (lambda, sigma1, sigma2) [plain dot products, no phi-average]:")
results = {}                                  # (lam,sg1,sg2) -> raw E
for lam in helicity:
    for sg1 in helicity:
        for sg2 in helicity:
            results[(lam, sg1, sg2)] = build_E(lam, sg1, sg2)
            print("  done  lambda=%s, sigma1=%s, sigma2=%s"
                  % (hname(lam), hname(sg1), hname(sg2)), flush=True)

# ----------------------------------------------------------------------
# phi-average  (closed-form even-moment averaging, same as Box code)
# ----------------------------------------------------------------------
def ddf(n):
    """double factorial n!! ;  (-1)!! = 1,  0!! = 1."""
    if n <= 0:
        return 1
    r = 1
    while n > 0:
        r *= n
        n -= 2
    return r

def mom_even(m, n):
    """(1/2pi) int cos(phi)^m sin(phi)^n dphi = (m-1)!!(n-1)!!/(m+n)!! for even m,n."""
    if m % 2 == 1 or n % 2 == 1:
        return 0
    return sp.Rational(ddf(m - 1) * ddf(n - 1), ddf(m + n))

def phi_average(expr):
    """Average over azimuth phi of a polynomial in Cphi,Sphi (algebraic)."""
    p = sp.Poly(sp.expand(expr), Cphi, Sphi)
    total = 0
    for (m, n), coeff in p.terms():
        if m % 2 == 0 and n % 2 == 0:
            total += coeff * mom_even(m, n)
    return sp.expand(total)

print("Computing phi-averages for all (lambda, sigma1, sigma2):")
avg_results = {}
for lam in helicity:
    for sg1 in helicity:
        for sg2 in helicity:
            avg_results[(lam, sg1, sg2)] = phi_average(results[(lam, sg1, sg2)])
            print("  averaged  lambda=%s, sigma1=%s, sigma2=%s"
                  % (hname(lam), hname(sg1), hname(sg2)), flush=True)

# ----------------------------------------------------------------------
# output helpers : rename sympy symbols to readable LaTeX / plain text
# ----------------------------------------------------------------------
name_map = {
    s1:  r'\sin\theta_1',   c1:  r'\cos\theta_1',
    z:   r'\cos\theta',     u:   r'\sin\theta',
    Cphi: r'\cos\phi',      Sphi: r'\sin\phi',
    q:   'q',               K:   'K',      ks: 'k_s',
    I:   'i',               beta: r'\beta',
}

def to_latex(expr):
    # plain expansion only; no nsimplify / simplify / phi-average
    latex = sp.latex(sp.expand(expr))
    for sym, rep in name_map.items():
        if sym == I:
            latex = latex.replace("i", rep)
        else:
            latex = latex.replace(sp.latex(sym), rep)
    return latex

def to_plain(expr):
    s = str(sp.expand(expr))
    for sym, rep in name_map.items():
        s = s.replace(str(sym), rep)
    return s

# ----------------------------------------------------------------------
# MAPPING TO LOOP SEEDS  J_r(m,n)
#    Each phi-averaged factor is mapped to
#        J_r(m,n) = f_r(x-m/2, x+n/2),   x = -i d \tilde\nu,
#    where d is the Triangle nonlocal-branch index used in the article.
#    Equal seeds are combined before output.  This is essential for the
#    Triangle: apparent J_4 terms in an uncollected monomial list cancel.
# ----------------------------------------------------------------------
def classify_term(term):
    """Split a monomial into const*q^m*K^(-n)*z^r*u^s*k_s^a."""
    factors = [term] if term.is_Atom else sp.Mul.make_args(term)
    const = sp.Integer(1)
    m = n = r = s = a = 0
    for fac in factors:
        if fac.is_Pow:
            base, exp = fac.as_base_exp()
            if base == q:            m += exp
            elif base == K:          n += -exp
            elif base == z:          r += exp
            elif base == u:          s += exp
            elif base == ks:         a += exp
            else:                    const *= fac
        elif fac.is_Symbol:
            if fac == q:            m += 1
            elif fac == K:          n += -1
            elif fac == z:          r += 1
            elif fac == u:          s += 1
            elif fac == ks:         a += 1
            else:                   const *= fac
        else:
            const *= fac
    return sp.expand(const), int(m), int(n), int(r), int(s), int(a)

def process_monomial(term):
    """Map one monomial to one or more (coeff,r,m,n) seed terms."""
    const, m, n, r, s, a = classify_term(term)
    if s % 2:
        raise ValueError("odd power of sin(theta) survived phi averaging: %s" % term)
    if m - n + a != 0:
        raise ValueError("term has unexpected k_s scaling: %s" % term)
    polyz = sp.Poly(sp.expand((1 - z**2)**(s // 2)), z)
    return [(sp.expand(const*cf), int(r+rz), m, n)
            for (rz,), cf in polyz.terms()]

def map_to_J(expr):
    """Return {(r,m,n): coeff}, combining identical seeds exactly."""
    combined = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        for coeff, r, m, n in process_monomial(term):
            key = (r, m, n)
            combined[key] = combined.get(key, sp.Integer(0)) + coeff
    return {key: coeff for key, value in combined.items()
            if (coeff := sp.expand(value)) != 0}

def Elabel(lam, sig1, sig2):
    return r"\overline{\mathcal{E}}^{%s}_{%s%s}" % (
        hname(lam), hname(sig1), hname(sig2))

def Ilabel(lam, sig1, sig2):
    return r"\mathcal{I}^{%s}_{%s%s;d}" % (
        hname(lam), hname(sig1), hname(sig2))

def Slabel(lam):
    return r"S(%s)" % hname(lam)

def seed_order(seed_dict):
    return sorted(seed_dict, key=lambda key: (key[2], key[1], key[0]))

def J_plain(coeff, r, m, n):
    c = to_plain(coeff)
    if sp.expand(coeff).is_Add:
        c = "( " + c + " )"
    return "%s * J_%d(%d,%d)" % (c, r, m, n)

def J_latex(coeff, r, m, n):
    c = to_latex(coeff)
    if sp.expand(coeff).is_Add:
        c = r"\left( %s \right)" % c
    return r"%s\,\mathcal{J}_{%d}\left(%d,%d\right)" % (c, r, m, n)

def seed_sum_plain(seed_dict):
    return " + ".join(J_plain(seed_dict[key], *key)
                      for key in seed_order(seed_dict)) or "0"

def seed_sum_latex(seed_dict):
    return " + ".join(J_latex(seed_dict[key], *key)
                      for key in seed_order(seed_dict)) or "0"

print("Mapping phi-averaged factors onto collected loop seeds J_r(m,n):")
mapped = {}                                    # (lam,sg1,sg2) -> seed dictionary
for lam in helicity:
    for sg1 in helicity:
        for sg2 in helicity:
            key = (lam, sg1, sg2)
            mapped[key] = map_to_J(avg_results[key])
            if any(r > 3 for r, m, n in mapped[key]):
                raise AssertionError("Triangle produced a non-cancelling J_r with r > 3")
            print("  mapped  lambda=%s, sigma1=%s, sigma2=%s (%d collected seeds)"
                  % (hname(lam), hname(sg1), hname(sg2), len(mapped[key])), flush=True)

# Polarization checks quoted in the Triangle section of the article.
def flip(sig):
    return sig if sig == 'L' else -sig

for sg1 in helicity:
    for sg2 in helicity:
        if sp.expand(avg_results[('L', sg1, sg2)]
                     - avg_results[('L', flip(sg1), flip(sg2))]) != 0:
            raise AssertionError("longitudinal hard-line helicity symmetry failed")
        if sp.expand(avg_results[(-1, sg1, sg2)]
                     - avg_results[(1, flip(sg1), flip(sg2))]) != 0:
            raise AssertionError("plus/minus hard-line conjugation symmetry failed")

# ----------------------------------------------------------------------
# Output all 27 individual factors.  The file name is retained for backward
# compatibility even though it now contains lambda=L,+,-, not only LambdaL.
# ----------------------------------------------------------------------
with open('triangle_loop_E.txt', 'w') as f:
    f.write("="*90 + "\n")
    f.write("Triangle factors  \\mathcal{E}^lambda_{sigma1 sigma2}\n")
    f.write("lambda, sigma1, sigma2 in {L,+,-}; 27 helicity combinations.\n")
    f.write("The Triangle branch label is d: x = -i d \\tilde\\nu.\n")
    f.write("="*90 + "\n\n")

    for lam in helicity:
        f.write("\n" + "="*70 + "\n")
        f.write("lambda = %s\n" % hname(lam))
        f.write("="*70 + "\n\n")

        f.write(">>> PLAIN FORM (no phi-average) <<<\n")
        for sg1 in helicity:
            for sg2 in helicity:
                key = (lam, sg1, sg2)
                f.write(Elabel(*key) + " =\n" + to_plain(results[key]) + "\n\n")

        f.write(">>> LATEX FORM (no phi-average) <<<\n")
        for sg1 in helicity:
            for sg2 in helicity:
                key = (lam, sg1, sg2)
                f.write(Elabel(*key) + " =\n" + to_latex(results[key]) + "\n\n")

        f.write(">>> PLAIN FORM (phi-averaged) <<<\n")
        for sg1 in helicity:
            for sg2 in helicity:
                key = (lam, sg1, sg2)
                f.write(Elabel(*key) + " =\n" + to_plain(avg_results[key]) + "\n\n")

        f.write(">>> LATEX FORM (phi-averaged) <<<\n")
        for sg1 in helicity:
            for sg2 in helicity:
                key = (lam, sg1, sg2)
                f.write(Elabel(*key) + " =\n" + to_latex(avg_results[key]) + "\n\n")

        f.write(">>> PLAIN FORM (collected J-seed expansion) <<<\n")
        for sg1 in helicity:
            for sg2 in helicity:
                key = (lam, sg1, sg2)
                f.write(Ilabel(*key) + " =\n" + seed_sum_plain(mapped[key]) + "\n\n")

        f.write(">>> LATEX FORM (collected J-seed expansion) <<<\n")
        for sg1 in helicity:
            for sg2 in helicity:
                key = (lam, sg1, sg2)
                f.write(Ilabel(*key) + " =\n" + seed_sum_latex(mapped[key]) + "\n\n")

print("All 27 Triangle factors written to triangle_loop_E.txt")

# ======================================================================
# BETA-WEIGHTED SOFT-HELICITY SUMS S(L), S(+), S(-)
# ======================================================================
def beta_weight(sig):
    return beta if sig == 'L' else sp.Integer(1)

def combined_S(lam):
    combined = {}
    for sg1 in helicity:
        for sg2 in helicity:
            weight = beta_weight(sg1)*beta_weight(sg2)
            for key, coeff in mapped[(lam, sg1, sg2)].items():
                combined[key] = combined.get(key, sp.Integer(0)) + weight*coeff
    return {key: coeff for key, value in combined.items()
            if (coeff := sp.expand(value)) != 0}

def weighted_I_sum_latex(lam):
    """Explicit nine-term sum before collecting equal J_r(m,n) seeds."""
    terms = []
    for sg1 in helicity:
        for sg2 in helicity:
            weight = beta_weight(sg1)*beta_weight(sg2)
            factor = "" if weight == 1 else to_latex(weight) + r"\,"
            terms.append(factor + Ilabel(lam, sg1, sg2))
    return " + ".join(terms)

def weighted_I_sum_plain(lam):
    """Plain-text counterpart of the explicit nine-term weighted sum."""
    terms = []
    for sg1 in helicity:
        for sg2 in helicity:
            weight = beta_weight(sg1)*beta_weight(sg2)
            factor = "" if weight == 1 else to_plain(weight) + " * "
            terms.append(factor + Ilabel(lam, sg1, sg2))
    return " + ".join(terms)

all_sums = {lam: combined_S(lam) for lam in helicity}
if all_sums[1] != all_sums[-1]:
    raise AssertionError("article relation S(+) = S(-) failed")

# Check S(+) against the three grouped formulas displayed in the article:
# I^+_LL, the four mixed-soft terms, and the four transverse-soft terms.
# The article uses s1^2+c1^2=1 in several coefficients; the generated output
# deliberately keeps s1 and c1 independent, as in box_loop_compute.py.
expected_plus = {
    (0, 0, 0): (1+c1**2)/4,
    (0, 0, 2): beta*(2+s1**2)/4,
    (2, 0, 2): beta**2*s1**2/2 - beta*(2+s1**2)/4 + (1+c1**2)/4,
    (1, 1, 2): beta**2*((1+c1**2)/4+s1**2/2)
                 + (1+c1**2)/4+s1**2/2,
    (3, 1, 2): beta**2*(-(1+c1**2)/4+s1**2/2)
                 + (1+c1**2)/4-s1**2/2,
    (0, 2, 2): beta**2*(1+c1**2)/4+s1**2/2,
    (2, 2, 2): beta**2*(-(1+c1**2)/4+s1**2/2)
                 + (1+c1**2)/4-s1**2/2,
}
if set(all_sums[1]) != set(expected_plus):
    raise AssertionError("S(+) seed support disagrees with the article")
for key, expected in expected_plus.items():
    difference = sp.expand((all_sums[1][key] - expected).subs(c1**2, 1-s1**2))
    if difference != 0:
        raise AssertionError("S(+) coefficient disagrees with the article at %s" % (key,))

transverse_sum = {
    key: sp.expand(all_sums[1].get(key, 0) + all_sums[-1].get(key, 0))
    for key in set(all_sums[1]) | set(all_sums[-1])
}
transverse_sum = {key: value for key, value in transverse_sum.items() if value != 0}

with open('triangle_loop_B.txt', 'w') as f:
    f.write("="*90 + "\n")
    f.write("S(lambda) = sum_{sigma1,sigma2} beta_sigma1 beta_sigma2 "
            "I^lambda_{sigma1 sigma2;d}\n")
    f.write("lambda in {L,+,-}; beta_L=beta and beta_+=beta_-=1.\n")
    f.write("Equal J_r(m,n) seeds are collected; d is the Triangle branch index.\n")
    f.write("="*90 + "\n\n")

    for lam in helicity:
        f.write(Slabel(lam) + "\n")
        f.write(">>> EXPLICIT NINE-TERM BETA-WEIGHTED DEFINITION <<<\n")
        f.write(Slabel(lam) + " = " + weighted_I_sum_latex(lam) + "\n")
        f.write(Slabel(lam) + " = " + weighted_I_sum_plain(lam) + "\n\n")
        f.write(">>> LATEX FORM <<<\n")
        f.write(Slabel(lam) + r" = \left[ " + seed_sum_latex(all_sums[lam])
                + r" \right]" + "\n\n")
        f.write(">>> PLAIN FORM <<<\n")
        f.write(Slabel(lam) + " = ( " + seed_sum_plain(all_sums[lam]) + " )\n\n")

    label_T = "S(+) + S(-)"
    f.write(label_T + "\n")
    f.write(">>> LATEX FORM <<<\n")
    f.write(label_T + r" = \left[ " + seed_sum_latex(transverse_sum)
            + r" \right]" + "\n\n")
    f.write(">>> PLAIN FORM <<<\n")
    f.write(label_T + " = ( " + seed_sum_plain(transverse_sum) + " )\n")

print("S(L), S(+), S(-), and S(+)+S(-) written to triangle_loop_B.txt")
