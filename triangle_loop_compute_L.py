# -*- coding: utf-8 -*-
"""
triangle_loop_compute_L.py
=====================================================================
Compute the azimuthally-averaged Triangle angular factors

              \overline{\mathcal{E}}^L_{\sigma1\sigma2}

for the physical Triangle DIAGRAM (one hard + two soft propagators),
for the fixed external helicity lambda = L  on the single hard line k1.

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
Angular factor (lambda = L)  -- THREE factor product:
    E^L_{σ1σ2} = A^L_{σ1} * B_{σ1σ2} * C^L_{σ2}
        A = e_{k1}^{(L)*} · e_K^{(σ1)}
        B = e_K^{(σ1)*}  · e_q^{(σ2)}
        C = e_q^{(σ2)*}  · e_{k1}^{(L)}

  (cf. the .tex:  E^L_{σ1σ2} = A_{Lσ1} B_{σ1σ2} C_{σ2 L} )

IMPORTANT (same convention as box_loop_compute.py):  EACH factor is kept
as a raw polynomial in Cphi, Sphi -- NOT reduced, NOT simplified, NOT
azimuthally averaged.  The phi-average is applied separately below.

---------------------------------------------------------------------
Only the longitudinal hard line (lambda = L) is requested, so this is
lambda1 fixed to 'L'.  The nine soft-helicity combinations
(σ1,σ2) in {L,+,-} x {L,+,-}  are computed.
The final weighted sum applies, to each soft line, one SYMBOLIC beta
factor for a longitudinal (L) line and none for transverse (+/-) lines:
    β_L = beta,   β_+ = β_- = 1,
with NO eta sign weights (helicity expansion all-plus, same as the Box
program).  The result is
    S(L) = sum_{σ1,σ2} β_{σ1} β_{σ2} \mathcal{I}^L_{σ1σ2},
where \mathcal{I}^L_{σ1σ2} is the unweighted J-seed expansion.
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

def conj(m):
    return sp.conjugate(m)

# ----------------------------------------------------------------------
# build E^L_{σ1σ2} = A * B * C   (lambda fixed to L)
# ----------------------------------------------------------------------
def build_E(sig1, sig2):
    """
    sig1 = soft helicity on the K-line  (in {'L', +1, -1})
    sig2 = soft helicity on the q-line   (in {'L', +1, -1})
    External hard line lambda = L is fixed.
    """
    # A = e_{k1}^{(L)*} · e_K^{(σ1)}    (k1 L is real -> conj no-op)
    eK_s1 = eL_K() if sig1 == 'L' else es_K(sig1)
    A = (conj(k1_L())).dot(eK_s1)

    # B = e_K^{(σ1)*} · e_q^{(σ2)}
    eq_s2 = eL_q() if sig2 == 'L' else es_q(sig2)
    B = (conj(eK_s1)).dot(eq_s2)

    # C = e_q^{(σ2)*} · e_{k1}^{(L)}    (k1 L real -> conj no-op)
    C = (conj(eq_s2)).dot(k1_L())

    # three-factor product -> polynomial in Cphi,Sphi
    return sp.expand(A*B*C)

helicity = ['L', 1, -1]                       # order: L, +, -

# ----------------------------------------------------------------------
# compute all 9 factors  (no phi-average yet)
# ----------------------------------------------------------------------
print("Computing (lambda = L)  [plain dot products, no phi-average]:")
results = {}
for sg1 in helicity:
    for sg2 in helicity:
        results[(sg1, sg2)] = build_E(sg1, sg2)
        print("  done  sigma1 = %s, sigma2 = %s" % (sg1, sg2), flush=True)

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

print("Computing phi-averages  \\overline{\\mathcal{E}}^L_{\\sigma1\\sigma2}:")
avg_results = {}
for sg1 in helicity:
    for sg2 in helicity:
        avg_results[(sg1, sg2)] = phi_average(results[(sg1, sg2)])
        print("  averaged  sigma1 = %s, sigma2 = %s" % (sg1, sg2), flush=True)

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
#    Each phi-averaged factor \bar{\mathcal{E}}^L_{σ1σ2} is expressed as
#    a linear combination of loop seeds
#        J_r(m,n) = f_r( x - m/2 , x + n/2 ),    x = -i d \tilde\nu
#    with f_0, f_1(=f_c), f_2(=f_cc) from the .tex loop-seed section.
#    sin(theta) eliminated via sin^2 = 1-cos^2;  z = cos(theta).
#    J_3 is kept explicitly (NOT reduced to J_0).
# ----------------------------------------------------------------------
def classify_term(term):
    """Split one monomial into (const, m, n, r, s, a)   [same as Box code]."""
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
    """One expanded monomial -> list of (coeff, r, m, n); u^s -> z expansion."""
    const, m, n, r, s, a = classify_term(term)
    if s == 0:
        return [(const, r, m, n)]
    polyz = sp.Poly(sp.expand((1 - z**2)**(s // 2)), z)
    out = []
    for (rz,), cf in polyz.terms():
        out.append((sp.expand(const * cf), int(r + rz), m, n))
    return out

def map_to_J(expr):
    """expr -> list of (coeff, r, m, n)  s.t. expr = sum coeff_i J_{r_i}(m_i,n_i)."""
    e = sp.expand(expr)
    terms_list = list(e.args) if e.is_Add else [e]
    out = []
    for t in terms_list:
        out.extend(process_monomial(t))
    return out

Ilabels = {
    ('L','L'):   r'\mathcal{I}^{L}_{LL;f}',
    ('L',1):     r'\mathcal{I}^{L}_{L+;f}',
    ('L',-1):    r'\mathcal{I}^{L}_{L-;f}',
    (1,'L'):     r'\mathcal{I}^{L}_{+L;f}',
    (1,1):       r'\mathcal{I}^{L}_{++;f}',
    (1,-1):      r'\mathcal{I}^{L}_{+-;f}',
    (-1,'L'):    r'\mathcal{I}^{L}_{-L;f}',
    (-1,1):      r'\mathcal{I}^{L}_{-+;f}',
    (-1,-1):     r'\mathcal{I}^{L}_{--;f}',
}

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

print("Mapping phi-averaged factors onto loop seeds J_r(m,n):")
mapped = {}
for sg1 in helicity:
    for sg2 in helicity:
        mapped[(sg1, sg2)] = map_to_J(avg_results[(sg1, sg2)])
        print("  mapped  sigma1 = %s, sigma2 = %s  (%d terms)"
              % (sg1, sg2, len(mapped[(sg1, sg2)])), flush=True)

# ----------------------------------------------------------------------
# write the NINE individual dot-product factors to file.
# Output order follows box_loop_compute.py EXACTLY:
#   1) polarization part        : plain then latex, raw A*B*C (no phi-average)
#   2) phi-averaged part        : plain then latex
#   3) mapped loop-seed part    : plain then latex  (J_r(m,n) expansions)
# ----------------------------------------------------------------------
label = {
    ('L','L'):   r'\overline{\mathcal{E}}^{L}_{LL}',
    ('L',1):     r'\overline{\mathcal{E}}^{L}_{L+}',
    ('L',-1):    r'\overline{\mathcal{E}}^{L}_{L-}',
    (1,'L'):     r'\overline{\mathcal{E}}^{L}_{+L}',
    (1,1):       r'\overline{\mathcal{E}}^{L}_{++}',
    (1,-1):      r'\overline{\mathcal{E}}^{L}_{+-}',
    (-1,'L'):    r'\overline{\mathcal{E}}^{L}_{-L}',
    (-1,1):      r'\overline{\mathcal{E}}^{L}_{-+}',
    (-1,-1):     r'\overline{\mathcal{E}}^{L}_{--}',
}

with open('triangle_loop_E_LambdaL.txt', 'w') as f:
    f.write("="*90 + "\n")
    f.write("Triangle dot-product factors  \\mathcal{E}^L_{\\sigma1\\sigma2}\n")
    f.write("lambda = L (single hard line k1),  sigma1,sigma2 in {L,+,-}\n")
    f.write("Plain A*B*C dot products (no phi-average, no K reduction).\n")
    f.write("Variables: ks, q, K, z(=cos th), u(=sin th), Cphi=cos ph, Sphi=sin ph,\n")
    f.write("           s1=sin th1, c1=cos th1;  K kept independent.\n")
    f.write("="*90 + "\n\n")

    # ---- 1) polarization part : raw dot products (no phi-average) ----
    f.write(">>> PLAIN FORM  (no phi-average) <<<\n")
    for sg1 in helicity:
        for sg2 in helicity:
            f.write(label[(sg1,sg2)] + " =\n")
            f.write(to_plain(results[(sg1,sg2)]))
            f.write("\n\n")

    f.write("\n>>> LATEX FORM  (no phi-average) <<<\n")
    for sg1 in helicity:
        for sg2 in helicity:
            f.write(label[(sg1,sg2)] + " =\n")
            f.write(to_latex(results[(sg1,sg2)]))
            f.write("\n\n")

    # ---- 2) phi-averaged part ----
    f.write("\n\n==============================================================\n")
    f.write("PHI-AVERAGED  \\overline{\\mathcal{E}}^L_{\\sigma1\\sigma2}\n")
    f.write("(phi integrated over [0,2*pi]; K = |q+k_s| kept independent)\n")
    f.write("==============================================================\n\n")

    f.write(">>> PLAIN FORM  (phi-averaged) <<<\n")
    for sg1 in helicity:
        for sg2 in helicity:
            f.write(label[(sg1,sg2)] + " =\n")
            f.write(to_plain(avg_results[(sg1,sg2)]))
            f.write("\n\n")

    f.write("\n>>> LATEX FORM  (phi-averaged) <<<\n")
    for sg1 in helicity:
        for sg2 in helicity:
            f.write(label[(sg1,sg2)] + " =\n")
            f.write(to_latex(avg_results[(sg1,sg2)]))
            f.write("\n\n")

    # ---- 3) mapped loop-seed part : J_r(m,n) expansions ----
    f.write("\n\n==============================================================\n")
    f.write("MAPPING TO LOOP SEEDS  J_r(m,n)\n")
    f.write("Each phi-averaged  \\bar{\\mathcal{E}}^L_{\\sigma1\\sigma2}  is expressed\n")
    f.write("as a linear combination of loop seeds\n")
    f.write("    J_r(m,n) = f_r(x - m/2 , x + n/2),   x = -i d \\tilde\\nu,\n")
    f.write("with f_0, f_c=f_1, f_cc=f_2, f_3 defined in the .tex.\n")
    f.write("sin(theta) eliminated via sin^2 = 1 - cos^2;  z = cos(theta).\n")
    f.write("J_3 is kept explicitly (NOT reduced to J_0).\n")
    f.write("This equals the dimensionless coefficient \\mathcal{I}^L_{\\sigma1\\sigma2}.\n")
    f.write("="*90 + "\n\n")

    f.write(">>> PLAIN FORM (J-seed expansion) <<<\n")
    for sg1 in helicity:
        for sg2 in helicity:
            f.write(Ilabels[(sg1,sg2)] + " =\n")
            f.write(" + ".join(J_plain(co, r, m, n)
                               for (co, r, m, n) in mapped[(sg1,sg2)]) or "0")
            f.write("\n\n")

    f.write("\n>>> LATEX FORM (J-seed expansion) <<<\n")
    for sg1 in helicity:
        for sg2 in helicity:
            f.write(Ilabels[(sg1,sg2)] + " =\n")
            f.write(" + ".join(J_latex(co, r, m, n)
                               for (co, r, m, n) in mapped[(sg1,sg2)]) or "0")
            f.write("\n\n")

print("Nine individual dot-product factors written to triangle_loop_E_LambdaL.txt")
print("(polarization -> phi-average -> J-seed mapping, matching the Box output order)")

# ======================================================================
#  BETA-WEIGHTED SUM  ->  S(L)   (same convention as box_loop_compute.py)
# ======================================================================
#  The main-text Triangle section writes  I^L = sum_{σ1,σ2} I^L_{σ1σ2},
#  where EACH individual  I^L_{σ1σ2}  carries the longitudinal beta-factor
#        β_L = (-1/2 + i \tilde\nu)^2 / (\tilde\nu^2 + 1/4),   β_+ = β_- = 1.
#  The helicity expansion is all-plus (no eta sign weights), exactly as in
#  box_loop_compute.py.  Each soft longitudinal (L) line carries one
#  SYMBOLIC beta factor (no nu dependence); transverse (+/-) lines carry
#  none.  The sum is done at the J-seed level (Scheme A).
#
#        S(L) = sum_{σ1,σ2 in {L,+,-}}  β_{σ1} β_{σ2} I^L_{σ1σ2}
#
#  NO sin^2+cos^2=1 reduction is applied to the coefficients.
# ======================================================================

def beta_weight(sig):
    """One beta factor for a soft L line; 1 for a soft +/- line."""
    return beta if sig == 'L' else sp.Integer(1)

# combined: key (r,m,n) -> summed coefficient  beta_{σ1} beta_{σ2}
combined = {}
for sg1 in helicity:
    for sg2 in helicity:
        w = beta_weight(sg1) * beta_weight(sg2)
        for (coeff, r, m, n) in mapped[(sg1, sg2)]:
            key = (r, m, n)
            combined[key] = combined.get(key, sp.Integer(0)) + w * coeff

order = sorted(combined.keys(), key=lambda k: (k[2], k[1], k[0]))

def Slabel():
    return r"S(L)"

with open('triangle_loop_B.txt', 'w') as f:
    f.write("="*90 + "\n")
    f.write("S(L) = sum_{sg1,sg2 in {L,+,-}}  beta_{sg1} beta_{sg2} I^{L}_{sg1,sg2;f}\n")
    f.write("Each soft longitudinal (L) line carries one beta factor;\n")
    f.write("transverse (+/-) lines carry none (beta_L = beta, beta_+/- = 1).\n")
    f.write("beta is kept symbolic (no nu dependence).  Sum done at the\n")
    f.write("J-seed level (Scheme A).  NO sin^2+cos^2=1 reduction applied.\n")
    f.write("S  =  sum_{(r,m,n)} C_{rmn} * J_r(m,n)\n")
    f.write("="*90 + "\n\n")

    f.write(">>> LATEX FORM <<<\n")
    terms_l = [J_latex(combined[k], k[0], k[1], k[2])
               for k in order if combined[k] != 0]
    f.write(Slabel() + r" = \left[ " + " + ".join(terms_l) + r" \right]")
    f.write("\n\n")

    f.write(">>> PLAIN FORM <<<\n")
    terms_p = [J_plain(combined[k], k[0], k[1], k[2])
               for k in order if combined[k] != 0]
    f.write(Slabel() + " = ( " + " + ".join(terms_p) + " )")
    f.write("\n\n")

print("beta-weighted sum S(L) written to triangle_loop_B.txt")
