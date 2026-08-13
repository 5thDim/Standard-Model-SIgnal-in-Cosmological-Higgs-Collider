# -*- coding: utf-8 -*-
"""
compute_box_loop_LambdaLL.py
=====================================================================
Compute the azimuthally-averaged Box angular factors

        \overline{\mathcal{E}}^{LL}_{\sigma1\sigma2}

for the physical Box DIAGRAM (two hard + two soft propagators),
for the fixed external helicity pair  lambda1 = lambda2 = L.

---------------------------------------------------------------------
Geometry (following the .tex conventions)
    k_s || z-axis          :  k_s = k_s (0,0,1)
    k_1 = k_1 (sinθ1, 0, cosθ1)              phi_1 = 0
    k_3 = k_3 (sinθ3 cosφ3, sinθ3 sinφ3, cosθ3)
    q   = q ( sinθ cosφ,   sinθ sinφ,   cosθ )
    K   = q + k_s                            (K is φ-independent)

Notation:  z    = cosθ ,  u = sinθ
           s1  = sinθ1 , c1 = cosθ1
           s3  = sinθ3 , c3 = cosθ3
           Cphi = cosφ , Sphi = sinφ        (φ = azimuth of q)
           cp  = cosφ3 , sp3 = sinφ3

---------------------------------------------------------------------
Polarization vectors
    soft on q :
        e^{(L)}_q = ( u Cphi, u Sphi, z )
        e^{(σ)}_q = 1/√2 ( z Cphi - σ i Sphi ,
                           z Sphi + σ i Cphi , -u )
    soft on K (from Bubble section of the tex) :
        e^{(L)}_K = 1/K ( q u Cphi, q u Sphi, q z + k_s )
        e^{(σ)}_K = 1/(√2 K) ( (q z+k_s) Cphi - σ i K Sphi,
                               (q z+k_s) Sphi + σ i K Cphi, -q u )
    hard k1 :
        e^{(L)}_{k1} = ( s1, 0, c1 )
    hard k3 :
        e^{(L)}_{k3} = ( s3 cp, s3 sp3, c3 )

---------------------------------------------------------------------
Angular factor (Appendix B form; physical, σ1,σ2 kept explicit)
    E^{LL}_{σ1σ2} = A^{L}_{σ1} * B^{L}_{σ1} * C^{L}_{σ2} * D^{L}_{σ2}
        A = e_{k1}^{(L)*} · e_K^{(σ1)}      (k1 L is real -> no conj)
        B = e_K^{(σ1)*}  · e_{k3}^{(L)}
        C = e_{k3}^{(L)*} · e_q^{(σ2)}      (k3 L is real -> no conj)
        D = e_q^{(σ2)*}  · e_{k1}^{(L)}

---------------------------------------------------------------------
IMPORTANT (as requested):  EACH factor is kept as a raw polynomial in
Cphi, Sphi (plus the phase-space variables) -- NOT reduced, NOT
simplified, NOT azimuthally averaged.  The phi-average is applied
separately below to produce \overline{\mathcal{E}}^{LL}_{\sigma1\sigma2}.
"""

import sympy as sp

# ----------------------------------------------------------------------
# symbols
# ----------------------------------------------------------------------
q    = sp.Symbol('q',   positive=True)
K    = sp.Symbol('K',   positive=True)      # |q + k_s|   kept as a real symbol
ks   = sp.Symbol('k_s', positive=True)
z    = sp.Symbol('z',   real=True)          # cos(theta)
u    = sp.Symbol('u',   real=True)          # sin(theta)
Cphi = sp.Symbol('Cphi', real=True)         # cos(phi)   (azimuth of q)
Sphi = sp.Symbol('Sphi', real=True)         # sin(phi)

s1   = sp.Symbol('s1', real=True)           # sin(theta1)
c1   = sp.Symbol('c1', real=True)           # cos(theta1)
s3   = sp.Symbol('s3', real=True)           # sin(theta3)
c3   = sp.Symbol('c3', real=True)           # cos(theta3)
cp   = sp.Symbol('cp', real=True)           # cos(phi3)
sp3  = sp.Symbol('sp3', real=True)          # sin(phi3)

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
    """longitudinal polarization of hard external line k1  (phi1 = 0)"""
    return sp.Matrix([s1, 0, c1])

def k3_L():
    """longitudinal polarization of hard external line k3"""
    return sp.Matrix([s3*cp, s3*sp3, c3])

def k1_s(sig):
    """transverse polarization of hard external line k1 (phi1=0); sig=+1/-1"""
    return sp.Matrix([c1, sig*I, -s1])/sp.sqrt(2)

def k3_s(sig):
    """transverse polarization of hard external line k3 (phi3); sig=+1/-1"""
    return sp.Matrix([
        (c3*cp - sig*I*sp3)/sp.sqrt(2),
        (c3*sp3 + sig*I*cp)/sp.sqrt(2),
        -s3/sp.sqrt(2)])

def conj(m):
    return sp.conjugate(m)

# ----------------------------------------------------------------------
# build E^{λ1λ2}_{σ1σ2} = A * B * C * D
# ----------------------------------------------------------------------
def build_E(lam1, lam2, sig1, sig2):
    """
    lam1, lam2 = external (hard) helicities on k1, k3  (in {'L', +1, -1})
    sig1, sig2 = soft helicities on the K-line, q-line (in {'L', +1, -1})
    """
    k1 = k1_L() if lam1 == 'L' else k1_s(lam1)
    k3 = k3_L() if lam2 == 'L' else k3_s(lam2)
    eK_s1 = eL_K() if sig1 == 'L' else es_K(sig1)
    eq_s2 = eL_q() if sig2 == 'L' else es_q(sig2)
    # A = e_{k1}^{(λ1)*} · e_K^{(σ1)}
    A = (conj(k1)).dot(eK_s1)
    # B = e_K^{(σ1)*} · e_{k3}^{(λ2)}
    B = (conj(eK_s1)).dot(k3)
    # C = e_{k3}^{(λ2)*} · e_q^{(σ2)}
    C = (conj(k3)).dot(eq_s2)
    # D = e_q^{(σ2)*} · e_{k1}^{(λ1)}
    D = (conj(eq_s2)).dot(k1)
    # four-factor product -> polynomial in Cphi,Sphi
    return sp.expand(A*B*C*D)

helicity = ['L', 1, -1]                       # order: L, +, -

def hname(x):
    """helicity name: 'L' -> 'L', +1 -> '+', -1 -> '-'."""
    return { 'L': 'L', 1: '+', -1: '-' }[x]

# ----------------------------------------------------------------------
# compute all 9 external-helicity combinations x 9 soft combinations
# ----------------------------------------------------------------------
ext = ['L', 1, -1]                            # external helicities, order: L, +, -

all_results = {}                              # (lam1,lam2,sg1,sg2) -> raw E
print("Computing all (lambda1, lambda2, sigma1, sigma2)  [plain dot products, no phi-average]:")
for lam1 in ext:
    for lam2 in ext:
        for sg1 in helicity:
            for sg2 in helicity:
                all_results[(lam1, lam2, sg1, sg2)] = build_E(lam1, lam2, sg1, sg2)
                print("  done  lam1=%s, lam2=%s, sigma1=%s, sigma2=%s"
                      % (hname(lam1), hname(lam2), hname(sg1), hname(sg2)), flush=True)

# ----------------------------------------------------------------------
# phi-average : produce  \overline{\mathcal{E}}^{LL}_{sigma1 sigma2}
#   For each plain dot product (a polynomial in Cphi, Sphi) we substitute
#   Cphi -> cos(phi), Sphi -> sin(phi) and integrate phi over [0, 2*pi],
#   dividing by 2*pi.  Since K = |q + k_s| depends only on the polar angle
#   z = cos(theta) of q (not on the azimuth phi), K factors out and the
#   trigonometric-polynomial integral is exact (all odd-power terms vanish).
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
    """(1/2pi) int_0^{2pi} cos(phi)^m sin(phi)^n dphi
       = (m-1)!!(n-1)!!/(m+n)!!  for even m,n ;
       0 if either m or n is odd (defensive, though the caller
       already filters to even,m even,n)."""
    if m % 2 == 1 or n % 2 == 1:
        return 0
    return sp.Rational(ddf(m - 1) * ddf(n - 1), ddf(m + n))

phi = sp.Symbol('phi', real=True)   # kept only for label clarity; not needed for averaging

def phi_average(expr):
    """
    Average over azimuth phi of a polynomial in Cphi, Sphi, WITHOUT integration:
    treat each monomial  Cphi^m Sphi^n  via the closed-form even-moment formula
    (0 if m or n is odd).  K = |q + k_s| is phi-independent and merely factors out
    as a constant coefficient, so the result is an ordinary algebraic expression.
    """
    p = sp.Poly(sp.expand(expr), Cphi, Sphi)
    total = 0
    for (m, n), coeff in p.terms():
        if m % 2 == 0 and n % 2 == 0:
            total += coeff * mom_even(m, n)
    return sp.expand(total)

print("Computing phi-averages for all (lambda1, lambda2, sigma1, sigma2):")
all_avg = {}                                  # (lam1,lam2,sg1,sg2) -> phi-averaged E
for lam1 in ext:
    for lam2 in ext:
        for sg1 in helicity:
            for sg2 in helicity:
                all_avg[(lam1, lam2, sg1, sg2)] = phi_average(all_results[(lam1, lam2, sg1, sg2)])
                print("  averaged  lam1=%s, lam2=%s, sigma1=%s, sigma2=%s"
                      % (hname(lam1), hname(lam2), hname(sg1), hname(sg2)), flush=True)

# ----------------------------------------------------------------------
# output helpers : rename sympy symbols to readable LaTeX / plain text
# ----------------------------------------------------------------------
name_map = {
    s1:  r'\sin\theta_1',   c1:  r'\cos\theta_1',
    s3:  r'\sin\theta_3',   c3:  r'\cos\theta_3',
    cp:  r'\cos\phi_3',     sp3: r'\sin\phi_3',
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
# write output file (box_loop_E.txt)
# ----------------------------------------------------------------------
def Elabel(lam1, lam2, sig1, sig2):
    return r"\overline{\mathcal{E}}^{%s%s}_{%s%s}" % (hname(lam1), hname(lam2), hname(sig1), hname(sig2))

with open('box_loop_E.txt', 'w') as f:
    f.write("==============================================================\n")
    f.write("Box dot-product factors  \\mathcal{E}^{lambda1 lambda2}_{sigma1 sigma2}\n")
    f.write("lambda1, lambda2, sigma1, sigma2 in {L,+, -}.\n")
    f.write("Plain A*B*C*D dot products (no phi-average, no K reduction).\n")
    f.write("Variables: ks, q, K, z(=cos th), u(=sin th), Cphi=cos ph, Sphi=sin ph,\n")
    f.write("           k1(s1,c1), k3(s3,c3,cp,sp3);  K kept independent.\n")
    f.write("==============================================================\n\n")

    for lam1 in ext:
        for lam2 in ext:
            f.write("\n" + "="*70 + "\n")
            f.write("lambda1 = %s ,  lambda2 = %s\n" % (hname(lam1), hname(lam2)))
            f.write("="*70 + "\n\n")

            f.write(">>> PLAIN FORM  (no phi-average) <<<\n")
            for sg1 in helicity:
                for sg2 in helicity:
                    f.write(Elabel(lam1, lam2, sg1, sg2) + " =\n")
                    f.write(to_plain(all_results[(lam1, lam2, sg1, sg2)]))
                    f.write("\n\n")

            f.write("\n>>> LATEX FORM  (no phi-average) <<<\n")
            for sg1 in helicity:
                for sg2 in helicity:
                    f.write(Elabel(lam1, lam2, sg1, sg2) + " =\n")
                    f.write(to_latex(all_results[(lam1, lam2, sg1, sg2)]))
                    f.write("\n\n")

            f.write("\n>>> PLAIN FORM  (phi-averaged) <<<\n")
            for sg1 in helicity:
                for sg2 in helicity:
                    f.write(Elabel(lam1, lam2, sg1, sg2) + " =\n")
                    f.write(to_plain(all_avg[(lam1, lam2, sg1, sg2)]))
                    f.write("\n\n")

            f.write("\n>>> LATEX FORM  (phi-averaged) <<<\n")
            for sg1 in helicity:
                for sg2 in helicity:
                    f.write(Elabel(lam1, lam2, sg1, sg2) + " =\n")
                    f.write(to_latex(all_avg[(lam1, lam2, sg1, sg2)]))
                    f.write("\n\n")

print("\nAll done (raw dot products + phi-averaged for all 9 external pairs). "
      "Results written to box_loop_E.txt")

# ======================================================================
#  MAPPING TO LOOP SEEDS  J_r(m,n)
# ======================================================================
#  Next step: take each phi-averaged factor  \bar{\mathcal{E}}^{LL}_{\sigma1\sigma2}
#  (an algebraic expression in q, K, k_s, z=cos th, u=sin th, s1,c1,s3,c3,cp,sp3)
#  and map it onto a linear combination of loop seeds
#
#        J_r(m,n)  =  f_r( x - m/2 ,  x + n/2 ),      x = -i f \tilde\nu
#
#  with  f_0=f_0, f_1=f_c, f_2=f_cc, f_3, f_4 exactly as defined in the .tex.
#
#  PROCEDURE (as requested):
#    (1) eliminate sin(theta)  :  the theta-dependence is polynomial in sin,cos
#        with even sin powers (0,2,4);  use sin^2 = 1 - cos^2  to convert
#        everything into cos(theta)=z.
#    (2) expand the result into monomials   coeff * q^m * K^{-n} * z^r     and
#        map each monomial to  coeff * J_r(m,n).
#    (3) keep the result in terms of J_0,J_1,J_2,J_3,J_4 -- do NOT reduce the
#        J_3,J_4 seeds back to J_0 via the closed forms.
#
#  The  f  branch index enters through x = -i f \tilde\nu, which is kept as a
#  symbolic seed label shared by f = + and f = -.  We do NOT substitute the
#  numerical \tilde\nu.  The output is the dimensionless coefficient
#  \mathcal{I}^{lambda1 lambda2}_{sigma1 sigma2;f} from the .tex.
#
#  IMPORTANT: we only APPEND ("a" mode) to box_loop_E.txt; all the
#  existing content above is left exactly unchanged.
# ======================================================================

def classify_term(term):
    """
    Split a single expanded monomial `term` into
        term  =  const * q^m * K^{-n} * z^r * u^s * (k_s)^a
    by classifying each multiplicative factor.
    const contains ONLY the angular parameters (s1,c1,s3,c3,cp,sp3) together
    with pure numbers and sympy I  -- i.e. it does NOT contain q, K, z, u or k_s.
    Returns  (const, m, n, r, s, a)  with  m,n,r,s,a integers.
    """
    factors = [term] if term.is_Atom else sp.Mul.make_args(term)
    const = sp.Integer(1)
    m = n = r = s = a = 0
    for fac in factors:
        if fac.is_Pow:
            base, exp = fac.as_base_exp()
            if base == q:            m += exp
            elif base == K:          n += -exp    # K^{+e} -> K^{-n}, n = -e
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
        else:                        # numbers, I, sqrt(2), etc.
            const *= fac
    return sp.expand(const), int(m), int(n), int(r), int(s), int(a)


def process_monomial(term):
    """
    Turn ONE expanded monomial into a list of  (coeff, r, m, n).
    The m,n (q and K exponents) are read off BEFORE touching theta, because
    q and K never mix with theta.  The sin/cos of theta (u = sin th, z = cos th)
    is then converted via  sin^2 = 1 - cos^2;  each u^s expands into several
    z^{2k}, so a single monomial gives SEVERAL different r (one per J_r seed).
    The angular/multiplying coefficient is returned by classify_term, and
    k_s is dropped: q, K and k_s are mutually dimension-matched, so the final
    J combination stays dimensionless (as you requested).
    """
    const, m, n, r, s, a = classify_term(term)
    if s == 0:
        return [(const, r, m, n)]
    # u^s = (1 - z^2)^(s/2),   s assumed even (0,2,4,...)
    polyz = sp.Poly(sp.expand((1 - z**2)**(s // 2)), z)
    out = []
    for (rz,), cf in polyz.terms():
        out.append((sp.expand(const * cf), int(r + rz), m, n))
    return out


def map_to_J(expr):
    """
    Return a list of  (coeff, r, m, n)  such that
        expr  =  sum_i coeff_i * J_{r_i}(m_i, n_i),
    where each monomial  q^{m} K^{-n} z^{r}  maps onto  J_{r}(m,n).
    coeff depends only on the angular parameters s1,c1,s3,c3,cp,sp3 (and pure
    numbers / I) -- it contains no q, K, z, u, or k_s.
    """
    e = sp.expand(expr)
    terms_list = list(e.args) if e.is_Add else [e]
    out = []
    for t in terms_list:
        out.extend(process_monomial(t))
    return out

def Ilabel(lam1, lam2, sig1, sig2):
    return r"\mathcal{I}^{%s%s}_{%s%s;f}" % (hname(lam1), hname(lam2), hname(sig1), hname(sig2))

def J_plain(coeff, r, m, n):
    """plain-text (coeff_str * J_r(m,n)); parenthesise the coeff if it is a sum"""
    c = to_plain(coeff)
    if sp.expand(coeff).is_Add:
        c = "( " + c + " )"
    return "%s * J_%d(%d,%d)" % (c, r, m, n)

def J_latex(coeff, r, m, n):
    """latex  coeff\,\mathcal{J}_{r}\left(m,n\right); parenthesise coeff if a sum"""
    c = to_latex(coeff)
    if sp.expand(coeff).is_Add:
        c = r"\left( %s \right)" % c
    return r"%s\,\mathcal{J}_{%d}\left(%d,%d\right)" % (c, r, m, n)

print("Mapping phi-averaged factors onto loop seeds J_r(m,n):")
all_mapped = {}                                  # (lam1,lam2,sg1,sg2) -> list of (coeff,r,m,n)
for lam1 in ext:
    for lam2 in ext:
        for sg1 in helicity:
            for sg2 in helicity:
                all_mapped[(lam1, lam2, sg1, sg2)] = map_to_J(all_avg[(lam1, lam2, sg1, sg2)])
                print("  mapped  lam1=%s, lam2=%s, sigma1=%s, sigma2=%s  (%d terms)"
                      % (hname(lam1), hname(lam2), hname(sg1), hname(sg2),
                         len(all_mapped[(lam1, lam2, sg1, sg2)])), flush=True)

# ---- APPEND J-seed expansions to box_loop_E.txt ----
with open('box_loop_E.txt', 'a') as f:
    f.write("\n\n" + "="*90 + "\n")
    f.write("MAPPING TO LOOP SEEDS  J_r(m,n)\n")
    f.write("Each phi-averaged  \\bar{\\mathcal{E}}^{lambda1 lambda2}_{sigma1 sigma2}  is expressed\n")
    f.write("as a linear combination of loop seeds\n")
    f.write("    J_r(m,n) = f_r(x - m/2 , x + n/2),   x = -i f \\tilde\\nu,\n")
    f.write("with f_0, f_c=f_1, f_cc=f_2, f_3, f_4 defined in the .tex.\n")
    f.write("sin(theta) eliminated via sin^2 = 1 - cos^2;  z = cos(theta).\n")
    f.write("J_3, J_4 are kept explicitly (NOT reduced to J_0).\n")
    f.write("This equals the dimensionless coefficient \\mathcal{I}^{lambda1 lambda2}_{sigma1 sigma2;f}.\n")
    f.write("="*90 + "\n\n")

    for lam1 in ext:
        for lam2 in ext:
            f.write("\n>>> PLAIN FORM (J-seed expansion)  lambda1=%s lambda2=%s <<<\n"
                    % (hname(lam1), hname(lam2)))
            for sg1 in helicity:
                for sg2 in helicity:
                    f.write(Ilabel(lam1, lam2, sg1, sg2) + " =\n")
                    f.write(" + ".join(J_plain(co, r, m, n)
                               for (co, r, m, n) in all_mapped[(lam1, lam2, sg1, sg2)]) or "0")
                    f.write("\n\n")

            f.write("\n>>> LATEX FORM (J-seed expansion)  lambda1=%s lambda2=%s <<<\n"
                    % (hname(lam1), hname(lam2)))
            for sg1 in helicity:
                for sg2 in helicity:
                    f.write(Ilabel(lam1, lam2, sg1, sg2) + " =\n")
                    f.write(" + ".join(J_latex(co, r, m, n)
                               for (co, r, m, n) in all_mapped[(lam1, lam2, sg1, sg2)]) or "0")
                    f.write("\n\n")

print("\nMapping done. J-seed expansions appended to box_loop_E.txt")

# ======================================================================
#  SUM OVER sigma1,sigma2 WITH BETA WEIGHTS  (Scheme A)  ->  box_loop_B.txt
# ======================================================================
#  S(lam1 lam2) = sum_{sg1,sg2 in {L,+,-}}  beta_{sg1} beta_{sg2} I^{lam1 lam2}_{sg1,sg2;f}
#  in the J-seed basis (Scheme A: sum the already-mapped terms).
#
#  The helicity expansion is all-plus (no sign weights).  Each soft
#  longitudinal (L) internal line carries one beta factor; transverse
#  (+/-) lines carry none.  The beta symbol is kept explicit (no nu
#  dependence), exactly as the J_r seeds are kept symbolic.
#
#  NO sin^2+cos^2=1 reduction is applied to the coefficients.
# ======================================================================

def beta_weight(sig):
    """One beta factor for a soft L line; 1 for a soft +/- line."""
    return beta if sig == 'L' else sp.Integer(1)

def Slabel(lam1, lam2):
    return r"S(%s%s)" % (hname(lam1), hname(lam2))

with open('box_loop_B.txt', 'w') as f:
    f.write("="*90 + "\n")
    f.write("S(lam1 lam2) = sum_{sg1,sg2 in {L,+,-}}  beta_{sg1} beta_{sg2} I^{lam1 lam2}_{sg1,sg2;f}\n")
    f.write("Each soft longitudinal (L) line carries one beta factor;\n")
    f.write("transverse (+/-) lines carry none (beta_L = beta, beta_+/- = 1).\n")
    f.write("beta is kept symbolic (no nu dependence).  Sum done at the\n")
    f.write("J-seed level (Scheme A).  NO sin^2+cos^2=1 reduction applied.\n")
    f.write("S  =  sum_{(r,m,n)} C_{rmn} * J_r(m,n)\n")
    f.write("="*90 + "\n\n")

    for lam1 in ext:
        for lam2 in ext:
            combined = {}
            for sg1 in helicity:
                for sg2 in helicity:
                    w = beta_weight(sg1) * beta_weight(sg2)
                    for (coeff, r, m, n) in all_mapped[(lam1, lam2, sg1, sg2)]:
                        key = (r, m, n)
                        combined[key] = combined.get(key, sp.Integer(0)) + w * coeff

            order = sorted(combined.keys(), key=lambda k: (k[2], k[1], k[0]))

            f.write("\n" + "="*70 + "\n")
            f.write(Slabel(lam1, lam2) + "\n")
            f.write("="*70 + "\n\n")

            f.write(">>> LATEX FORM <<<\n")
            terms_l = [J_latex(combined[k], k[0], k[1], k[2])
                       for k in order if combined[k] != 0]
            f.write(Slabel(lam1, lam2) + r" = \left[ " + " + ".join(terms_l) + r" \right]")
            f.write("\n\n")

            f.write(">>> PLAIN FORM <<<\n")
            terms_p = [J_plain(combined[k], k[0], k[1], k[2])
                       for k in order if combined[k] != 0]
            f.write(Slabel(lam1, lam2) + " = ( " + " + ".join(terms_p) + " )")
            f.write("\n\n")

print("Scheme-A beta-weighted sums written to box_loop_B.txt")

# ======================================================================
#  PARTIAL SUMS OF B  (equal +/- hard external-mode contributions)
# ======================================================================
#  Since the + and - transverse hard external modes give equal left/right
#  subdiagram contributions, we combine:
#     P1  =  S(L+) + S(L-)
#     P2  =  S(+L) + S(-L)
#     P3  =  S(++) + S(--) + S(+-) + S(-+)
#  The sums are done at the J-seed level and appended to box_loop_B.txt.
# ======================================================================

def combined_S(lam1, lam2):
    """J-seed coefficient dictionary  {(r,m,n): coeff}  for S(lam1 lam2)."""
    combined = {}
    for sg1 in helicity:
        for sg2 in helicity:
            w = beta_weight(sg1) * beta_weight(sg2)
            for (coeff, r, m, n) in all_mapped[(lam1, lam2, sg1, sg2)]:
                key = (r, m, n)
                combined[key] = combined.get(key, sp.Integer(0)) + w * coeff
    return combined

def add_S(pairs):
    """Sum the J-seed coefficient dictionaries of several (lam1, lam2) pairs."""
    total = {}
    for lam1, lam2 in pairs:
        for key, val in combined_S(lam1, lam2).items():
            total[key] = total.get(key, sp.Integer(0)) + val
    return total

partial_sums = [
    ("S(L+) + S(L-)",                 [('L', 1), ('L', -1)]),
    ("S(+L) + S(-L)",                 [(1, 'L'), (-1, 'L')]),
    ("S(++) + S(--) + S(+-) + S(-+)", [(1, 1), (-1, -1), (1, -1), (-1, 1)]),
]

with open('box_loop_B.txt', 'a') as f:
    f.write("\n\n" + "="*90 + "\n")
    f.write("PARTIAL SUMS OF B  (equal +/- hard external modes combined)\n")
    f.write("P1 = S(L+) + S(L-)\n")
    f.write("P2 = S(+L) + S(-L)\n")
    f.write("P3 = S(++) + S(--) + S(+-) + S(-+)\n")
    f.write("Summed at the J-seed level (Scheme A).  beta kept symbolic.\n")
    f.write("="*90 + "\n\n")

    for label, pairs in partial_sums:
        total = add_S(pairs)
        order = sorted(total.keys(), key=lambda k: (k[2], k[1], k[0]))

        f.write("\n" + "="*70 + "\n")
        f.write(label + "\n")
        f.write("="*70 + "\n\n")

        f.write(">>> LATEX FORM <<<\n")
        terms_l = [J_latex(total[k], k[0], k[1], k[2])
                   for k in order if total[k] != 0]
        f.write(label + r" = \left[ " + " + ".join(terms_l) + r" \right]")
        f.write("\n\n")

        f.write(">>> PLAIN FORM <<<\n")
        terms_p = [J_plain(total[k], k[0], k[1], k[2])
                   for k in order if total[k] != 0]
        f.write(label + " = ( " + " + ".join(terms_p) + " )")
        f.write("\n\n")

print("Partial sums of B written to box_loop_B.txt")
