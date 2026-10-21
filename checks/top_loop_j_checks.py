"""Recompute the Appendix B.2 soft loop through the manuscript's J moments.

Checks the azimuthal polynomials against the original vector numerator,
reduces the J finite sums symbolically, and compares independent Gamma
evaluations with the previous closed coefficients. No hard/time integral
or full box amplitude is tested.
"""
from functools import lru_cache
from math import factorial
from pathlib import Path
import re

import mpmath as mp
import sympy as s
from sympy.parsing.sympy_parser import (
    implicit_multiplication_application, parse_expr, standard_transformations,
)

q, z, x = s.symbols("q z x")
half = s.Rational(1, 2)
w = q*z + half
r2 = q**2*(1-z**2)
ell2 = r2+w*w
polys = {
    "M": ell2**2+ell2/4-w*w/2,
    "P": (ell2+s.Rational(1, 4))*r2/2,
    "Q": (ell2+s.Rational(1, 4))*(w*w-r2/2)-w*w/2,
    "H0": r2*(r2-w*w+s.Rational(1, 4))/4,
}
polys["H2"] = w*w*r2-polys["H0"]
polys["H3"] = w**4-3*w*w*r2-w*w/4+polys["H0"]
polys["HB"] = (4*w*w*r2-r2/4
               -(ell2-s.Rational(1, 4))*(w*w+r2/2))/4
polys["HD"] = r2*(s.Rational(1, 4)-w*w)/4
polys = {k: s.expand(v) for k, v in polys.items()}


def check_azimuth():
    # Form the original bilinear numerator before angular averaging.
    a, b, c1, s1, c3, s3, cp, sp = s.symbols("a b c1 s1 c3 s3 cp sp")
    ell = s.Matrix([a, b, w])
    n1, n3 = s.Matrix([s1, 0, c1]), s.Matrix([s3*cp, s3*sp, c3])
    axis = s.Matrix([0, 0, 1])

    def numerator(left, right):
        return (2*left.dot(ell)*right.dot(ell)
                -left.dot(axis)*right.dot(axis)/2
                -left.dot(right)*(ell.dot(ell)-s.Rational(1, 4)))

    def average(expr):
        result = 0
        for (i, j), coeff in s.Poly(s.expand(expr), a, b).terms():
            if i % 2 or j % 2:
                continue
            k = (i+j)//2
            moment = s.Rational(factorial(i)*factorial(j),
                                4**k*factorial(i//2)*factorial(j//2)*factorial(k))
            result += coeff*moment*r2**k
        return s.expand(result).subs(sp**2, 1-cp**2).subs(
            s1**4, (1-c1**2)**2).subs(s1**2, 1-c1**2).subs(
            s3**4, (1-c3**2)**2).subs(s3**2, 1-c3**2).expand()

    angular = (polys["H0"]+polys["H2"]*(c1**2+c3**2)
               +polys["H3"]*c1**2*c3**2
               +polys["HB"]*4*c1*c3*s1*s3*cp
               +polys["HD"]*(1-c1**2)*(1-c3**2)*(2*cp**2-1))
    cases = [
        (ell, ell, polys["M"]),
        (ell, n3*n3.dot(ell), polys["P"]+c3**2*polys["Q"]),
        (n1*n1.dot(ell), ell, polys["P"]+c1**2*polys["Q"]),
        (n1*n1.dot(ell), n3*n3.dot(ell), angular),
    ]
    for left, right, target in cases:
        assert s.expand(average(numerator(left, right))-target) == 0
    print("PASS four original vector numerators vs azimuthal polynomials")


@lru_cache(None)
def shifted_f0_ratio(m, n):
    # Gamma recurrence, relative to f0(x+1/2,x+1/2).
    a0 = x+half
    return (s.rf(s.Rational(3, 2)-a0, m)*s.rf(s.Rational(3, 2)-a0, n)
            *s.rf(1-a0, m)*s.rf(1-a0, n)
            /s.rf(s.Rational(5, 2)-2*a0, m+n)/s.rf(3-2*a0, m+n))


@lru_cache(None)
def j_ratio(power, order):
    # q^(power-1) K^-1 z^order -> J_order(power-1,1;x).
    assert (power-order) % 2 == 0
    result = 0
    for u in range(order+1):
        for v in range(order-u+1):
            h = order-u-v
            coeff = s.Rational(factorial(order)*(-1)**(v+h),
                               2**order*factorial(u)*factorial(v)*factorial(h))
            m, n = (power-order)//2+v, u
            assert m >= 0 and n >= 0
            result += coeff*shifted_f0_ratio(m, n)
    return s.factor(result)


def reduce_poly(poly):
    return s.factor(sum(coeff*j_ratio(power, order)
                        for (power, order), coeff in s.Poly(poly, q, z).terms()))


def f0(a, b):
    return (mp.gamma(a+b-mp.mpf("1.5"))*mp.gamma(mp.mpf("1.5")-a)
            *mp.gamma(mp.mpf("1.5")-b)
            /mp.gamma(3-a-b)/mp.gamma(a)/mp.gamma(b))


def j_numeric(power, order, xv):
    # Evaluate every shifted Gamma function directly, without recurrence.
    a, b = xv-mp.mpf(power-1)/2, xv+mp.mpf(".5")
    return mp.fsum(
        mp.mpf(factorial(order)*(-1)**(v+order-u-v))
        /(2**order*factorial(u)*factorial(v)*factorial(order-u-v))
        *f0(a+mp.mpf(order)/2-v, b-u)
        for u in range(order+1) for v in range(order-u+1)
    )


def check_reduction():
    rho = x+half  # Used only to compare against the pre-rewrite manuscript.
    den = 32*(rho-3)*(rho-2)*(4*rho-7)*(4*rho-5)
    old = {
        "P": (rho-1)*(2*rho-3)/(8*(rho-2)*(4*rho-7)*(4*rho-5)),
        "Q": 0,
        "H0": (rho-1)*(2*rho-3)*(2*rho-5)/den,
        "H2": -(rho-1)*(2*rho-3)**2/den,
        "HB": -2*(rho-1)*(rho-2)*(2*rho-3)/den,
        "HD": (rho-2)*(2*rho-5)*(2*rho-3)/den,
    }
    old["M"], old["H3"] = 3*old["P"], -3*old["H2"]
    reduced = {}
    for name, poly in polys.items():
        reduced[name] = reduce_poly(poly)
        assert s.factor(reduced[name]-old[name]) == 0, name
        print(f"{name}: polynomial = {poly}")
        print(f"{name}: reduced / f0 = {reduced[name]}")
    print("PASS all eight J combinations vs previous closed coefficients")
    mp.mp.dps = 70
    maxerr = mp.mpf(0)
    for mu in (mp.mpf(".2"), mp.mpf(".73"), mp.mpf("2"), mp.mpf("5")):
        xv = -1j*mu
        for name, poly in polys.items():
            actual = mp.fsum(
                mp.mpf(str(coeff))*j_numeric(power, order, xv)
                for (power, order), coeff in s.Poly(poly, q, z).terms())
            expected = s.lambdify(x, reduced[name], "mpmath")(xv)*f0(xv+half, xv+half)
            err = abs(actual-expected)/max(abs(f0(xv+half, xv+half)), abs(expected))
            maxerr = max(maxerr, err)
            assert err < mp.mpf("1e-60"), (mu, name, err)
    print(f"PASS independent Gamma sums at four masses: scaled residual {mp.nstr(maxerr, 4)}")


def check_nu_form():
    nu = s.symbols("nu", real=True)
    common = -(1+s.I*nu)/(4*(5+2*s.I*nu)*(3+2*s.I*nu)
                         *(5+4*s.I*nu)*(3+4*s.I*nu))
    result = {
        "P": -(1+s.I*nu)*(1+2*s.I*nu)/(
            4*(3+2*s.I*nu)*(5+4*s.I*nu)*(3+4*s.I*nu)),
        "H0": common*(2+s.I*nu)*(1+2*s.I*nu),
        "H2": -common*(1+s.I*nu)*(1+2*s.I*nu),
        "HB": -common*(1+2*s.I*nu)*(3+2*s.I*nu)/2,
        "HD": common*(2+s.I*nu)*(3+2*s.I*nu),
    }
    for name, expression in result.items():
        assert s.cancel(expression-reduce_poly(polys[name]).subs(x, -s.I*nu)) == 0
    print("PASS final nu-only rational factors multiplying the common Gamma function")


def check_manuscript_j_sums():
    """Check the detailed note and the compact SMS expressions."""
    root = Path(__file__).resolve().parents[1]
    tex = (root/"Top Quark Signal.tex").read_text()
    blocks = []
    for label in ("top-loop-J-isotropic", "top-loop-J-angular"):
        end = tex.index(r"\label{eq:"+label+"}")
        start = tex.rfind(r"\begin{align}", 0, end)
        blocks.append(tex[start:end])
    name_map = {"M": "M", "P": "P", "Q": "Q", "0": "H0",
                "2": "H2", "3": "H3", "B": "HB", "D": "HD"}
    variables = {}

    def j_symbol(match):
        order, power, second = map(int, match.groups())
        assert second == 1
        name = f"j{order}m{power}n{second}"
        variables[name] = s.Symbol(name)
        return f"({name})"

    found = set()
    for block in blocks:
        entries = list(re.finditer(r"\\mathcal U_([MPQ023BD])=\{\}&", block))
        for i, match in enumerate(entries):
            stop = entries[i+1].start() if i+1 < len(entries) else len(block)
            rhs = block[match.end():stop]
            rhs = rhs.replace(r"\nonumber\\", "").replace("&", "").strip().rstrip(",.")
            rhs = re.sub(r"\\mathcal J_(\d)\((\d+),(\d+)\)", j_symbol, rhs)
            rhs = re.sub(
                r"\\frac(\{[0-9]+\}|[0-9])(\{[0-9]+\}|[0-9])",
                lambda m: f"({m[1].strip('{}')}/{m[2].strip('{}')})", rhs)
            rhs = rhs.replace("[", "(").replace("]", ")")
            rhs = " ".join(rhs.split())
            parsed = parse_expr(rhs, local_dict=variables, transformations=(
                standard_transformations+(implicit_multiplication_application,)))
            name = name_map[match[1]]
            expected = sum(
                coeff*s.Symbol(f"j{order}m{power-1}n1")
                for (power, order), coeff in s.Poly(polys[name], q, z).terms())
            assert s.expand(parsed-expected) == 0, (name, parsed, expected)
            found.add(name)
    assert found == set(polys)
    print("PASS all eight J sums in Top Quark Signal")

    sms = (root/"SMS in CHC.tex").read_text()
    end = sms.index(r"\label{eq:top-loop-reduction}")
    start = sms.rfind(r"\begin{align}", 0, end)
    block = sms[start:end]
    lhs12 = r"(4\pi)^{3/2}\mathcal I^{12}_{\mathrm{top}}"
    lhs22 = r"(4\pi)^{3/2}\mathcal I^{22}_{\mathrm{top}}"
    first, second = block.index(lhs12), block.index(lhs22)
    terms = {
        "P": block[first+len(lhs12):second],
        "H": block[second+len(lhs22):],
    }
    angular = {
        r"(c_1^2+c_3^2-3c_1^2c_3^2)": "angularA",
        r"\sin(2\theta_1)\sin(2\theta_3)\cos\phi_3": "angularB",
        r"\sin^2\theta_1\sin^2\theta_3\cos(2\phi_3)": "angularD",
    }
    for symbol in angular.values():
        variables[symbol] = s.Symbol(symbol)

    def formal_j(poly):
        return sum(coeff*s.Symbol(f"j{order}m{power-1}n1")
                   for (power, order), coeff in s.Poly(poly, q, z).terms())

    expected = {
        "P": formal_j(polys["P"]),
        "H": formal_j(polys["H0"])
        +variables["angularA"]*formal_j(polys["H2"])
        +variables["angularB"]*formal_j(polys["HB"])
        +variables["angularD"]*formal_j(polys["HD"]),
    }
    for name, rhs in terms.items():
        rhs = rhs.replace("&=", "").replace(r"\nonumber\\", "").replace("&", "")
        rhs = rhs.replace(r"\quad", "").replace(r"\Big", "")
        rhs = re.sub(r"\\hspace\{[^}]+\}", "", rhs)
        rhs = " ".join(rhs.split())
        for latex, symbol in angular.items():
            rhs = rhs.replace(latex, f"({symbol})")
        rhs = re.sub(r"\\mathcal J_(\d)\((\d+),(\d+)\)", j_symbol, rhs)
        rhs = re.sub(
            r"\\frac(\{[0-9]+\}|[0-9])(\{[0-9]+\}|[0-9])",
            lambda m: f"({m[1].strip('{}')}/{m[2].strip('{}')})", rhs)
        rhs = rhs.replace("[", "(").replace("]", ")").strip().rstrip(",.")
        parsed = parse_expr(rhs, local_dict=variables, transformations=(
            standard_transformations+(implicit_multiplication_application,)))
        assert s.expand(parsed-expected[name]) == 0, name
    print("PASS compact I(J) expressions in SMS against the detailed polynomials")


if __name__ == "__main__":
    check_azimuth()
    check_reduction()
    check_nu_form()
    check_manuscript_j_sums()
