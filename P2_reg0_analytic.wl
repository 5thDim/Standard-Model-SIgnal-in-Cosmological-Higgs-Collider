(* Run this after the notebook definitions of P2L and P2PM have been evaluated. *)

ClearAll[
  reg0, nu0, hPM, hL, P2PMReg0, P2LReg0,
  p2PMExpression, p2LExpression, p2PMHyper, p2LHyper,
  p2PMDummy, p2LDummy, p2PMFiniteCoefficient, p2LFiniteCoefficient,
  P2PMPoleCoefficients, P2LPoleCoefficients
];

p2PMExpression = Expand[P2PM[nu0, reg0]];
p2LExpression = Expand[P2L[nu0, reg0]];

(* Store every hypergeometric object before replacing it by an inert one-variable
   placeholder.  This stops SeriesCoefficient from trying to prove general
   parameter identities for 3F2(1). *)
p2PMHyper = Cases[
   p2PMExpression,
   _HypergeometricPFQ | _HypergeometricPFQRegularized,
   Infinity
];

p2LHyper = Cases[
   p2LExpression,
   _HypergeometricPFQ | _HypergeometricPFQRegularized,
   Infinity
];

Module[{j = 0},
  p2PMDummy = p2PMExpression /.
    (_HypergeometricPFQ | _HypergeometricPFQRegularized) :>
      (++j; hPM[j][reg0])
];

Module[{j = 0},
  p2LDummy = p2LExpression /.
    (_HypergeometricPFQ | _HypergeometricPFQRegularized) :>
      (++j; hL[j][reg0])
];

(* Exact Laurent coefficients.  No numerical fitting is used. *)
p2PMFiniteCoefficient = SeriesCoefficient[p2PMDummy, {reg0, 0, 0}];
p2LFiniteCoefficient = SeriesCoefficient[p2LDummy, {reg0, 0, 0}];

P2PMPoleCoefficients =
  Table[SeriesCoefficient[p2PMDummy, {reg0, 0, k}], {k, -3, -1}];

P2LPoleCoefficients =
  Table[SeriesCoefficient[p2LDummy, {reg0, 0, k}], {k, -3, -1}];

(* Restore the exact hypergeometric functions.  Derivative[hPM[i]][0] and
   Derivative[hL[i]][0] therefore mean exact parameter derivatives along the
   regulator path appearing in the original P2PM/P2L definitions. *)
hPM[i_Integer][z_] := p2PMHyper[[i]] /. reg0 -> z;
hL[i_Integer][z_] := p2LHyper[[i]] /. reg0 -> z;

P2PMReg0[nu_] := p2PMFiniteCoefficient /. nu0 -> nu;
P2LReg0[nu_] := p2LFiniteCoefficient /. nu0 -> nu;

(* Useful checks:

   N[P2PMReg0[3/5], 20]
   N[P2PM[3/5, 10^-8], 20]

   The first pair agrees in the inspected notebook:
   2.8700491344849296112 + 4.0086036532547537370 I
   2.8700491208826029669 + 4.0086035869218648045 I

   P2LReg0 is exact but slower to numericalize because it contains 18
   hypergeometric parameter-derivative objects.
*)
