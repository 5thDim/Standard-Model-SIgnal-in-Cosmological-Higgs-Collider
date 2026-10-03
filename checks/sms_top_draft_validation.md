# SMS top-quark chapter validation

Date: 2026-10-04

## Changes

- Existing gauge-boson text: coefficient species labels only (`g`), with the
  original equation content and order retained.
- Existing appendices: the four sections became A.1--A.4 in their original
  order; their subsections became subsubsections. No material was moved between
  the gauge-boson main text and appendix.
- Added Section 3 (five subsections) and Appendix B (four subsections).
- The Whittaker seed is `mathcal K`; two hard moments express the time
  coefficients directly through it, eliminating the note's intermediate
  J/S/R moment names. Long seed and coefficient formulas are in Appendix B.
- Added the two fermion references already present in the local note.
- Existing plot assets were retained. Equations and captions use the species
  labels. The new top-quark cut diagram is drawn in LaTeX.
- `Top Quark Signal.tex` and its PDF retain their pre-edit SHA-256 hashes.

## Preservation and build

The first, gauge-only edit was reversed mechanically and recovered the original
SMS TeX byte for byte. After adding the new material, removing Section 3,
Appendix B and the new preamble definitions recovered that gauge-only version
up to whitespace at the insertion boundaries. All 87 original labels remain
in their original order. There are 67 new labels, with no duplicates or
unresolved references. The gauge-boson main equations remain (1)--(37).

Build command:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -outdir=/tmp/sms-top-draft/build 'SMS in CHC.tex'
```

The resulting PDF has 31 pages. No undefined references/citations, duplicate
labels, overfull boxes or underfull boxes were reported. The contents page,
gauge-boson box, Appendix A opening, top-quark setup, cancellation, signal,
seed and explicit-coefficient pages were inspected in the rendered PDF.
Section 3 is on pages 15--19; Appendix B starts on page 25.

## Formula checks

1. `PYTHONDONTWRITEBYTECODE=1 python3 -u checks/direct_trace_checks.py` passed.
   The independent first-order hard blocks were checked at masses 0.2, 0.73
   and 2 for every SK pair. The four separate traces agree with the scalar
   expression for all 16 SK assignments, including unequal exact Higgs factors
   (maximum absolute residual below 9e-19). This also checks that the
   cancellation is between actual external orderings.
2. The direct K-only formulas for H0, H1 and the radial derivative were
   compared symbolically with the original moment formulas, treating the
   Whittaker seeds as arbitrary functions. All four SK pairs gave exact zero
   differences. The derivative comparison includes all explicit powers of r
   and both shifted-seed prefactors.
3. The soft-loop coefficient matrix, including `8 c_+^2 16^(2 i nu)`, was
   compared with the note's final angular-coefficient evaluator at masses
   0.2 and 0.73 and three angular configurations. Maximum relative residual:
   6.1e-26 at 25-digit precision.
4. The two mixed soft branches, four orderings, Pauli trace and soft mode
   normalization reproduce the prefactor of the nonoscillatory signal.

The attempted direct finite-regulator numerical seed comparison was stopped
because of slow special-function evaluation; the exact symbolic identity
check above replaces that notation-reduction check. No result from the
interrupted run is counted as a pass.

These checks validate the manuscript transfer, condensed seed algebra and
factorized leading signal. They are not a new integration of the full
unexpanded box or a validation of analytic contact-term renormalization.

Original snapshots, notation map, build log, validation script and outputs:
`/tmp/sms-top-draft/`.
