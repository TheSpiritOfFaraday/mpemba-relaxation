# Independent finite-time certificate audit

Checked and rerun 11 October 2026: `python3 research/finite_time_certificate.py`.

Verdict: all exact checks passed. The certificate proves at least one total-variation crossing in (0.10,0.11) and at least one KL crossing in (0.14,0.16). It does not certify uniqueness, first-crossing status, or a more precise crossing time.

The code was reviewed line by line. All acceptance computations use `fractions.Fraction`; `float` is used only to print diagnostic approximations. Rational rounding floors lower endpoints and ceilings upper endpoints at a 50-decimal grid. Gibbs-root input intervals are revalidated by exact powers. The uniformization matrix is checked entrywise nonnegative with exact unit column sums.

For x=253t, the exponential-series remainder after degree N is bounded by the first omitted term divided by 1-x/(N+2). Subsequent term ratios are no larger than x/(N+2), and that denominator is positive at every certified time. In the omitted matrix series, each actual probability coordinate is at most one. Multiplication by the reciprocal exponential enclosure therefore encloses the exact evolved state.

The logarithm uses log(a)=2 sum_{k>=0} z^(2k+1)/(2k+1), z=(a-1)/(a+1), with absolute remainder bounded by 2|z|^(2N+1)/[(2N+1)(1-z²)]. Positivity and |z|<1 are asserted. Monotonicity of log converts endpoint enclosures into interval enclosures. Absolute-value bounds, KL products, and sums remain valid despite correlations between interval variables; losing those correlations enlarges the enclosure rather than invalidating it.

The exact endpoint signs in the JSON output, together with continuity of the Markov evolution and the distances, prove the stated crossing brackets. Plot root estimates may be reported separately as high-precision numerical values.
