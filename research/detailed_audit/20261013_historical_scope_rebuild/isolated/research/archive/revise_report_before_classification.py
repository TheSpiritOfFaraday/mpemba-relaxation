"""Idempotently build revised report from archived initial report and figure data."""
from pathlib import Path
import csv,json
ROOT=Path(__file__).resolve().parents[1]
s=(ROOT/'research/archive/research_report_before_physics_revision.tex').read_text()
s=s.replace('\\usepackage{amsmath,amssymb,amsthm,booktabs,longtable,array,hyperref}','\\usepackage{amsmath,amssymb,amsthm,booktabs,longtable,array,hyperref}\n\\usepackage{pgfplots}\n\\usepgfplotslibrary{groupplots}\n\\pgfplotsset{compat=1.18}')
s=s.replace('We also prove a corrected sufficient condition for absence using the generator acting on adjacent differences, and derive a finite-temperature theorem for the entire all-$\\kappa\\leq1$ subfamily, for any number of states.','We connect the triplet inequalities to established stochastic monotonicity, obtaining corrected sufficient conditions for absence. The supporting order argument is presented as an application of existing theory, with explicit comparison to Ohga, Hayakawa and Ito. Overlap plots and a certified finite-time distance crossing show the physical effect, whose amplitude advantage is about 0.1102\\%.')
s=s.replace('The main result is a mathematical counterexample, not an inference from an unsuccessful search or an approximate spectral plot.','The proposed triplet test can miss a genuine four-state relaxation reversal. In the example, the hotter preparation begins farther from equilibrium, yet later becomes slightly closer. The reversal is mathematically rigorous and physically small.')
s=s.replace('The decisive conclusion here is eventual relaxation ordering, not a prescribed finite crossing time.','The asymptotic conclusion is supplemented below by certified finite-time crossing brackets.')
s=s.replace('For definitions we use',r'''For supporting no-Mpemba results we additionally checked \href{https://arxiv.org/html/2410.06623v2}{Ohga, Hayakawa and Ito, arXiv:2410.06623v2}, version-history date 19 October 2024, including Result 2 and SI B--D, on 11 October 2026. Their theorem is a prior application of monotone coupling to absence of mean-energy crossings. Section~\ref{sec:priororder} compares the assumptions and observables precisely. We do not claim a new general monotonicity principle.

For definitions we use''')
# Move full algebraic expansions out of the central proof.
begin=s.index('Define four polynomials\n')
end=s.index('The nonvanishing denominator is certified',begin)
old=s[begin:end]
polys=old[old.index('\\begin{align*}'):old.index('Their common nonzero factor')]
s=s[:begin]+r'''Let $P_i(x)$ be the explicit cofactor polynomials listed in \hyperref[app:cofactors]{Appendix A}. Their common factor corresponds to the first column of $\operatorname{adj}(M^\top-xI)$. Put $w_i=P_i(r)/P_1(r)$. '''+s[end:]
begin=s.index('\\subsection{An exact rational slow-mode companion}')
end=s.index('\\section{Mathematical explanation and corrected theorem}',begin)
companion=s[begin:end].replace('\\subsection{An exact rational slow-mode companion}',r'\section*{Appendix B: A rational slow-mode companion}')
s=s[:begin]+s[end:]
# Figures use inlined data, keeping the LaTeX source standalone for Codex preview.
def coords(rows,x,y,stride=1):return ' '.join(f'({r[x]:.10g},{r[y]:.12g})' for r in (rows[::stride]+[rows[-1]]))
over=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((ROOT/'output/figures/overlap_data.csv').open())]
rel=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((ROOT/'output/figures/relaxation_data.csv').open())]
for r in rel:r['ratio']=100*(r['KL_hot']/r['KL_cold']-1)
fig1=r'''\begin{figure}[htbp]
\centering
\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.6cm},width=0.43\textwidth,height=5.3cm,tick label style={font=\scriptsize},label style={font=\small},title style={font=\small},xlabel={Initial inverse temperature $\beta$},scaled ticks=false]
\nextgroupplot[title={A: Full cooling interval},xmin=0,xmax=1,ymin=0,ymax=.1,ylabel={$|b_2(\beta)|$}]
\addplot[black,thick] coordinates {FULL};
\addplot[only marks,red,mark=*] coordinates {(.001,.0937454279941)};
\addplot[only marks,blue,mark=*] coordinates {(.01,.0938488356221)};
\nextgroupplot[title={B: High-temperature enlargement},xmin=0,xmax=.07,xtick={0,.02,.04,.06},xticklabel style={/pgf/number format/fixed,/pgf/number format/precision=2},ymin=-.55,ymax=.22,ylabel={Amplitude increase (\%)},legend style={font=\scriptsize,draw=none,at={(.02,.02)},anchor=south west}]
\addplot[black,thick,forget plot] coordinates {ZOOM};
\addplot[gray,forget plot] coordinates {(0,0)(.07,0)};
\addplot[only marks,red,mark=*] coordinates {(.001,0)};\addlegendentry{$T_h=1000$}
\addplot[only marks,blue,mark=*] coordinates {(.01,.1103068468)};\addlegendentry{$T_c=100$}
\end{groupplot}
\end{tikzpicture}
\caption{Thermal overlap of the exact primary generator, with $w_1=1$ and $T_b=1$. The broad curve conceals a small high-temperature maximum; the enlargement displays it. The blue preparation has a larger slow amplitude than the hotter red preparation. The enlarged ordinate uses the hot amplitude as its reference; the quoted 0.1102\% reduction uses the cold amplitude. Numerical curves use 80-digit arithmetic; the two selected overlaps are independently certified exactly.}
\label{fig:overlap}
\end{figure}
'''.replace('FULL',coords(over[:350],'beta','absolute_overlap',3)).replace('ZOOM',coords(over[350:],'beta','enhancement_percent',2))
# include last endpoint omitted by sampling to preserve full curve span
fig2=r'''\begin{figure}[htbp]
\centering
\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.6cm},width=0.43\textwidth,height=5.4cm,tick label style={font=\scriptsize},label style={font=\small},title style={font=\small},xlabel={Time $t$ (units $a_{12}^{-1}$)},xmin=0,xmax=.5]
\nextgroupplot[title={A: Relaxation distance},ymode=log,ymin=1e-11,ymax=2,ytick={1,1e-2,1e-4,1e-6,1e-8,1e-10},ylabel={$D_{\rm KL}(p(t)\Vert\pi_b)$},legend style={font=\scriptsize,draw=none,at={(.98,.98)},anchor=north east}]
\addplot[red,thick] coordinates {HOT};\addlegendentry{Hot preparation}
\addplot[blue,dashed,thick] coordinates {COLD};\addlegendentry{Cool preparation}
\nextgroupplot[title={B: Relative ordering},ylabel={$100(D_h/D_c-1)$ (\%)},ymin=-.6,ymax=6.2]
\addplot[black,thick] coordinates {RATIO};
\addplot[gray] coordinates {(0,0)(.5,0)};
\addplot[gray,dotted] coordinates {(0,-.2202492067)(.5,-.2202492067)};
\addplot[only marks,violet,mark=*] coordinates {(.1496723312,0)};
\node[font=\scriptsize,anchor=west] at (axis cs:.21,1) {Crossing near $0.14967$};
\draw[->,gray] (axis cs:.22,.9) -- (axis cs:.1496723312,.04);
\end{groupplot}
\end{tikzpicture}
\caption{Actual evolution $p(t)=e^{tM}\pi(\beta)$ for the same two preparations. Relative entropy decreases for each, and the hotter state's distance becomes smaller. The ratio panel resolves the crossing that is difficult to see on the logarithmic plot. The horizontal dotted line is the asymptotic ratio difference, $-0.2202492\%$. Rational uniformization and logarithm bounds certify opposite KL orderings at $t=0.14$ and $0.16$; the displayed root location is numerical, with no claim of uniqueness.}
\label{fig:crossing}
\end{figure}
'''.replace('HOT',coords(rel,'time','KL_hot',3)).replace('COLD',coords(rel,'time','KL_cold',3)).replace('RATIO',coords(rel,'time','ratio',2))
physics=r'''\subsection{Visible reversal and its size}
Figures~\ref{fig:overlap} and \ref{fig:crossing} separate two facts: the hotter initial state has a smaller slow-mode amplitude, and the full time-dependent distance actually changes ordering. These use the primary generator in Theorem~1, with the bath and generator fixed throughout. The high-temperature amplitude reaches its numerical maximum near $\beta=0.0237422$.

FIGONE

Initially $D_h(0)\simeq1.210137886$ and $D_c(0)\simeq1.155489557$, so the hotter preparation is farther away in relative entropy. The two distances cross numerically near $t=0.14967233$. Exact rational bounds prove a crossing in $(0.14,0.16)$: at the endpoints,
\begin{align*}
D_h(0.14)-D_c(0.14)&\in(5.9361291891,5.9361291893)\,10^{-8},\\
D_h(0.16)-D_c(0.16)&\in(-1.3100324358,-1.3100324356)\,10^{-8}.
\end{align*}
Total variation also reverses, with a certified crossing in $(0.10,0.11)$ and numerical location $0.10436410$. The difference between these times is expected: finite-time crossings depend on the distance. Continuity plus the certified signs proves existence; no uniqueness or earliest-crossing theorem is asserted.

FIGTWO

The relevant fractional amplitude reduction is
\[
\epsilon_b=1-\frac{|b_h|}{|b_c|}\simeq0.0011018531=0.11018531\%.
\]
For relative entropy the leading term is quadratic, so
\[
1-\lim_{t\to\infty}\frac{D_h(t)}{D_c(t)}
=1-\frac{b_h^2}{b_c^2}\simeq0.22024921\%.
\]
Both preparations have the same decay exponent $\lambda_2$. This is an ordinary effect, not the exponentially faster strong effect. For an increasingly small common distance threshold, the asymptotic time advance is
\[
\Delta t=\frac{\log(|b_c|/|b_h|)}{|\lambda_2|}
\simeq5.5488\,10^{-5}\,a_{12}^{-1}.
\]
The slow relaxation time is $1/|\lambda_2|\simeq0.0503313\,a_{12}^{-1}$, so this advance is only about $0.11025\%$ of one slow relaxation time. It is not a percentage reduction in the full cooling time, which depends on the chosen threshold. The time unit has no assigned experimental value, and $T_h/T_b=1000$, $T_c/T_b=100$. These data establish a strict mathematical failure of the triplet inference; they do not establish a practically useful cooling protocol.

'''.replace('FIGONE',fig1).replace('FIGTWO',fig2)
s=s.replace('\\subsection{A corrected sufficient condition}',physics+'\\subsection{A corrected sufficient condition}')
# precise attribution and all-time consequence after original spectral theorem
needle='A second subfamily consists of downward gap laws'
idx=s.index(needle)
end=s.index('\\section{Approaches explored}',idx)
s=s[:idx]+s[idx:end].replace('These order arguments are independently derived here; priority relative to the wider stochastic-monotonicity literature has not been established.','The following comparison locates this argument within established order theory; no general monotonicity novelty is claimed.')+s[end:]
prior=r'''\subsection{Relation to earlier no-Mpemba theorems}\label{sec:priororder}
\href{https://arxiv.org/html/2410.06623v2}{Ohga, Hayakawa and Ito (OHI), Result 2 and SI D} already use monotone coupling to exclude Mpemba behavior. Their SI D2 explicitly treats a total order. Our matrix condition is precisely a finite-chain expression of stochastic monotonicity, not a separate physical principle. The relevant distinctions are:
\begin{center}\small
\begin{tabular}{p{.17\textwidth}p{.36\textwidth}p{.36\textwidth}}
\toprule &OHI Result 2&The preceding spectral criterion\\\midrule
State order&Finite partial order; energy respects it&Total order by distinct energies\\
Dynamics&Irreducible and Gibbs-stationary; reversibility unnecessary&Positive reversible rates\\
Order condition&Comparable jumps and a monotone coupling&$K=BLH_0$ is Metzler\\
Bath range&Hypotheses throughout a temperature interval&One fixed finite positive bath\\
Spectrum&No simplicity requirement&Simple slowest nonzero eigenvalue\\
Conclusion&No cooling mean-energy crossing at any time&No asymptotic cooling/heating reversal in fixed norms or KL\\\bottomrule
\end{tabular}
\end{center}
OHI SI B extends their Result 1, the microscopic characterization, to KL divergence; it does not state Result 2 as an all-time KL theorem. Thus the different observable and time scopes must not be conflated.

Indeed $BL=KB$ gives $Be^{tL}=e^{tK}B$. The Metzler condition is equivalent to $e^{tL}$ preserving increasing functions: sufficiency follows from $e^{tK}\geq0$, and necessity by differentiation at zero on coordinate vectors of adjacent differences. This also gives a monotone coupling. Choose $r$ large enough that $I+L/r$ is stochastic and $I+K/r$ is entrywise nonnegative. Its transition rows are stochastically ordered. A common Poisson clock of rate $r$ and common-uniform inverse-CDF updates preserve the order of two copies. On a total order, energy compatibility and comparability of jumps are automatic.

For fixed bath, $e^{tL}E$ is then increasing, so
\[
\partial_\beta\{\pi(\beta)^\top e^{tL}E\}
=-\Cov_{\pi(\beta)}(e^{tL}E,E)\leq0.
\]
This gives all-time energy ordering without reversibility or spectral simplicity. It explains the direct relation to OHI while avoiding an inference from a single observable to all possible distances.

The likelihood-ratio and stochastic-order methods below also have precedent in \href{https://arxiv.org/abs/1109.6075}{Fill and Kahn, Ann. Appl. Probab. 23 (2013), 1778--1816}, particularly their comparison framework and Corollary 3.3. Their setup compares different kernels under suitable initial laws; the thermal preparations here differ instead. We provide the following consequence for clarity, without a priority claim.

\begin{corollary}[All-time divergence ordering under reversible monotonicity]
Let a finite irreducible generator be reversible with positive Gibbs equilibrium, on a total energy order, and let $K$ be Metzler. For any convex $\Phi$ finite on positive arguments, with $\Phi(1)=0$, define
\[
D_\Phi(p\Vert\pi_b)=\sum_i\pi_{bi}\Phi(p_i/\pi_{bi}).
\]
For same-side Gibbs preparations, the farther initial temperature has at least as large a divergence at every $t\geq0$. This includes KL and total variation, and requires no simple eigenvalue. For KL the inequality is strict at every finite time for distinct initial temperatures.
\end{corollary}
\begin{proof}
Write $p\succeq_{\rm st}q$ when $\sum_i p_i f_i\geq\sum_i q_i f_i$ for every increasing $f$. Gibbs laws are stochastically ordered by temperature, by the pairwise covariance identity. A monotone semigroup preserves this order. Thus for cooling $p_h(t)\succeq_{\rm st}p_c(t)\succeq_{\rm st}\pi_b$.

Reversibility implies $p_c(t)/\pi_b=e^{tL}[\pi(\beta_c)/\pi_b]$ componentwise. The initial ratio is proportional to $e^{(\beta_b-\beta_c)E_i}$ and is increasing; it stays increasing. Choose a nondecreasing subgradient $g$ of $\Phi$. Convexity and stochastic order give
\[
D_\Phi(p_h\Vert\pi_b)-D_\Phi(p_c\Vert\pi_b)
\geq\sum_i(p_{hi}-p_{ci})g(p_{ci}/\pi_{bi})\geq0.
\]
For heating, $p_c\preceq_{\rm st}p_h\preceq_{\rm st}\pi_b$ and $p_h/\pi_b$ is decreasing. Apply the same supporting-line inequality at $p_h/\pi_b$ to obtain $D_\Phi(p_c\Vert\pi_b)\geq D_\Phi(p_h\Vert\pi_b)$. KL is strictly convex, and invertibility of $e^{tM}$ keeps two distinct initial distributions distinct at finite $t$, proving strictness.
\end{proof}
This does not assert all-time ordering for every arbitrary norm. The original spectral proof supplies the every-fixed-norm conclusion asymptotically with a simple slow mode. Our counterexample violates monotonicity explicitly: $K_{21}=a_{12}-a_{13}=-27$. It therefore challenges neither OHI's theorem nor this corollary.

The defensible contribution claim is the exact counterexample to the unrestricted triplet inference, its robust family and mechanism, and the explicit triplet-rate connection to an established sufficient condition. We do not present the general order-preservation mechanism or the elementary divergence consequence as newly discovered principles. The all-$\kappa\leq1$ statement is an application; its broader publication priority remains unestablished.

'''
s=s.replace('\\section{Approaches explored}',prior+'\\section{Approaches explored}')
# Verification additions
s=s.replace('\\subsection{Executed discovery calculations}',r'''\subsection{Figure and finite-time verification}
Run \path{research/physics_plots.py} for the figures and data, and \path{research/finite_time_certificate.py} for exact finite-time bounds. Plot values use 80-digit mpmath arithmetic. Direct matrix exponentiation and reversible symmetric spectral propagation agree within $6\,10^{-81}$ at ten checks spanning $t=0$ to $0.4$. The plot CSV files and numerical metrics are saved under \path{output/figures/}. Those numerical comparisons are distinct from the exact certificates.

The finite-time certificate uses uniformization with $\nu=253$ and the exactly stochastic matrix $P=I+M/\nu$:
\[
e^{tM}p(0)=e^{-\nu t}\sum_{k=0}^{\infty}\frac{(\nu t)^k}{k!}P^kp(0).
\]
It retains terms through $k=250$. Each omitted probability component is at most one; the scalar positive-series tail is bounded by its first omitted term divided by $1-\nu t/252$. The same positive Taylor bounds give a rational enclosure of $e^{-\nu t}$. Initial Gibbs probabilities use the already-certified algebraic brackets. For logarithms it uses
\[
\log r=2\sum_{k=0}^{59}\frac{z^{2k+1}}{2k+1}+R,
\quad z=\frac{r-1}{r+1},\quad
|R|\leq\frac{2|z|^{121}}{121(1-z^2)}.
\]
All interval operations round outward onto a rational grid of mesh $10^{-50}$. Exact rational endpoint signs prove the displayed KL and TV crossing brackets. The script uses floating point only to print readable values. The root estimates shown on the figures are numerical; their enclosing crossing intervals are exact.

\subsection{Executed discovery calculations}''')
s=s.replace('The report proves an asymptotic effect and does not compute a specific finite crossing time or characterize every transient distance crossing.','The report proves an asymptotic effect and certifies finite crossing brackets for KL and total variation. It does not prove uniqueness of those crossings or characterize every transient distance crossing.')
s=s.replace('The corrected order-preserving conditions are sufficient, not necessary.','The corrected order-preserving conditions are sufficient, not necessary, and are explicitly attributed to established monotonicity methods.')
s=s.replace('For external review, the immediate useful package is this report plus the two exact certificate scripts and their JSON outputs.','For external review, the immediate useful package is this report, its figures, the two spectral/overlap certificates, the finite-time certificate, and the hypothesis-by-hypothesis literature comparison.')
s=s.replace('The corrected theorem is stated only with a simple slowest nonzero eigenvalue; no unstated limit to the degenerate case is used.','The spectral theorem assumes a simple slowest nonzero eigenvalue; no unstated limiting argument is used. The all-time divergence corollary separately covers degenerate spectra under its reversible monotonicity assumptions.')
s=s.replace('\\end{document}',r'''\clearpage
\section*{Appendix A: Explicit slow-mode cofactor polynomials}\label{app:cofactors}
The primary counterexample's exact left eigenvector is $w_i=P_i(r)/P_1(r)$, where $r$ is the simple slow root specified in Theorem~1 and
'''+polys+'\n'+companion+'\n\\end{document}')
s=s.replace(r'\path{theorem_notes.md}',r'\path{research/theorem_notes.md}').replace(r'\path{cone_search.py}',r'\path{research/cone_search.py}').replace(r'\path{adversarial_search.py}',r'\path{research/adversarial_search.py}').replace(r'\path{adversarial_notes.md}',r'\path{research/adversarial_notes.md}')
s=s.replace('The independent verifier checked the primary candidate', 'A separate model-assisted verification checked the primary candidate')
# Explicitly identify the published inference and its exact logical negation.
s=s.replace('The four-state implication specified in the assignment is false.', 'We disprove the four-state case of the finite-temperature triplet conjecture of Avitan, Factor and Gelbwaser-Klimovsky (arXiv:2603.04567v1).')
s=s.replace("The assignment's transcription is correct.", "This is a direct transcription of the published test, using destination-first rate indices.")
s=s.replace('The target was not changed.', 'The source-to-counterexample correspondence is made explicit in Section~\\ref{sec:claimscope}.')
s=s.replace(r'\subsection{What the cooling witness proves}', (ROOT/'research/claim_scope.tex').read_text()+'\n'+r'\subsection{What the cooling witness proves}\label{sec:cooling}')
s=s.replace('\\begin{theorem}\nThe four-triplet absence inference is false for four distinct energies, a positive finite bath temperature, strictly positive reversible rates and a simple slowest nonzero eigenvalue.', r'''\begin{theorem}[Counterexample to the published finite-temperature inference]\label{thm:counterexample}
There exist four distinct ordered energies, a finite positive bath temperature and a strictly positive reversible generator with a simple slowest nonzero eigenvalue such that all four original-rate predicates $P_{ijk}$ fail, yet two finite cooling preparations obey $|b_2(\beta_h)|<|b_2(\beta_c)|$. Thus the antecedent of $(\mathcal C_4)$ holds and its conclusion is false.''')
s=s.replace('There is no remaining mathematical gap in the specified four-state implication: the counterexample disproves it.', 'The explicit witness contradicts the published unrestricted finite-temperature implication $(\\mathcal C_4)$ under its stated original-rate interpretation. This conclusion concerns neither the separate zero-temperature theorem nor the simulation boxes discussed in Section~\\ref{sec:claimscope}.')
# A different exact verification path for the same primary generator.
s=s.replace(r'\subsection{Figure and finite-time verification}', (ROOT/'research/independent_verification.tex').read_text()+'\n'+r'\subsection{Figure and finite-time verification}')
s=s.replace('Rational interval arithmetic certifies the thermal-overlap inequality.', 'Rational interval arithmetic certifies the thermal-overlap inequality. A second standalone verification uses exact inertia counts and spectral residual bounds.')
s=s.replace('For external review, the immediate useful package is this report, its figures, the two spectral/overlap certificates, the finite-time certificate, and the hypothesis-by-hypothesis literature comparison.', 'For external review, the immediate useful package is this report, the standalone inertia/residual checker and its adversarial tests, the original spectral/overlap certificates, the finite-time certificate, figures, and the hypothesis-by-hypothesis literature comparison. Section~\\ref{sec:separatecheck} gives the alternative mathematical argument and distinguishes algorithmic independence from external peer review.')
# Broader novelty assessment: retain proofs while separating applications from priority claims.
s=s.replace('The supporting order argument is presented as an application of existing theory, with explicit comparison to Ohga, Hayakawa and Ito.', 'The counterexample is the central contribution; the order and convex-divergence results are supporting applications of established stochastic-monotonicity and thermomajorization methods.')
s=s.replace('We do not claim a new general monotonicity principle.\n\nFor definitions', 'The broader assessment in Section~\\ref{sec:novelty} also checks Fill--Kahn (v3, 6 September 2013), Baake--Sumner (v3, 25 September 2022), and Vu--Hayakawa (v2, 6 January 2025), including their relevant full-text arguments. Supporting applications carry no independent novelty claim.\n\nFor definitions')
s=s.replace(r'\section{Mathematical explanation and corrected theorem}', r'\section{Mathematical explanation and supporting applications}')
s=s.replace(r'\subsection{A corrected sufficient condition}', r'\subsection{Triplet-rate applications of stochastic monotonicity}')
s=s.replace('[Order-preserving sufficient condition]', '[Application of stochastic monotonicity]')
s=s.replace('[All-time divergence ordering under reversible monotonicity]', '[Application: all-time convex-divergence ordering]')
s=s.replace('particularly their comparison framework and Corollary 3.3. Their setup compares different kernels under suitable initial laws; the thermal preparations here differ instead.', 'specifically Proposition 3.2, Corollary 3.3 and the density-ratio identity (3.2). They compare different kernels from a common initial law and obtain TV, separation and, with additional assumptions, weighted $L^2$ comparisons. Here one kernel acts on different Gibbs preparations. Their argument is a methodological precedent, not a verbatim statement of the corollary below.')
old='The defensible contribution claim is the exact counterexample to the unrestricted triplet inference, its robust family and mechanism, and the explicit triplet-rate connection to an established sufficient condition. We do not present the general order-preservation mechanism or the elementary divergence consequence as newly discovered principles. The all-$\\kappa\\leq1$ statement is an application; its broader publication priority remains unestablished.'
assert old in s
s=s.replace(old,(ROOT/'research/novelty_assessment.tex').read_text().strip())
s=s.replace('The September paper also studies the fully degenerate reset point; general novelty of that point is not claimed.', 'The September paper also studies the fully degenerate reset point. Section~\\ref{sec:novelty} traces the general algebra to established equal-input generators; only its problem-specific use is needed here.')
# Predictive one-parameter family, with standalone inline figure data.
mc=list(csv.DictReader((ROOT/'output/figures/mechanism_curves.csv').open()))
mo=list(csv.DictReader((ROOT/'output/figures/mechanism_onset.csv').open()))
def mechanism_coords(rows,x,y):
    return ' '.join('('+r[x]+','+r[y]+')' for r in rows)
mf=r'''\begin{figure}[htbp]
\centering
\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.65cm},width=.43\textwidth,height=5.6cm,tick label style={font=\scriptsize},label style={font=\small},title style={font=\small},scaled ticks=false]
\nextgroupplot[title={A: A thermal maximum emerges},xlabel={Initial inverse temperature $\beta$},ylabel={$10^4(|b_z(\beta)|-|b_z(0)|)$},xmin=0,xmax=.05,ymin=-11.5,ymax=2.4,xtick={0,.01,.02,.03,.04,.05},xticklabel style={/pgf/number format/fixed,/pgf/number format/precision=2},legend style={font=\scriptsize,draw=none,at={(.02,.02)},anchor=south west}]
\addplot[blue,thick] coordinates {M025};\addlegendentry{$z=-0.025$}
\addplot[orange,thick] coordinates {MCRIT};\addlegendentry{$z=z_*$}
\addplot[green!60!black,thick] coordinates {M012};\addlegendentry{$z=-0.012$}
\addplot[red,thick] coordinates {M004};\addlegendentry{$z=-0.004$}
\addplot[gray,forget plot] coordinates {(0,0)(.05,0)};
\nextgroupplot[title={B: The onset boundary},xlabel={Highest-state mode value $z$},ylabel={$\beta_m$},xmin=-.03,xmax=-.004,ymin=0,ymax=.027,xtick={-.03,-.02,-.01},ytick={0,.01,.02},xticklabel style={/pgf/number format/fixed,/pgf/number format/precision=2},yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=2},legend style={font=\scriptsize,draw=none,at={(.02,.98)},anchor=north west}]
\addplot[fill=orange!25,draw=none,forget plot] coordinates {PHASE (-.004,0)} --cycle;
\addplot[black,thick] coordinates {PHASE};\addlegendentry{Unique maximum}
\addplot[gray,dashed,forget plot] coordinates {(-.0179401112234,0)(-.0179401112234,.027)};
\end{groupplot}
\end{tikzpicture}
\caption{Predictive family of Theorem~\ref{thm:family}. All four triplet predicates fail throughout. Above $z_*$ the amplitude initially increases with inverse temperature and has exactly one cooling maximum; below $z_*$ it decreases throughout cooling. The shaded region has $\partial_\beta|b_z|>0$, so any two finite inverse temperatures inside one vertical slice give a witness. The $\beta=0$ limit is only a reference for the left ordinate. Curves are numerical; the onset and single-turn classification have a continuous-parameter proof.}
\label{fig:mechanism}
\end{figure}
'''
for token,column in [('M025','z_minus025'),('MCRIT','critical'),('M012','z_minus012'),('M004','z_minus004')]:
    mf=mf.replace(token,mechanism_coords(mc,'beta',column))
mf=mf.replace('PHASE',mechanism_coords(mo,'z','peak_beta'))
mechanism=(ROOT/'research/predictive_mechanism.tex').read_text().replace('MECHANISM_FIGURE',mf)
s=s.replace(r'\subsection{Visible reversal and its size}',mechanism+'\n'+r'\subsection{Visible reversal and its size}')
s=s.replace('The example has a simple slow mode and strict margins, hence belongs to an open family of counterexamples.', 'The example has a simple slow mode and strict margins. A controlled family additionally gives an exact onset inequality, a unique thermal maximum when the effect occurs, and onset laws giving linear growth of the inverse-temperature window and quadratic growth of the amplitude advantage.')
s=s.replace('No assertion about a strong Mpemba effect is needed:', 'Section~\\ref{sec:mechanism} extends the companion into a parameter family that crosses a proved cooling threshold while all triplet tests continue to fail. No assertion about a strong Mpemba effect is needed:')
s=s.replace('What remains open in this work is a sharp replacement characterization,', 'Section~\\ref{sec:mechanism} now gives a sharp onset characterization for a controlled family and a single-turn eigenvector geometry. What remains open is an unrestricted replacement characterization,')
start=s.index('A concrete next calculation is to fix the rational bath weights above,')
stop=s.index('\n\nFor external review,',start)
s=s[:start]+r'''A concrete next calculation is to vary $q_4$ and $p$ in the construction of Section~\ref{sec:mechanism}, solving the same two rate equations and checking where positivity, triplet failure and spectral separation persist. The single-turn lemma would continue to classify cooling wherever its two energy/coefficient inequalities hold. This would test how the onset boundary and effect size change beyond the one-parameter family, rather than merely searching for another isolated witness.'''+s[stop:]
assert 'MECHANISM_FIGURE' not in s
assert 'assignment' not in s
(ROOT/'research_report.tex').write_text(s)
print('Revised standalone LaTeX written, including three inline-data figures.')
