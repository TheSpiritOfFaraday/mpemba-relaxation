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

# Extended classification and a surviving all-dimension triplet implication.
new_sections='\n'.join((ROOT/'research'/name).read_text() for name in [
    'ordered_conductance_theorem.tex','tridiagonal_contrast_entropy.tex',
    'thermal_ratio_obstruction.tex','near_bath_entropy_window.tex',
    'energy_threshold.tex','amplitude_bound.tex',
    'universal_classification.tex','degenerate_extension.tex','commensurate_classification.tex',
    'finite_chi_classification.tex','mode_realizability.tex',
    'quarter_limit_family.tex','full_triplet_amplitude_bound.tex','finite_span_amplitude.tex',
    'full_triplet_endpoint_criterion.tex','one_fast_mode_classification.tex',
    'rate_endpoint_criterion.tex','rational_log_endpoint.tex',
    'higher_dimension_extension.tex','high_dimension_amplitude.tex',
    'high_dimension_bound_reduction.tex','reset_trace_obstruction.tex',
    'general_n_arithmetic_lp.tex'])
s=s.replace(r'\section{Approaches explored}',new_sections+'\n'+r'\section{Approaches explored}')
quarter_rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((ROOT/'output/figures/quarter_overlap_data.csv').open())]
quarter_max=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((ROOT/'output/figures/quarter_gain_data.csv').open())]
quarter_figure=r'''\begin{figure}[htbp]
\centering
\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.6cm},width=.43\textwidth,height=5.4cm,tick label style={font=\scriptsize},label style={font=\small},title style={font=\small},ymin=0,ymax=26]
\nextgroupplot[title={A: Finite-temperature amplitude advantage},xmin=0,xmax=.35,xlabel={Initial inverse temperature $\beta$},ylabel={Advantage (\%)},legend style={font=\scriptsize,draw=none,at={(.98,.98)},anchor=north east}]
\addplot[blue,thick] coordinates {CURVE3};\addlegendentry{$\varepsilon=10^{-3}$}
\addplot[orange,thick] coordinates {CURVE10};\addlegendentry{$\varepsilon=10^{-10}$}
\addplot[green!50!black,thick] coordinates {CURVE50};\addlegendentry{$\varepsilon=10^{-50}$}
\nextgroupplot[title={B: Effect size and energy separation},xlabel={Bath-scaled energy span},ylabel={Maximum advantage (\%)},xmin=0,xmax=240]
\addplot[black,thick] coordinates {MAXIMA};
\addplot[gray,dashed] coordinates {(0,25)(240,25)};
\node[font=\scriptsize,anchor=south east] at (axis cs:237,25) {Proved $25\%$ limit};
\end{groupplot}
\end{tikzpicture}
\caption{The unperturbed family of Theorem~\ref{thm:quarterlimit}. Panel A plots $100[1-|b_2(\beta_h)|/|b_2(\beta)|]$ with finite $\beta_h=[\log(1/\varepsilon)]^{-2}$. Panel B follows its thermal maximum as the energy span $2\log(1/\varepsilon)$ grows. The substantial limiting advantage is accompanied by increasing energy separation; the plot does not claim practical cooling gains at modest spans. Curves and maxima use 80-digit arithmetic and illustrate the analytic theorem. The separate certified $22.617\%$ witness uses a different fixed hot preparation. The strict-triplet perturbation preserves the limiting result but is not plotted.}
\label{fig:quarter}
\end{figure}
'''
for exponent in [3,10,50]:
    quarter_figure=quarter_figure.replace('CURVE'+str(exponent),coords([r for r in quarter_rows if int(r['exponent'])==exponent],'beta','gain',2))
quarter_figure=quarter_figure.replace('MAXIMA',coords(quarter_max,'energy_span','maximum_gain',2))
s=s.replace(r'\section{The sharp four-state amplitude bound under all triplet tests}',quarter_figure+'\n'+r'\section{The sharp four-state amplitude bound under all triplet tests}')
s=s.replace(r'\title{A counterexample to the four state Mpemba triplet conjecture}',r'\title{Mpemba relaxation beyond triplet tests}')
a=s.index(r'\begin{abstract}');b=s.index(r'\end{abstract}',a)
s=s[:a]+r"""\begin{abstract}
We disprove the four-state finite-temperature triplet conjecture of Avitan, Factor and Gelbwaser-Klimovsky. An explicit positive reversible generator fails the prescribed necessary-condition test on every original-rate triplet but exhibits ordinary cooling between two finite positive initial temperatures. Exact arithmetic certifies the slow mode and overlap reversal; a controlled family gives a sharp onset threshold. We then prove what survives of the triplet inference in any dimension: its $G$ inequalities exclude all strong effects and asymptotic heating. The proof makes the generator a symmetric Metzler matrix in a weighted prefix-contrast basis. It also excludes finite-time heating reversals in weighted Euclidean distance, including degenerate slow spaces. A sharp energy-and-bath threshold characterizes universal absence of ordinary cooling across this larger rate family. For arbitrary four-state generators with a simple slow mode, we give a necessary-and-sufficient asymptotic cooling/heating test using endpoint signs and at most one uniquely determined additional temperature. An explicit Wronskian identity reduces the problem to a strictly convex curve and implies the sharp $n-2$ bound on stationary points of a Gibbs expectation. With a slow eigenspace of dimension $r$, the number of distinct nonbath strong preparations is at most $n-r-1$, and the bound is attainable. Degenerate four-state modes admit a separate finite classification for asymptotic relative entropy and weighted Euclidean distance. The results separate ordinary cooling, strong cancellation, heating, and transient distance-dependent behavior. Established order theory, Chebyshev systems, and spectral realization are attributed explicitly.
"""+s[b:]
needle=r'\section{Problem and literature status}'
s=s.replace(needle,r"""Sections~\ref{sec:fulltripletbound} and~\ref{sec:tripletendpoint} establish the sharp $25\%$ amplitude supremum and a single-endpoint cooling criterion under the original triplet premise. The new exclusion theorem in Section~\ref{sec:orderedconductances} shows that the counterexample has a restricted mechanism: it reverses ordinary cooling while retaining a nonzero slow overlap. The same triplet premise rules out strong effects in every dimension and rules out asymptotic heating. Section~\ref{sec:universal} supplies a full finite asymptotic classification for arbitrary simple four-state slow modes; Section~\ref{sec:degenerateextension} states the distance-specific extension to other slow multiplicities. These conclusions have mathematical proofs rather than being inferred from the discovery searches.

"""+needle)
s=s.replace('The slowest nonzero eigenvalue is additionally assumed simple and is certified below.', 'The counterexample and scalar-overlap classification additionally assume a simple slowest nonzero eigenvalue. Section~\\ref{sec:degenerateextension} treats the other multiplicities for specified asymptotic distances, while Theorem~\\ref{thm:orderedk} does not require simplicity.')
a=s.index(r'\section{Remaining gap}');b=s.index(r'\section*{Appendix A:',a)
s=s[:a]+r"""\section{Remaining gap}
The original ordinary-cooling implication is false. Its strong-effect and asymptotic-heating exclusions survive in every dimension by Theorem~\ref{thm:orderedk}, under weaker hypotheses. Theorem~\ref{thm:finiteclassification} classifies all simple four-state modes once the full slow eigenvector is known, with a separate finite prescription for degenerate modes and specified distances.

The inverse problem has a separate exact solution at fixed rational bath weights and rational slow-mode direction: Theorem~\ref{thm:inversefeasibility} reduces realizability under every triplet failure to eight half-planes intersecting an open disk. A continuously certified family approaches the sharp $25\%$ amplitude-advantage supremum, with a uniform spectral gap and a perturbation making all triplet inequalities strict. The new upper bound follows from a restriction on the first three components of the slow eigenfunction. Within this full triplet-failure family, ordinary cooling occurs exactly when the uniform energy--mode covariance is positive; any cooling maximum is unique. The broader $G$-only contrast family retains its different sharp bound $1-2/n$.

What remains open here is a comparably sharp criterion directly in the original rates that avoids computing the full slow eigenspace, an efficient certified implementation for every transcendental equality case, and a classification of transient reversals for distances such as relative entropy and total variation. The finite-time heating exclusion proved here is for weighted Euclidean distance. Neither the scalar classification nor the sixteen-term degenerate prescription classifies all transient crossings. Supporting classical tools and the rank-one realization are not claimed as new; publication priority for the problem-specific consequences requires further expert review.

\section{Best next step}
The most useful next calculation is to test whether the contrast-cone heating exclusion extends from weighted Euclidean distance to finite-time relative entropy, or whether a rigorous transient heating example separates them. A completed 1200-second search evaluated 6,740,275,200 sampled derivatives in dimensions 4, 5, 6 and 8 and found none negative; this is search evidence, not a proof. The full four-state amplitude supremum and the location of its possible anomalous interval are now settled. A separate next direction is a sharp amplitude bound for larger state counts under all original $F,G$ inequalities. These questions exploit the new structure rather than searching again for strong effects that the theorem now excludes.

For external review, the package includes explicit generators, two distinct primary counterexample checks, finite-time certificates, the all-dimension contrast proof and symbolic checks, the finite classifier with polynomial-reference audits, and exact sharpness constructions. These are independent calculations within one project, not external human peer review or proof-assistant formalization.

"""+s[b:]

s=s.replace('The results separate ordinary cooling, strong cancellation, heating, and transient distance-dependent behavior.', 'An exact inverse test reduces realization of a prescribed four-state slow mode under the triplet premise to eight half-planes intersecting a disk. A uniformly gapped family gives cooling amplitude advantages approaching $25\\%$, including after perturbation into strict triplet failure; the larger contrast family obeys a sharp dimension-dependent amplitude bound. The results separate ordinary cooling, strong cancellation, heating, and transient distance-dependent behavior.')
s=s.replace('For four states there are at most six distinct exponential terms. Chebyshev/variation arguments could bound turning points, but cannot exclude the certified turning behavior. No new general zero-count theorem is claimed. This route became explanatory, not decisive.', 'A direct pairwise expansion has at most six terms in four states, but the Gibbs Wronskian representation is sharper: it gives at most $n-2$ stationary points and makes the four-state derivative curve strictly convex. This is the decisive route to the finite classification in Section~\\ref{sec:universal}. The Chebyshev-system machinery is classical; the explicit identity and its application are proved here.')
s=s.replace(r'\subsection{Boundary audit}',r'''\subsection{Verification of the expanded results}
The general contrast identity is checked symbolically for dimensions two through seven in
\texttt{ordered\_conductance\_certificate.py}; the displayed all-dimension proof does not
depend on those finite checks. A separate numerical audit compares 5,000 heating
preparations across 1,000 generators in dimensions 3, 4, 5, 7 and 10, with no violation
of the proved weighted-Euclidean ordering.

The classifier audit records 20,000 side classifications, 600 comparisons against
exact polynomial Sturm tests, and 40 high-precision Wronskian checks, with no remaining
discrepancy. The tests use independent covariance evaluations with exponential scaling
to avoid underflow. They validate an implementation; the general classification rests
on the analytic proof. The energy-threshold and degenerate-mode scripts give exact
rational sharpness examples. Floating evaluation alone does not certify equality cases.

The exact commensurate classifier passes sixteen boundary, multiplicity and algebraic
embedding checks. A cross-audit of 100 systems compares 200 side decisions with a
separate geometric classifier: 70-digit symmetric diagonalization supplies modes,
which are then converted to floating-point inputs for that independent implementation.
All comparisons agree, with no unresolved boundary decision. This cross-audit exposed
and led to repair of a scaled-algebraic-root representation bug; its first failure log
is retained. The tests support the implementation and do not replace its exact Sturm
argument. The pure-transient pair and five-state extension have separate exact
certificates, including two crossing roots and ten strict triplet failures, respectively.

The exact inverse constructor is audited with ten checks, including affine invariance
of the proposed mode, rejection of excluded strong-mode directions and invalid inputs,
and independent rational LDL inertia of the resulting generators. The inverse
constructor's zero-slack setting is complete for the original inequalities; positive
slack strengthens them and cannot certify original infeasibility. The quarter-limit
family is certified over an entire parameter interval by 30 exact Bernstein sign
certificates. Its finite thermal witness uses integer-power root enclosures, with an
analytic perturbation bound for the strict-triplet variant. The amplitude upper bound
and its sharpness construction have direct analytic proofs. Commands, dependencies and
individual scope limitations are listed in \texttt{RUNME.txt}.

\subsection{Boundary audit}''')

a=s.index(r'\begin{abstract}');b=s.index(r'\end{abstract}',a)
s=s[:a]+r'''\begin{abstract}
We disprove the four-state finite-temperature triplet conjecture of Avitan, Factor and Gelbwaser-Klimovsky. Positive reversible generators can fail every prescribed original-rate triplet test yet exhibit ordinary cooling. We prove that the slow-amplitude advantage is strictly below $25\%$, and construct a uniformly gapped family approaching this sharp supremum, including with every triplet inequality strict. Within the triplet-failure family, ordinary cooling occurs exactly when one endpoint covariance is positive, and its thermal maximum is unique. We prove a corrected implication in any dimension: the triplet premise excludes strong effects and asymptotic heating, and excludes finite-time heating in weighted Euclidean distance. The proof uses positivity in a weighted prefix-contrast basis and gives a sharp energy/bath threshold for universal cooling absence in a larger rate family. For arbitrary four-state reversible systems, we give a finite asymptotic classification: an explicit Gibbs Wronskian identity reduces simple slow modes to endpoint signs and at most one additional temperature, while degenerate modes have a separate distance-specific test. Commensurate energies and rational rates admit an exact executable classification, including algebraic eigenvectors and equality cases. For any fixed thermal pair, a finite-time weighted-Euclidean reversal is decided by endpoints and at most one additional time. Finally, realization of a prescribed rational slow mode under the original triplet premise reduces to eight half-planes intersecting an open disk. The results distinguish local rate tests, global spectral geometry, ordinary cooling, strong cancellation, and transient behavior. Classical spectral, order, and zero-counting tools, and prior scalar strong-temperature bounds, are explicitly attributed.
'''+s[b:]
s=s.replace('We prove a corrected implication in any dimension:', 'A dimension-extension theorem gives open sets of counterexamples for every state count $n\\geq4$. We prove a corrected implication in any dimension:')
s=s.replace('an efficient certified implementation for every transcendental equality case, and a classification of transient reversals for distances such as relative entropy and total variation.', 'an efficient certified implementation for arbitrary incommensurate transcendental equality cases, and a classification of transient reversals for distances such as relative entropy and total variation. Sections~\\ref{sec:exactclassifier} and~\\ref{sec:finitechi} settle the exact commensurate implementation and the specified-pair finite-time weighted-Euclidean test, respectively.')

s=s.replace('Within the triplet-failure family, ordinary cooling occurs exactly when one endpoint covariance is positive, and its thermal maximum is unique.',
    'Within the triplet-failure family, ordinary cooling occurs exactly when one endpoint covariance is positive, and its thermal maximum is unique. An equivalent scalar rate test requires one cubic root. Every repeated-slow four-state generator in this family is a reset plus mixing on one adjacent pair, and admits no convex-divergence reversal at any time.')
s=s.replace('These conclusions have mathematical proofs rather than being inferred from the discovery searches.',
    'Theorem~\\ref{thm:tripletdegenerate} completes the full-premise degenerate boundary: reset plus one adjacent pair is necessary, and excludes all-time convex-divergence reversals. The same classification holds in every dimension $n\\geq4$ when there is just one faster mode. Corollary~\\ref{cor:rateendpoint} implements the simple-mode cooling decision without an eigenvector. These conclusions have mathematical proofs rather than being inferred from the discovery searches.')
s=s.replace('What remains open here is a comparably sharp criterion directly in the original rates that avoids computing the full slow eigenspace, an efficient certified implementation',
    'The scalar rate test in Section~\\ref{sec:rateendpoint} avoids computing the full slow eigenvector, but retains the largest nonzero cubic root. A comparably sharp criterion using local original-rate inequalities alone is not established. Other remaining tasks are an efficient certified implementation')
s=s.replace('The finite-time heating exclusion proved here is for weighted Euclidean distance.',
    'The general contrast-cone finite-time heating exclusion is for weighted Euclidean distance. The explicitly classified degenerate full-premise systems separately exclude reversals for every convex divergence at every time.')
s=s.replace('The full four-state amplitude supremum and the location of its possible anomalous interval are now settled.',
    'A subsequent 3602-second targeted four-state optimizer made 66,396,528 objective calls and found no negative score. Its best candidate remained positive under an independent 100-digit direct matrix exponential; neither search proves absence. The full four-state amplitude supremum and the location of its possible anomalous interval are now settled.')
s=s.replace('individual scope limitations are listed in \\texttt{RUNME.txt}.',
    'individual scope limitations are listed in \\texttt{RUNME.txt}. The one-fast-mode boundary classification has 1,694 exact symbolic checks covering its generic identities and converse pairs in dimensions four through eight; its all-dimension support proof is separate from those finite checks. The scalar rate test is a direct spectral-projector identity, not a new general projector theorem.')
s=s.replace('Commensurate energies and rational rates admit an exact executable classification, including algebraic eigenvectors and equality cases.',
    'Commensurate energies and rational rates admit an exact executable classification, including algebraic eigenvectors and equality cases. Under the full triplet premise, Baker\'s theorem additionally gives a terminating exact rational-input decision without commensurate energies.')
s=s.replace('an efficient certified implementation for arbitrary incommensurate transcendental equality cases',
    'an efficient certified implementation for general nonrational input and its transcendental equality cases')
s=s.replace('settle the exact commensurate implementation and the specified-pair finite-time weighted-Euclidean test, respectively.',
    'settle the exact commensurate implementation and the specified-pair finite-time weighted-Euclidean test, respectively. Section~\\ref{sec:rationallogendpoint} removes the commensurability restriction for exact rational inputs under the full triplet premise, using Baker\'s logarithmic independence theorem.')
s=s.replace('A dimension-extension theorem gives open sets of counterexamples for every state count $n\\geq4$.',
    'For every $n\\geq4$, a construction with $n-3$ high-energy states attains amplitude advantages approaching $1-3/n$, so no dimension-independent bound below one survives; a finite five-state witness exceeds $37.67\\%$.')
s=s.replace('A separate next direction is a sharp amplitude bound for larger state counts under all original $F,G$ inequalities.',
    'For larger state counts, Theorem~\\ref{thm:highdimamplitude} supplies the achievable limit $1-3/n$ under all original tests, whereas the general contrast upper bound remains $1-2/n$. A sharp full-premise bound for $n>4$ remains open between these bounds.')
s=s.replace('A sharp full-premise bound for $n>4$ remains open between these bounds.',
    'A sharp full-premise bound for $n>4$ remains open between these bounds. Proposition~\\ref{prop:highdimboundreduction} reduces a matching upper bound to the slow-mode inequality $w_1+w_2\\geq2w_3$. Section~\\ref{sec:resettraceobstruction} gives an exact obstruction to one natural proof route: even after subtraction of its minimum reset, a generator satisfying every original triplet inequality can have a three-state trace violating $F_{123}\\leq0$. The arithmetic implication itself remains unresolved in higher dimension.')
s=s.replace('The most useful next calculation is to test whether the contrast-cone heating exclusion extends from weighted Euclidean distance to finite-time relative entropy, or whether a rigorous transient heating example separates them.',
    'The principal remaining transient question is finite-time relative-entropy heating order for general contrast-Metzler generators outside the quantified bath neighborhood. Theorems~\\ref{thm:contrastchainentropy} and~\\ref{thm:nearbatheentropy} now settle, respectively, all heating temperatures for tridiagonal contrast generators and both bath sides within $|\\beta-\\beta_b|(E_n-E_1)\\leq2/\\sqrt5$. The general case is still open here. Section~\\ref{sec:thermalratioobstruction} proves that the monotone-ratio step of the tridiagonal proof cannot simply be assumed for general admissible dynamics.')
s=s.replace('The general contrast-cone finite-time heating exclusion is for weighted Euclidean distance.',
    'The unrestricted contrast-cone finite-time heating exclusion is for weighted Euclidean distance. Finite-time KL ordering is additionally proved on a uniform neighborhood of the bath, on both sides, while tridiagonal contrast generators preserve heating order for every convex divergence at all finite temperatures.')
s=s.replace('An equivalent scalar rate test requires one cubic root.',
    'A finite-span refinement gives necessary energy and preparation-temperature separations for a specified advantage. An equivalent scalar rate test requires one cubic root.')
s=s.replace('The proof uses positivity in a weighted prefix-contrast basis and gives a sharp energy/bath threshold for universal cooling absence in a larger rate family.',
    'The proof uses positivity in a weighted prefix-contrast basis and gives a sharp energy/bath threshold for universal cooling absence in a larger rate family. Finite-time relative entropy has a uniform normal-ordering window on both bath sides; tridiagonal contrast generators preserve heating order for every convex divergence.')
s=s.replace(r'\newtheorem{corollary}{Corollary}',
    r'\newtheorem{corollary}{Corollary}'+'\n'+r'\newtheorem{proposition}{Proposition}')
a=s.index(r'\begin{abstract}');b=s.index(r'\end{abstract}',a)
s=s[:a]+r'''\begin{abstract}
We disprove the four-state finite-temperature triplet conjecture of Avitan, Factor and Gelbwaser-Klimovsky with a strictly positive reversible generator and an exact certificate. Under the conjecture's original triplet premise, the ordinary cooling amplitude advantage has sharp supremum $25\%$. A uniformly gapped family approaches this limit with all triplet inequalities strict. Cooling occurs exactly when an endpoint covariance is positive; an equivalent rate test uses one cubic root. Every repeated-slow four-state generator under the premise is a reset plus mixing on one adjacent pair and has normal convex-divergence ordering at every time. In arbitrary dimension, the triplet premise excludes strong effects and asymptotic heating, yet permits ordinary cooling advantages approaching $1-3/n$. Thus no dimension-independent amplitude bound below one survives. A weighted prefix-contrast representation also yields finite-time heating exclusion in weighted Euclidean distance, a uniform finite-time relative-entropy ordering window on both bath sides, and all-temperature heating order for every convex divergence in a structured generator family. For arbitrary simple four-state slow modes, an explicit Gibbs Wronskian identity gives a finite asymptotic classification; degenerate modes have a separate distance-specific test. Exact executable classifications cover commensurate rational systems and, under the triplet premise, rational rates without commensurate energies. We distinguish these proved results from unresolved general finite-time questions and the higher-dimensional sharp bound. Classical order, spectral, total-positivity, and transcendence tools are explicitly attributed.
'''+s[b:]
s=s.replace('This is an independent exact realization of the same thermal mechanism, constructed backward from a prescribed mode.', '')
# Keep the Lean coverage statement in the reproducible manuscript build.
s=s.replace(r'\date{11 October 2026}', r'\date{12 October 2026}')
s=s.replace('Classical order, spectral, total-positivity, and transcendence tools are explicitly attributed.',
    'Lean verifies a rational-mode counterexample, its eventual Pearson chi-square relaxation reversal, and the sharp four-state amplitude supremum; the formalization has explicitly delimited coverage and reproducible build instructions. Classical order, spectral, total-positivity, and transcendence tools are explicitly attributed.')
outcome_marker=r'\section{Problem and literature status}'
assert s.count(outcome_marker)==1
s=s.replace(outcome_marker, r'''The Lean companion verifies the rational-mode example in \hyperref[app:rational]{Appendix B}, an initial-to-eventual reversal in Pearson chi-square distance under its master equation, and the sharp four-state amplitude supremum. Section~\ref{sec:lean} maps these claims to named theorems and gives reproduction commands. Its sharpness proof uses a simpler family with two triplet equalities; the report's stronger uniform-gap and all-strict perturbation claims are outside that formal coverage.

'''+outcome_marker)
lean_marker=r'\subsection{Figure and finite-time verification}'
assert s.count(lean_marker)==1
s=s.replace(lean_marker,(ROOT/'research/lean_verification.tex').read_text()+'\n'+lean_marker)
s=s.replace('These are independent calculations within one project, not external human peer review or proof-assistant formalization.',
    r'The package now also includes the Lean verification described in Section~\ref{sec:lean}. The remaining analytic proofs and executable checks are not thereby formalized. All work belongs to the same project; no external human peer review is claimed.')
s=s.replace(r'\section*{Appendix B: A rational slow-mode companion}',
    r'''\section*{Appendix B: A rational slow-mode companion}\label{app:rational}
The data and strict overlap reversal below are verified in Lean, together with an initial-to-eventual Pearson chi-square relaxation reversal. Section~\ref{sec:lean} states the exact coverage and distinguishes the Lean enclosures from the tighter Python-certified intervals printed here.
''')
# Preserve the new physical research priority in every report rebuild.
priority_marker=r'\section{Best next step}'
assert s.count(priority_marker)==1
s=s.replace(priority_marker,
    (ROOT/'research/physical_finite_time_goal.tex').read_text()
    + '\n' + (ROOT/'research/physical_finite_time_results.tex').read_text()
    + '\n' + (ROOT/'research/physical_robustness_certificate.tex').read_text()
    + '\n' + (ROOT/'research/physical_finite_time_figure.tex').read_text()
    + '\n' + (ROOT/'research/thermal_collision_realization.tex').read_text()
    + '\n' + (ROOT/'research/activated_intermediate_study.tex').read_text()
    + '\n' + r'\subsection{Other open mathematical questions}')
s=s.replace('The principal remaining transient question is finite-time relative-entropy heating order',
    'A separate open mathematical question is finite-time relative-entropy heating order')
s=s.replace(r'\title{Mpemba relaxation beyond triplet tests}',
    r'\title{Mpemba relaxation beyond triplet tests\\[0.4em]\large Comprehensive research report}')
# Add the teaching course without renumbering any retained technical result.
# Its starred headings and unnumbered displays preserve the original counters.
course=(ROOT/'research/foundational_course.tex').read_text()
assert course.count('% PHYSICAL_TEACHING_BRIDGE')==1
course=course.replace('% PHYSICAL_TEACHING_BRIDGE', (ROOT/'research/physical_teaching_bridge.tex').read_text())
anchor=r'\maketitle'
assert s.count(anchor)==1
s=s.replace(anchor, anchor+'\n'+course, 1)
(ROOT/'research_report.tex').write_text(s)
print('Revised standalone LaTeX written, including inline-data figures and the physical-model study.')
