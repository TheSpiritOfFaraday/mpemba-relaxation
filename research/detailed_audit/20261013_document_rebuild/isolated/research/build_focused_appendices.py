"""Build technical appendices from the preserved comprehensive report.

The report remains authoritative for retained results. This script only
changes document framing, section numbering and internal appendix names.
"""
from pathlib import Path
import re
import json

ROOT=Path(__file__).resolve().parents[1]
s=(ROOT/'research_report.tex').read_text()
preamble=s[:s.index(r'\title{')]
body=s[s.index(r'\section{Problem and literature status}'):s.index(r'\end{document}')]
original_labels=set(re.findall(r'\\label\{([^}]+)\}',body))
# Starred end sections become true lettered appendices with correct targets.
body=body.replace(r'\section*{Appendix A: Explicit slow-mode cofactor polynomials}',
                  r'\section{Explicit slow-mode cofactor polynomials}')
body=body.replace(r'\section*{Appendix B: A rational slow-mode companion}',
                  r'\section{The rational slow-mode companion}')
# The earlier report used these literal locations, which now have new letters.
body=body.replace(r'\hyperref[app:cofactors]{Appendix A}',r'Appendix~\ref{app:cofactors}')
body=body.replace(r'\hyperref[app:rational]{Appendix B}',r'Appendix~\ref{app:rational}')
body=body.replace('Appendix B',r'Appendix~\ref{app:rational}')
body=body.replace('Appendix A',r'Appendix~\ref{app:cofactors}')
body=body.replace(r'\section{Best next step}',r'\section{Further research directions}')
# Add a supplementary appendix for the distinct family now used in Lean.
sharp=r'''
\section{Explicit rates for the formal sharpness family}\label{app:formalsharp}
This family is the one used in the focused paper's sharp-supremum proof
and in \path{formal/Mpemba/SharpFamily.lean}. It is distinct from the
uniformly gapped construction in Appendix~\ref{sec:quarterlimit}.
For sufficiently small $\varepsilon>0$, use
\[
q=(1,1-\varepsilon/3,3/5,\varepsilon^2),\qquad
w=(1,-4/5,-1/3-4\varepsilon/9+5\varepsilon^3/3,-\varepsilon).
\]
Set $a_{12}=1$ and define
\begin{align*}
P&=15\varepsilon^4-5\varepsilon^3+69\varepsilon^2-20\varepsilon+68,\\
R&=15\varepsilon^3-4\varepsilon-12,\qquad
S=75\varepsilon^3-20\varepsilon+21,\\
a_{13}&=\frac{150\varepsilon^6+640\varepsilon^5+190\varepsilon^4+
459\varepsilon^3+138\varepsilon^2-252\varepsilon-1224}{RP},\\
a_{14}&=\frac{30\varepsilon^3+125\varepsilon^2-42\varepsilon+240}{3P},\\
a_{23}&=\frac{5(3-\varepsilon)(225\varepsilon^5-185\varepsilon^4+
425\varepsilon^3-113\varepsilon^2-60\varepsilon+204)}{SP},\\
a_{24}&=\frac{5(3-\varepsilon)(9\varepsilon^2-2\varepsilon+16)}{3P},
\qquad a_{34}=\frac{3(9\varepsilon^2-2\varepsilon+16)}{P}.
\end{align*}
Detailed balance determines the upward rates and the diagonals conserve
probability. Rational identities give $q^\top w=0$, $M^\top w=-3w$,
$F_{124}=0$ and $G_{234}=0$. To make the continuity argument precise,
write $(a,b,c,d,e,f)=(a_{12},a_{13},a_{14},a_{23},a_{24},a_{34})$ and
$(t,u,v)=(q_2,q_3,q_4)$, and define the cleared expressions
\begin{align*}
\widehat F_{123}&=bt(2+u)-d(2t+u)-2at(1-t),\\
\widehat F_{134}&=cu(2+v)-f(2u+v)-2bu(1-u),\\
\widehat F_{234}&=eu(2t+v)-ft(2u+v)-2du(t-u).
\end{align*}
These have the same signs as the corresponding original $F$ expressions.
The following thirteen rational functions are continuous and strictly
positive at $\varepsilon=0$:
\[
P,-R,S,b,c,d,e,f,d-bt,-\widehat F_{123},-\widehat F_{134},
-\widehat F_{234},\frac{2(\varepsilon-3)(15\varepsilon^2-5\varepsilon-6)}{9P}.
\]
They remain positive on a nonempty interval of positive $\varepsilon$.
The last expression multiplied by $\varepsilon$ equals $e-ct$, giving
$G_{124}<0$ there. The equality $eu=ft$ and $e-ct>0$ imply $cu<f$,
so $G_{134}<0$ as well. All original triplets therefore fail.
The positive contrast coordinates of $w$ and two connecting positive
off-diagonal contrast entries establish that $-3$ is the simple largest
contrast eigenvalue. The formal proof transfers this to the full
generator; it does not assume a numerical eigensystem or a uniform gap.
Appendix~\ref{sec:lean} gives the formal theorem names and reproduction
instructions. The limits with $\varepsilon=e^{-L^2}$,
$\beta_h=L^{-3}$ and $\beta_c=L^{-1}$ complete the sharpness proof.
'''
assert original_labels <= set(re.findall(r'\\label\{([^}]+)\}',body))
intro=r'''
\title{Technical appendices to\\Collective Mpemba cooling beyond triplet tests}
\author{}
\date{12 October 2026}
\begin{document}
\maketitle
These appendices accompany the focused paper on the counterexample,
mechanism and sharp $25\%$ amplitude bound. They preserve the detailed
proofs and additional results of the comprehensive research report.
Their theorem, equation and figure numbers are internal to this document.
The main paper leads with the rational witness; the retained original
primary witness below has two different rates and its own numerical plots.
Neither witness's finite-time data should be substituted for the other's.

\paragraph{Reading guide.}
The source-predicate correspondence is in Appendix~\ref{sec:claimscope};
the controlled mechanism is in Appendix~\ref{sec:mechanism};
the eigenfunction restriction and sharp upper bound are in
Appendix~\ref{sec:fulltripletbound}. The uniformly gapped family is in
Appendix~\ref{sec:quarterlimit}, the rational witness in
Appendix~\ref{app:rational}, and the simpler formal sharpness rates in
Appendix~\ref{app:formalsharp}. Lean coverage is specified in
Appendix~\ref{sec:lean}. Classifications, higher-dimensional results,
finite-time extensions and their limitations are retained in the other
appendices, rather than being additional claims of the focused main text.
\clearpage
\tableofcontents
\clearpage
\appendix
'''
out=preamble+intro+body+sharp+'\n\\end{document}\n'
(ROOT/'focused_appendices.tex').write_text(out)
receipt={'source':'research_report.tex','retained_labels':len(original_labels),
 'all_original_labels_retained':True,'added_appendix':'Explicit rates for the formal sharpness family',
 'scope':'All technical report sections retained from Problem and literature status onward.'}
(ROOT/'research/focused_appendices_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
