"""Finite-test classification of simple-mode four-state thermal overlaps.

The mathematical theorem is recorded separately. This module is a floating-point
implementation for exploration, NOT a rigorous interval decision procedure.
Near-zero decisions are returned as 'boundary_or_uncertain', never rounded into
an asserted strict inequality. The bath defaults to inverse temperature one.
"""
import numpy as np
from scipy.optimize import brentq

def coordinates(E,beta):
    E=np.asarray(E,float)-E[0];a,b,c=E[1:]
    if np.isinf(beta):return -1.,0.
    f2,f3,f4=np.exp(-beta*E[1:]);r3=np.exp(-beta*(b-a));r4=np.exp(-beta*(c-a))
    den=a+b*r3+c*r4
    u=(-a+(b-a)*f3+(c-a)*f4)/den
    v=r3*(-b+(a-b)*f2+(c-b)*f4)/den
    return float(u),float(v)

def slope_constants(E):
    E=np.asarray(E,float)-E[0];a,b,c=E[1:]
    return a*b*(b-a),a*c*(c-a),b*c*(c-b),(c-a)/(c-b)

def slope_x(E,x):
    c0,c1,c2,gamma=slope_constants(E)
    return (c2*x**gamma-c0)/(c0+c1*x)

def classify_mode(E,w,bath=1.,tol=1e-10):
    E=np.asarray(E,float);w=np.asarray(w,float)
    if E.shape!=(4,) or w.shape!=(4,) or not np.all(np.diff(E)>0) or bath<=0:raise ValueError('Ordered four-state energies and positive bath required')
    E=E-E[0]
    A,B,C=w[:3]-w[3];scale=max(abs(A),abs(B),abs(C))
    if not scale:return {'cooling':False,'heating':False,'constant_mode':True}
    A,B,C=A/scale,B/scale,C/scale
    def phi(beta):
        u,v=coordinates(E,beta);return float(A+B*u+C*v)
    vertex=None;target=None;s0=slope_x(E,1.)
    if C!=0:
        target=-B/C
        if -1<target<s0:
            # Solve in y=-log(x), avoiding an absolute x-tolerance that could
            # round a very low-temperature root into the endpoint x=0.
            c0,c1,c2,gamma=slope_constants(E)
            logtarget=np.log1p(target)
            def equation(y):
                return (np.logaddexp(np.log(c2)-gamma*y,np.log(c1)-y)
                        -np.logaddexp(np.log(c0),np.log(c1)-y)-logtarget)
            upper=1.
            while equation(upper)>0:upper*=2
            root=brentq(equation,0.,upper,xtol=1e-13)
            vertex=root/(E[3]-E[2])
    endpoints={'zero':phi(0.),'bath':phi(bath),'infinity':float(A-B)}
    cool=[endpoints['zero'],endpoints['bath']];heat=[endpoints['bath'],endpoints['infinity']]
    if vertex is not None:
        val=phi(vertex)
        if vertex<bath:cool.append(val)
        elif vertex>bath:heat.append(val)
    def both_signs(values):
        low,high=min(values),max(values)
        if low < -tol and high > tol:return True
        if low > tol or high < -tol:return False
        # A zero endpoint does not itself establish a sign change; preserve ambiguity.
        return 'boundary_or_uncertain'
    return {'cooling':both_signs(cool),'heating':both_signs(heat),'cooling_values':cool,'heating_values':heat,'normalized_phi_endpoints':endpoints,'vertex_beta':vertex,'slope_target':target,'slope_at_zero':s0,'scope':'ordinary asymptotic effect for a simple slow mode; floating-point classification'}

if __name__=='__main__':
    import json
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    data=json.loads((root/'research/adversarial_rational_mode_certificate.json').read_text())
    from fractions import Fraction
    q=np.array([float(Fraction(x)) for x in data['q']]);w=np.array([float(Fraction(x[0])) for x in data['left_intervals']])
    print(json.dumps(classify_mode(-np.log(q),w),indent=2))
