"""Check the additive tutorial revision and its elementary worked algebra.

This is scoped validation, not a re-proof of every advanced theorem.
"""
from pathlib import Path
import hashlib
import json
import re
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'research/detailed_audit/pedagogical_revision'
BASELINE = 'dd5b3062ab2320b062862cec650c4d68a6a52cac973912e3500c590d15280d81'
# The source/approach and historical-scope reviews made explicit prose edits.
# Their saved bounded diffs preserve the original baseline and pin the result;
# do not silently accept arbitrary current technical edits as a new baseline.
REVIEWED_TECHNICAL = '13935e1690d44fd4fa467f2ee10a154db79cc76aa12b9358c3c3bd003ea656fb'
checks = []


def require(ok, label):
    if not ok:
        raise ValueError(label)
    checks.append(label)


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    before = (OUT / 'research_report_before.tex').read_text()
    course = (ROOT / 'research/foundational_course.tex').read_text()
    require(course.count('% PHYSICAL_TEACHING_BRIDGE') == 1,
            'Exactly one physical teaching bridge marker')
    course = course.replace('% PHYSICAL_TEACHING_BRIDGE',
                            (ROOT / 'research/physical_teaching_bridge.tex').read_text())
    after = (ROOT / 'research_report.tex').read_text()
    require(hashlib.sha256(before.encode()).hexdigest() == BASELINE,
            'Preserved report matches fixed pre-teaching baseline')
    require(after.count(course) == 1, 'Exactly one copy of the teaching material')
    anchor = r'\begin{abstract}'
    pre_edit = (ROOT / 'research/detailed_audit/before_source_approach_research_report.tex').read_text()
    require(pre_edit[pre_edit.index(anchor):] == before[before.index(anchor):],
            'Original technical body was byte-preserved through the 40-lesson expansion')
    require(hashlib.sha256(after[after.index(anchor):].encode()).hexdigest() == REVIEWED_TECHNICAL,
            'Current technical body matches the explicitly reviewed source/approach and historical-scope edits')
    require(after[:after.index('\\maketitle')] == before[:before.index('\\maketitle')],
            'Original preamble and title are preserved')
    require(re.findall(r'\\section\*\{Lesson (\d+)\.', course) ==
            list(map(str, range(1, 41))), 'Lessons 1 through 40 appear in order')
    require(not re.search(r'\\begin\{(?:theorem|lemma|corollary|proposition|equation|align|figure)\}', course),
            'No numbered theorem, equation, align or figure counter is consumed')
    reviews = json.loads((ROOT / 'research/detailed_audit/claim_reviews.json').read_text())
    statements = []
    for m in re.finditer(r'\\begin\{(theorem|lemma|corollary|proposition)\}(?:\[([^\]]+)\])?', after):
        end = after.index('\\end{' + m[1] + '}', m.end())
        block = after[m.start():end + len(m[1]) + 6]
        statements.append(hashlib.sha256(block.encode()).hexdigest())
    for key, review in reviews.items():
        if key.startswith('research_report.tex:'):
            require(review['source_sha256'] == digest(ROOT / 'research_report.tex') and
                    review['statement_sha256'] == statements[int(key.split(':')[1]) - 1],
                    'Prior review still matches unchanged statement: ' + key)
    R = s.Rational
    q = [R(1), R(199, 200), R(2, 5), R(1, 400)]
    require(sum(q) == R(959, 400), 'Primary weight sum 959/400')
    require([x/sum(q) for x in q] == [R(400,959),R(398,959),R(160,959),R(1,959)],
            'Primary normalized probabilities')
    p, e = [R(1,2),R(1,3),R(1,6)], [0,1,2]
    mean = sum(x*y for x,y in zip(p,e))
    require(mean == R(2,3) and sum(x*y*y for x,y in zip(p,e))-mean**2 == R(5,9),
            'Lesson 4 mean and variance')
    p, e = [R(4,7),R(2,7),R(1,7)], [0,1,2]
    mean = sum(x*y for x,y in zip(p,e))
    require(mean == R(4,7) and sum(x*y*y for x,y in zip(p,e))-mean**2 == R(26,49),
            'Lesson 7 log-two-scaled mean and variance')
    M = s.Matrix([[-1,3],[1,-3]])
    pi, w, v = s.Matrix([R(3,4),R(1,4)]), s.Matrix([1,-3]), s.Matrix([R(1,4),-R(1,4)])
    require(M*pi == s.zeros(2,1) and M*v == -4*v and M.T*w == -4*w and (w.T*v)[0] == 1,
            'Two-state stationary, left, right and normalization identities')
    t, z = s.symbols('t z', real=True)
    trajectory = R(1,4)+(z-R(1,4))*s.exp(-4*t)
    require(s.simplify(s.diff(trajectory,t)-(1-4*trajectory)) == 0 and trajectory.subs(t,0) == z,
            'Two-state explicit trajectory solves ODE and initial condition')
    diff = [R(1,2)-pi[i] for i in range(2)]
    require(sum(abs(x) for x in diff)/2 == R(1,4) and sum(diff[i]**2/pi[i] for i in range(2)) == R(1,3),
            'TV and chi-square worked example')
    require(sum([R(199,200),R(56,5),R(249,4000)]) == R(49029,4000),
            'Primary first-column exit sum')
    require(R(1799,8)*q[2]/q[1] == R(17990,199), 'Primary reverse rate a32')
    A,B,C,tv,uv = R(249,10),R(251,25),R(28),R(2,5),R(1,160)
    require(A*tv-B == -R(2,25) and (A+B)/(C*(1+tv)) == R(1747,1960) and
            A*(1+tv*uv/2)-B*(1+uv/2)-C*(1-tv) == -R(7761,4000),
            'All three hand-calculated triplet 134 quantities')
    A,B,C,tv,uv = s.symbols('A B C t u')
    F=A*(1+tv*uv/2)-B*(1+uv/2)-C*(1-tv)
    require(s.factor(F-(2/(1+tv)+uv/2)*(A*tv-B)-(1-tv)/(1+tv)*(A+B-C*(1+tv))) == 0,
            'Exact branch identity F=hG+dH')
    a,b,tv,y=s.symbols('a b t y')
    f=(1-tv)/(1+tv+y)
    derivative=-a*tv*s.diff(f,tv)-b*y*s.diff(f,y)
    require(s.factor(derivative-(2*a*tv+b*y-(b-a)*tv*y)/(1+tv+y)**2) == 0,
            'Conditional first contrast derivative in quarter bound')
    Ph,Pc,mh,mc,w4=s.symbols('Ph Pc mh mc w4')
    bh=Ph*mh+(1-Ph)*w4;bc=Pc*mc+(1-Pc)*w4
    require(s.factor(bc-(Pc/Ph)*bh-Pc*(mc-mh)-(1-Pc/Ph)*w4) == 0,
            'Conditional-average cancellation identity')
    require(R('0.0938488356')-R('0.0937454280') == R('0.0001034076'),
            'Printed primary overlap interval margin')
    require(1-R(3,4)**2 == R(7,16) and 1-R(21,100)/R(23,100) == R(2,23),
            'Percentage and interval exercise conversions')
    M2=s.Matrix([[-R(1,2),2],[R(1,2),-2]])
    require(M2*s.Matrix([R(4,5),R(1,5)]) == s.zeros(2,1) and
            M2*s.Matrix([1,-1]) == -R(5,2)*s.Matrix([1,-1]),
            'Exercise 3 generator and spectral solution')
    x=s.symbols('x')
    require(s.series((1+x)*s.log(1+x),x,0,3).removeO() == x+x*x/2,
            'KL local expansion through second order')
    require(R(4,10)+R(3,10)-R(2,10)*R(7,2) == 0 and
            R(4,10)+R(3,10)+R(2,10)*R(7,2)**2 == R(63,20),
            'Exercise 11 prefix mean and squared length')
    require(s.Matrix([1,2])*s.Matrix([[3,4]])/3 == s.Matrix([[1,R(4,3)],[2,R(8,3)]]),
            'Exercise 12 Schur outer-product correction')
    clock_M=s.Matrix([[-2,1],[2,-1]])
    require(s.eye(2)+clock_M/2 == s.Matrix([[0,R(1,2)],[1,R(1,2)]]),
            'Lesson 34 uniformized two-state transition matrix')
    first=R(1,3)+R(2,3)*s.exp(-3*t)
    require(s.simplify(s.diff(first,t)-(1-3*first))==0 and first.subs(t,0)==1,
            'Lesson 34 two-state evolution and initial condition')
    require(s.simplify(s.exp(-1)/6/(1-R(1,4))-2*s.exp(-1)/9)==0,
            'Lesson 35 Poisson tail example')
    require(R(9,10)*R(9,1000)-R(8,1000)==R(1,10000),
            'Lesson 36 conservative cell margin')
    u,v,b=s.symbols('u v b',positive=True)
    trace=s.diag(-2*u,-2*v)-s.Matrix([b,b])*s.Matrix([[2*u,2*v]])/(-2*b)
    require(trace==s.Matrix([[-u,v],[u,-v]]),
            'Lesson 39 exact trace generator including diagonal entries')
    result={'passed':True,'checks':checks,'count':len(checks),
            'source_sha256':digest(ROOT/'research_report.tex'),
            'course_sha256':digest(ROOT/'research/foundational_course.tex'),
            'reviewed_technical_sha256':REVIEWED_TECHNICAL,
            'scope':'40-lesson insertion, original technical preservation before explicitly reviewed later edits, pinned current technical body, source-bound statements and worked algebra. Does not replace advanced proof review or visual PDF inspection.'}
    (OUT/'course_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'checks':len(checks)},indent=2))


if __name__ == '__main__':
    main()
