"""Evidence audit for the preserved finite-time research objective.

Run after the complete coverage replay and document builds. This checks
receipts and artifact consistency; it does not replace their mathematical
proofs, the independent exponential checks, or visual review.
"""
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import sys
from pypdf import PdfReader
from certify_physical_nominal import physical_box
from thermal_collision_model import algebra
from check_publication_consistency import verify_typography_revision

ROOT=Path(__file__).resolve().parents[1]
def read(path): return json.loads((ROOT/path).read_text())
def digest(path): return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

def main():
    parameters=read('research/physical_rate_refined.json')['result']['parameters']
    exact=json.loads((ROOT/'research/physical_rate_refined.json').read_text(),parse_float=Decimal)['result']['parameters']
    cert=read('research/robust_finite_time_certificate/certificate.json')
    assert cert['complete'] is True and cert['parameters']==parameters
    assert cert['input_sha256']==digest('research/physical_rate_refined.json')
    assert cert['runner_sha256']==digest('research/run_robust_finite_time_certificate.py')
    assert cert['relative_radius']==.05
    assert set(cert['obligations'])==set(map(str,range(1,9)))
    for name,build in cert['builds'].items():
        path='research/'+name+'.cpp'
        assert digest(path)==build['source_sha256']
        assert hashlib.sha256(Path(build['command'][-1]).read_bytes()).hexdigest()==build['binary_sha256']
        embedded=re.search(r'const char\*params\[\]=\{([^}]+)\}',(ROOT/path).read_text()).group(1)
        assert [Decimal(x) for x in re.findall(r'"([^"]+)"',embedded)]==exact
    for mode,record in cert['obligations'].items():
        assert record['complete'] is True and record['pending']==0 and record['enclosure_errors']==0
        assert record['mode']==int(mode) and record['radius']>=.05
        assert record['start_u'][0]<=.1<=record['start_u'][1] and record['end_u']==[2,2]
        log=next((ROOT/'research/robust_finite_time_certificate').glob(mode+'_*.log'))
        assert hashlib.sha256(log.read_bytes()).hexdigest()==record['log_sha256']
    ranges=physical_box(parameters)
    assert ranges==read('research/physical_box_ranges_recheck.json')
    collision=read('research/thermal_collision_audit.json')
    assert collision['parameters']==parameters and algebra(parameters)==collision['algebra']
    assert collision['event_audit']['all_diagnostic_checks_passed'] is True
    assert collision['event_audit']['trajectories_per_preparation']==500000
    independent=read('research/temperature_reduction_independent_audit.json')
    assert independent['all_passed'] is True
    assert independent['temperature_derivative_value_checks']==6744 and independent['gradient_checks']==360
    assert all(read('research/threshold_constraint_logic_audit.json').values())
    nominal=read('research/physical_rate_refined_interval_certificate.json')
    assert nominal['parameters']==parameters and nominal['certified_min_time_saving_lower']>.0812
    assert nominal['cold_distance_floor_lower']>.001 and nominal['minimum_KL_inequality_margin']>0
    point=read('research/physical_rate_refined_audit.json')
    assert point['parameters']==parameters and point['independent']['all_original_triplet_tests_fail'] is False
    assert sum(not x['failure'] for x in point['independent']['triplets'])==2
    assert point['robustness']['failing_cases']==0
    rows=list(csv.DictReader((ROOT/'output/figures/physical_finite_time.csv').open()))
    assert len(rows)==151 and all(k in rows[0] for k in ['TV_cold','TV_hot','chi2_cold','chi2_hot'])
    assert float(rows[0]['KL_hot'])>float(rows[0]['KL_cold'])>.01
    expected={
      'focused_paper.tex':'f49ed22fa5394bdf1f1bb455dbdc8ed1b0c5643f7c334b17ef569c025293bb1e',
      'output/pdf/focused_paper.pdf':'d32fbbc7aaf947f764dccc8ae1cfe204568812178ba06bbf31e43f4d83105c75',
      'research/archive/research_report_before_focused_paper.tex':'346bceb86c912d45af05a5dbcb0216510988adb140306ca873bac5c701a29569',
      'research/archive/research_report_before_focused_paper.pdf':'d8052d5049348c30197b3d2092cc8a5c895bf2326ee1d53dbdb4a30d2617719c'}
    # Retain the old hashes as the baseline; permit only the independently
    # checked one-character equation-15 typography correction and its PDF.
    focused_revision = None
    if digest('focused_paper.tex') != expected['focused_paper.tex']:
        focused_revision = verify_typography_revision()
        expected['focused_paper.tex'] = focused_revision['source_after_sha256']
        expected['output/pdf/focused_paper.pdf'] = focused_revision['pdf_after_sha256']
    for path,value in expected.items(): assert digest(path)==value,path
    pages={}
    for name in ['focused_paper','research_report','focused_appendices']:
        pdf=PdfReader(ROOT/'output/pdf'/f'{name}.pdf');pages[name]=len(pdf.pages)
        text='\n'.join(p.extract_text() for p in pdf.pages)
        assert '??' not in text
        if name!='focused_paper':
            assert '135,721' in text and '65,045' in text and 'energy-conserving' in text
            assert '5% time-saving claim remains open' not in text
    assert 8<=pages['focused_paper']<=12
    source=(ROOT/'research_report.tex').read_text()
    assert 'Comprehensive research report' in source
    for fragment in ['physical_finite_time_goal','physical_finite_time_results','physical_robustness_certificate','thermal_collision_realization','activated_intermediate_study']:
        assert (ROOT/f'research/{fragment}.tex').read_text() in source
    # Explicitly maintain physical and statistical scope in the delivered text.
    for phrase in ['idealized stochastic collision model','not a Lean formalization or external peer review',
                   'not a certificate','all-triplet-failure premise','preparation costs']:
        assert phrase in source,phrase
    requirements=[
      {'requirement':'Physically motivated reversible model with shared geometry',
       'evidence':['research/physical_model_design.txt','research/thermal_collision_realization.tex','research/thermal_collision_audit.json'],
       'finding':'Fixed two-basin geometry, symmetric saddles, and exact energy-conserving canonical collision realization; idealized stochastic assumptions explicit.'},
      {'requirement':'Moderate energies, barriers, rates and temperatures across the robustness domain',
       'evidence':['research/physical_box_ranges_recheck.json'], 'finding':ranges},
      {'requirement':'Finite-time distance advantage at resolved distances',
       'evidence':['research/robust_finite_time_certificate/certificate.json'],
       'finding':'Mode 1 covers full +/-5% box and u in [0.4,0.8]: Dh<=0.9Dc, Dc>=0.001. Interval length0.4tau exceeds0.25tau.'},
      {'requirement':'Time savings at every threshold and all simultaneous +/-5% perturbations',
       'evidence':['research/physical_robustness_certificate.tex','research/robust_finite_time_certificate/certificate.json'],
       'finding':'Modes2,3,6 bracket all roots; mode4 proves cold concavity; mode5 proves hot-root monotonicity; modes7,8 prove both sufficient endpoint cases. Analytic reduction yields th<=0.95tc for every threshold in[.001,.01].'},
      {'requirement':'Independent verification, explicit numerical error and reproducibility',
       'evidence':['research/temperature_reduction_independent_audit.json','research/threshold_constraint_logic_audit.json','research/run_robust_finite_time_certificate.py','RUNME.txt'],
       'finding':'Fresh source-bound coverage replay, 70-digit direct-exponential and derivative cross-check, false-claim controls, and separate nominal uniformization certificate. Reviewable interval proof, not Lean or external peer review.'},
      {'requirement':'Report methods, metrics, mechanism, physical limits and search history',
       'evidence':['research_report.tex','output/pdf/research_report.pdf','output/figures/physical_finite_time.csv','research/activated_intermediate_study.tex'],
       'finding':'KL distance versus threshold time, absolute and normalized clocks, alternative-distance diagnostics, rejected candidates, exact stochastic bath and separate coarse-graining sensitivity recorded. No material, experimental, global-optimum or new-triplet-counterexample claim.'},
      {'requirement':'Preserve focused paper and archived original; retain appendices',
       'evidence':list(expected)+['output/pdf/focused_appendices.pdf'],
       'finding':'Preserved baseline hashes or the verified equation-15 typography revision match; focused paper remains10 pages. Updated report and appendices compile and were rendered for separate visual inspection.'}]
    evidence=set(expected)
    for item in requirements:evidence.update(item['evidence'])
    result={'requirements_checked':True,'pages':pages,'requirements':requirements,
            'focused_typography_revision':focused_revision,
            'artifact_sha256':{p:digest(p) for p in sorted(evidence)},
            'scope':'Research goal achieved for an idealized reversible stochastic rate model at the predeclared quantitative benchmarks. No experimentally demonstrated cooling, preparation-cost advantage, universal metric claim, or material validation.',
            'visual_review':'Recorded separately by the authoring agent; this script checks text/artifact consistency, not page pixels.'}
    (ROOT/'research/physical_finite_time_completion_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'requirements_checked':True,'pages':pages,'proof_modes':8,'independent_value_checks':6744,'independent_gradient_checks':360},indent=2))

if __name__=='__main__':main()
