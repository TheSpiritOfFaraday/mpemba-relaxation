"""Build the offline guide and a source-indexed audit scope; bind completion to reviewed hashes."""
from pathlib import Path
import datetime
import hashlib
import html
import json
import re

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/detailed_audit'

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    OUT.mkdir(exist_ok=True)
    inventory=[]
    artifacts=[]
    for path in sorted(ROOT.rglob('*')):
        rel=path.relative_to(ROOT)
        # These are generated after this inventory; delivery receipts bind
        # their final hashes without creating a self-referential snapshot.
        if str(rel) in {'freshman_guide.html', 'output/pdf/freshman_guide.pdf'}: continue
        if any(x in rel.parts for x in ['.lake','.tools','__pycache__','tmp','.git']): continue
        if 'detailed_audit' in rel.parts: continue # own receipts are linked separately
        if path.is_file() and path.suffix in {'.tex','.py','.cpp','.lean','.sh','.md','.txt','.html'}:
            inventory.append({'path':str(rel),'sha256':digest(path),
                              'role':'historical' if 'archive' in rel.parts else 'current-or-supporting'})
        elif path.is_file() and path.suffix in {'.json','.csv','.pdf','.png','.log'}:
            artifacts.append({'path':str(rel),'sha256':digest(path),
                              'role':'historical' if 'archive' in rel.parts else 'data-output-or-receipt; interpretation requires review'})
    claims=[]
    reviews_path=OUT/'claim_reviews.json'
    reviews=json.loads(reviews_path.read_text()) if reviews_path.exists() else {}
    for filename in ['focused_paper.tex','research_report.tex','focused_appendices.tex']:
        source=(ROOT/filename).read_text()
        for i,m in enumerate(re.finditer(r'\\begin\{(theorem|lemma|corollary|proposition)\}(?:\[([^\]]+)\])?', source),1):
            end=source.find(r'\end{'+m[1]+'}',m.end())
            block=source[m.start():end+len(m[1])+6]
            label=re.search(r'\\label\{([^}]+)\}',block)
            claim={'id':filename+':'+str(i),'file':filename,
                           'line':source[:m.start()].count('\n')+1,'kind':m[1],
                           'title':m[2] or 'Untitled '+m[1], 'label':label[1] if label else None,
                           'statement_sha256':hashlib.sha256(block.encode()).hexdigest(),
                           'review_status':'pending detailed proof audit'}
            review=reviews.get(claim['id'])
            if review and review['source_sha256']==digest(ROOT/filename) and review['statement_sha256']==claim['statement_sha256']:
                claim.update(review_status=review['status'],review_evidence=review['evidence'])
            claims.append(claim)
    commands=[]
    for i,line in enumerate((ROOT/'RUNME.txt').read_text().splitlines(),1):
        if re.match(r'^  (?:python3|OPENBLAS_NUM_THREADS|sh |clang\+\+|\./|lake |cd )',line):
            commands.append({'line':i,'command':line.strip(),'coverage':'not automatically inferred from a passing suite'})
    receipts=[]
    for path in sorted(OUT.glob('*/receipt.json')):
        item=json.loads(path.read_text()); item['receipt_path']=str(path.relative_to(ROOT)); receipts.append(item)
    completion_path=OUT/'final_completion.json'
    completion=json.loads(completion_path.read_text()) if completion_path.exists() else {}
    expected=completion.get('source_hashes', {})
    current={item['path']:item['sha256'] for item in inventory}
    completed=(completion.get('complete') is True and current==expected and
               all(c['review_status'].startswith('reviewed:') for c in claims))
    scope={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
           'objective':'Double check all work step by step in extreme detail and create a college-freshman HTML guide.',
           'overall_complete':completed,
           'completion_evidence':'research/detailed_audit/final_completion.json',
           'completion_scope':'Current report claims, derivations, advertised computations and teaching deliverables; exploratory archives are assessed as provenance, not certified as alternative proof software.',
           'inventory':inventory,'data_output_inventory':artifacts,
           'inventory_excludes':['Downloaded toolchains and dependency caches','Temporary directory','This audit directory and generated HTML/PDF guide; final hashes recorded in delivery receipts'],
           'numbered_claims':claims,'reproduction_commands':commands,
           'additional_required_reviews':[
               'Source correspondence and literature attribution, including original predicate branches',
               'Unnumbered arguments, computations, figures, tables, limits and scope statements',
               'Eight continuous physical-box obligations and all interval implementation dependencies',
               'Negative controls and failure/unresolved handling for exact classifiers and interval engines',
               'Historical searches distinguished from proofs; no inference from failure to find an example',
               'Generated report/appendix consistency, PDF layout and beginner-guide browser verification',
               'Fresh replay of remaining advertised verification commands and dependency audit'],
           'replay_receipts':receipts}
    (OUT/'scope.json').write_text(json.dumps(scope,indent=2)+'\n')
    rows=[]
    for item in receipts:
        for c in item['checks']:
            rows.append('<tr><td>'+html.escape(c['script'])+'</td><td>'+html.escape(c['status'])+
                        '</td><td><a href="'+html.escape(c['log'])+'">Run log</a></td></tr>')
    claims_html=''.join('<tr><td>'+html.escape(c['file'])+'<br><span class="small">line '+str(c['line'])+'</span></td><td>'+html.escape(c['title'])+
                        '</td><td>'+html.escape(c['review_status'])+'</td></tr>' for c in claims)
    arithmetic=json.loads((OUT/'freshman_arithmetic.json').read_text())
    template=(ROOT/'research/freshman_guide_template.html').read_text()
    for name in ['mechanism', 'contrast', 'shape', 'sharpness', 'broader', 'algorithms', 'entropy', 'dimensions', 'intervals', 'physical_support', 'sources_asymptotics', 'figures', 'discovery']:
        marker='<!--DETAIL_'+name.upper()+'-->'
        if template.count(marker)!=1:
            raise ValueError('Expected exactly one detail marker: '+marker)
        template=template.replace(marker, (ROOT/('research/guide_'+name+'_detail.html')).read_text())
    guide=template.replace('<!--REPLAY_ROWS-->', ''.join(rows)).replace('<!--CLAIM_ROWS-->',claims_html)
    guide=guide.replace('<!--AUDIT_STAMP-->',html.escape(scope['generated_utc']))
    guide=guide.replace('/*ARITHMETIC_DATA*/',json.dumps(arithmetic).replace('</','<\\/'))
    (ROOT/'freshman_guide.html').write_text(guide)
    print(f'Built freshman_guide.html; inventoried {len(inventory)} source/support files, {len(claims)} numbered claims, {len(commands)} command lines.')

if __name__=='__main__': main()
