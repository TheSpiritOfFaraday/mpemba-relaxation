"""Check literal-data correspondence and the bounded focused-paper correction.

This does not replace mathematical proof review. The old focused-paper hashes
remain fixed: the only permitted change is the missing backslash in equation 15.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import re
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
OLD_SOURCE='f49ed22fa5394bdf1f1bb455dbdc8ed1b0c5643f7c334b17ef569c025293bb1e'
OLD_PDF='d32fbbc7aaf947f764dccc8ae1cfe204568812178ba06bbf31e43f4d83105c75'

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def require(test,label):
    if not test: raise ValueError(label)

def verify_typography_revision():
    source=(ROOT/'focused_paper.tex').read_text()
    corrected=r'\text{ordinary asymptotic cooling exists}\quad\Longleftrightarrow\quad'
    original=r'\text{ordinary asymptotic cooling exists}quad\Longleftrightarrow\quad'
    require(source.count(corrected)==1, 'exactly one corrected onset formula')
    reconstructed=source.replace(corrected,original)
    require(hashlib.sha256(reconstructed.encode()).hexdigest()==OLD_SOURCE,
            'source differs from original only by the intended one-character correction')
    before=ROOT/'research/detailed_audit/typography_revision/focused_paper_before.pdf'
    after=ROOT/'output/pdf/focused_paper.pdf'
    require(digest(before)==OLD_PDF,'original PDF preserved byte-for-byte')
    old,new=PdfReader(before),PdfReader(after)
    require(len(old.pages)==len(new.pages)==10,'ten pages preserved')
    changed=[]
    for i,(a,b) in enumerate(zip(old.pages,new.pages),1):
        aa,bb=a.extract_text(),b.extract_text()
        if aa!=bb:
            changed.append(i)
            require(aa.count('cooling exists quad')==1,'expected literal typesetting defect')
            require(aa.replace('cooling exists quad','cooling exists')==bb,
                    'extracted text changes only by removal of the stray word quad')
    require(changed==[6], 'only page six text changed')
    return {'source_before_sha256':OLD_SOURCE,'pdf_before_sha256':OLD_PDF,
            'source_after_sha256':digest(ROOT/'focused_paper.tex'),
            'pdf_after_sha256':digest(after),'changed_pages':changed,
            'change':'Restore one missing LaTeX backslash in equation 15.',
            'visual_check':'This function checks extracted text only. Hash-bound visual inspection is recorded separately in browser_and_pdf_review.json.'}

def main():
    arithmetic=json.loads((ROOT/'research/detailed_audit/freshman_arithmetic.json').read_text())
    lean=(ROOT/'formal/Mpemba/Data.lean').read_text()
    block=re.search(r'def M .*?!!\[(.*?)\]',lean,re.S).group(1)
    matrix=[[Fraction(x.strip()) for x in row.split(',')] for row in block.split(';')]
    expected=[[Fraction(x) for x in row] for row in arithmetic['generator']]
    require(matrix==expected,'all 16 Lean matrix literals equal reconstructed focused generator')
    report=(ROOT/'research_report.tex').read_text()
    appendix=(ROOT/'focused_appendices.tex').read_text()
    fragments=['physical_finite_time_goal','physical_finite_time_results',
               'physical_robustness_certificate','thermal_collision_realization','activated_intermediate_study']
    for name in fragments:
        fragment=(ROOT/f'research/{name}.tex').read_text()
        require(fragment in report,f'{name} included in report')
        require(fragment in appendix,f'{name} included in appendices')
    result={'passed':True,'scope':'Literal-data and bounded publication consistency only.',
            'lean_matrix_entries_compared':16,'fragments_compared_in_both_documents':fragments,
            'typography_revision':verify_typography_revision()}
    (ROOT/'research/detailed_audit/publication_consistency.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
