#!/usr/bin/env python3
"""Rebuild Markdown catalog and source index without network access.

Run before rebuilding the package manifest. Generated navigation indexes and
integrity reports are excluded from the catalog to avoid circular inventories.
"""
from pathlib import Path
from urllib.parse import quote
import json
import math
import os

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / 'docs/catalog'
PAGE_SIZE = 120
GROUPS = [
    ('reading-and-release-guides', 'Reading and release guides', 'Manuscripts, documentation, reproduction instructions, and release tools.'),
    ('research-programs', 'Research programs', 'Scientific checks, discovery programs, and plotting scripts outside the detailed audit snapshots.'),
    ('mathematical-fragments', 'Mathematical fragments', 'LaTeX statements and arguments incorporated into the manuscripts.'),
    ('research-results-and-notes', 'Research results and notes', 'Saved certificates, numerical outputs, logs, and notes outside detailed audit snapshots.'),
    ('lean-formal-proofs', 'Lean formal proofs', 'Proof modules, theorem audits, and pinned dependency definitions.'),
    ('papers-and-guides', 'Papers and guides', 'Exported PDFs and their historical compilation logs.'),
    ('figures-and-data', 'Figures and data', 'Figure images, vector PDFs, tables, and plotting receipts.'),
    ('research-audit-records', 'Research audit records', 'Historical reviews, replays, and isolated snapshots; these are provenance, not a recommended first reading route.'),
    ('historical-discovery-archive', 'Historical discovery archive', 'Earlier discovery files preserved for provenance.'),
    ('guide-screenshots', 'Guide screenshots', 'Historical visual checks of the interactive teaching guide.'),
]
REPORTS = {'PACKAGE_MANIFEST.json', 'PRIVACY_REPORT.json', 'PACKAGE_CHECKS.json', 'NAVIGATION_CHECKS.json', 'ANONYMIZATION.json'}


def group(rel):
    prefix = rel.parts[:2]
    mapping = {('research', 'detailed_audit'): 'research-audit-records',
               ('research', 'archive'): 'historical-discovery-archive',
               ('output', 'playwright'): 'guide-screenshots',
               ('output', 'pdf'): 'papers-and-guides',
               ('output', 'figures'): 'figures-and-data'}
    if prefix in mapping:
        return mapping[prefix]
    if rel.parts[0] == 'formal':
        return 'lean-formal-proofs'
    if rel.parts[0] == 'research':
        if rel.suffix in {'.py', '.cpp', '.sh'}:
            return 'research-programs'
        if rel.suffix == '.tex':
            return 'mathematical-fragments'
        return 'research-results-and-notes'
    return 'reading-and-release-guides'


def link(origin, target):
    return quote(os.path.relpath(target, origin.parent).replace(os.sep, '/'), safe='/.-_')


def size(n):
    return f'{n / 1048576:.1f} MiB' if n >= 1048576 else f'{max(1, math.ceil(n / 1024))} KiB'


def main():
    CATALOG.mkdir(parents=True, exist_ok=True)
    groups = {key: [] for key, _, _ in GROUPS}
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file():
            continue
        rel = p.relative_to(ROOT)
        if any(x in {'.git', '__pycache__', '.venv', '.lake', '.tools'} for x in rel.parts):
            continue
        if rel.parts[:2] == ('docs', 'catalog') or rel.as_posix() == 'docs/REFERENCES.md':
            continue
        if len(rel.parts) == 1 and rel.name in REPORTS:
            continue
        groups[group(rel)].append(p)
    total = sum(map(len, groups.values()))
    overview = ['# Complete file catalog', '', '[Home](../../README.md) · [Repository guide](../README.md) · [Research map](../RESEARCH_MAP.md)', '',
                'Use the research map for a curated route. This catalog links to every research and supporting file, including the historical snapshots.', '',
                f'**{total:,} indexed files**, arranged by purpose. Each page has at most {PAGE_SIZE} entries. Generated navigation indexes and integrity reports are excluded; the [package manifest](../../PACKAGE_MANIFEST.json) covers every distributed file except itself.', '']
    generated = set()
    overview += ['| Collection | Files | Open |', '| --- | ---: | --- |']
    for key, title, _ in GROUPS:
        overview.append(f'| [{title}](#{key}) | {len(groups[key]):,} | [First page]({key}-01.md) |')
    overview.append('')
    for key, title, description in GROUPS:
        members = groups[key]
        pages = [members[i:i + PAGE_SIZE] for i in range(0, len(members), PAGE_SIZE)]
        overview += ['## ' + title, '', description, '', f'{len(members):,} files.', '']
        for index, chunk in enumerate(pages, 1):
            name = f'{key}-{index:02d}.md'
            generated.add(name)
            path = CATALOG / name
            previous = f'[Previous]({key}-{index-1:02d}.md)' if index > 1 else ''
            following = f'[Next]({key}-{index+1:02d}.md)' if index < len(pages) else ''
            nav = ' · '.join(x for x in ['[Home](../../README.md)', f'[Catalog](README.md#{key})', previous, following] if x)
            lines = [f'# {title} — {index} of {len(pages)}', '', nav, '', description, '', '| File | Type | Size |', '| --- | --- | --- |']
            for p in chunk:
                rel = p.relative_to(ROOT).as_posix()
                label = rel.replace('|', '\\|').replace('[', '\\[').replace(']', '\\]')
                lines.append(f'| [{label}]({link(path, p)}) | {p.suffix.lstrip(".") or "text"} | {size(p.stat().st_size)} |')
            lines += ['', nav, '']
            path.write_text('\n'.join(lines))
            overview.append(f'- [Page {index}]({name}) — {len(chunk)} files; `{chunk[0].name}` through `{chunk[-1].name}`')
        overview.append('')
    (CATALOG / 'README.md').write_text('\n'.join(overview))
    for p in CATALOG.glob('*.md'):
        if p.name != 'README.md' and p.name not in generated and any(p.name.startswith(key + '-') for key, _, _ in GROUPS):
            p.unlink()
    sources = json.loads((ROOT / 'THIRD_PARTY_SOURCES.json').read_text())['files']
    by_url = {}
    for source in sources:
        by_url.setdefault(source['source_url'], []).append(source['former_path'])
    refs = ['# Reference literature', '', '[Home](../README.md) · [Repository guide](README.md) · [Research map](RESEARCH_MAP.md)', '',
            'These links come from the preserved research records. The manuscripts give full scholarly citations and explain how each source is used. This index is not a new literature review.', '',
            'Downloaded third-party full text is omitted. Former local paths help interpret historical receipts; obtain papers from the cited publishers or repositories under their terms.', '',
            'The [machine-readable source inventory](../THIRD_PARTY_SOURCES.json) retains the reference-copy checksums.', '', '| Public source | Former local copies |', '| --- | --- |']
    for url, paths in sorted(by_url.items()):
        labels = '<br>'.join('`' + p + '`' for p in sorted(paths))
        refs.append(f'| [{url.removeprefix("https://")}]({url}) | {labels} |')
    refs += ['', 'See [rights and attribution](../RIGHTS.md) for reuse status.', '']
    (ROOT / 'docs/REFERENCES.md').write_text('\n'.join(refs))
    print(json.dumps({'indexed_files': total, 'catalog_pages': len(generated) + 1, 'reference_links': len(by_url)}, indent=2))


if __name__ == '__main__':
    main()
