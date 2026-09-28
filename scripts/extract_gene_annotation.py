#!/usr/bin/env python3
"""Create the three-column annotation table from GENCODE gene records."""
from pathlib import Path
import re
root = Path(__file__).resolve().parents[1]
source = root / 'reference/gencode.v48.primary_assembly.annotation.gtf'
out = root / 'metadata/gencode_v48_gene_annotation.tsv'
rows = {}
with source.open() as stream:
    for line in stream:
        if line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) != 9 or cols[2] != 'gene':
            continue
        attrs = dict(re.findall(r'(\w+) "([^"]*)"', cols[8]))
        gene = attrs['gene_id']
        row = (gene, attrs.get('gene_name', gene), attrs['gene_type'])
        if gene in rows and rows[gene] != row:
            raise ValueError(f'Conflicting annotation for {gene}')
        rows[gene] = row
if not rows:
    raise ValueError('No gene annotations found')
# Keep supplied annotation when equivalent; fail on unexpected differences.
if out.exists():
    supplied = {tuple(x.split('\t')) for x in out.read_text().splitlines()}
    if supplied != set(rows.values()):
        raise ValueError('Downloaded annotation differs from supplied table; investigate before replacing it')
    print('Supplied annotation matches downloaded GTF')
else:
    out.write_text(''.join('\t'.join(row) + '\n' for row in rows.values()))
