#!/usr/bin/env python3
"""Regenerate STAR summary tables from original STAR logs and GeneCounts."""
from pathlib import Path
import csv
root = Path(__file__).resolve().parents[1]
runs = (root/'metadata/SRR_accessions.txt').read_text().split()
qc, strands = [], []
keys = ['Number of input reads','Uniquely mapped reads %','% of reads mapped to multiple loci',
        '% of reads mapped to too many loci','% of reads unmapped: too short','% of reads unmapped: other']
for run in sorted(runs):
    prefix = root/'aligned'/run/run
    metrics = {}
    for line in Path(str(prefix)+'_Log.final.out').read_text().splitlines():
        if '|' in line:
            k, v = line.split('|', 1)
            metrics[k.strip()] = v.strip()
    qc.append([run] + [metrics[k] for k in keys])
    totals = [0, 0, 0]
    for line in Path(str(prefix)+'_ReadsPerGene.out.tab').read_text().splitlines():
        cols = line.split('\t')
        if cols[0].startswith('N_'):
            continue
        totals = [a + int(b) for a, b in zip(totals, cols[1:4])]
    strands.append([run] + totals + [f'{totals[2]/totals[1]:.2f}' if totals[1] else 'NA'])
(root/'results').mkdir(exist_ok=True)
for name, header, rows in [
    ('STAR_alignment_QC.tsv', ['sample','input_reads','uniquely_mapped_pct','multi_mapped_pct','too_many_loci_pct','unmapped_too_short_pct','unmapped_other_pct'], qc),
    ('STAR_strandedness_check.tsv', ['sample','unstranded_counts','forward_stranded_counts','reverse_stranded_counts','reverse_to_forward_ratio'], strands)]:
    with (root/'results'/name).open('w') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)
