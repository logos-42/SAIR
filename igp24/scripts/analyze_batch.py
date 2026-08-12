#!/usr/bin/env python3
"""Analyze a verified batch against current targets. Usage:
  python3 analyze_batch.py <batch_status.json> <submitted_polys.txt> [--all]
Prints: verification stats, distinct (t,r), NEW non-baseline hits (vs new_targets.json),
baseline-pair disc-improvement candidates, and per-family breakdown.
"""
import json, sys, csv
from collections import Counter

def load_new_targets():
    try:
        return set((it['t'], it['r']) for it in json.load(open('data/new_targets.json')))
    except Exception:
        return None

def load_baseline():
    base = set()
    try:
        with open('data/lmfdb_baseline.csv') as f:
            for row in csv.DictReader(f):
                base.add((int(row['label'][3:]), int(row['r'])))
    except Exception:
        pass
    return base

def main():
    status_path, poly_path = sys.argv[1], sys.argv[2]
    d = json.load(open(status_path))['data']
    polys = d.get('polynomials', [])
    ok = [p for p in polys if p.get('status') == 'ok']
    lines = [l.strip() for l in open(poly_path) if l.strip()]
    print(f"batch {d.get('batchStatus')}: ok={d.get('verifiedCount')} failed={d.get('failedCount')} total={d.get('totalCount')}")
    tc = Counter((p['label'], p['r']) for p in ok)
    print(f"distinct (t,r): {len(tc)}")
    print("top 15:", tc.most_common(15))
    newt = load_new_targets()
    base = load_baseline()
    if newt is not None:
        hits = [(p['t'], p['r'], p['label'], p.get('scoreable')) for p in ok
                if (p['t'], p['r']) in newt]
        print(f"\n*** NEW non-baseline hits: {len(set((t,r) for t,r,_,_ in hits))} ***")
        for t, r, label, sc in sorted(set(hits)):
            print(f"  HIT {label} r={r} scoreable={sc}")
    # scoreable baseline-disc candidates
    disc = [(p['t'], p['r'], p['label'], p.get('fieldDiscAbs'), p.get('baselineDiscAbs'))
            for p in ok if p.get('inBaseline') and p.get('scoringStatus') in ('scoreable', 'pending')
            and p.get('baselineUnlocked')]
    if disc:
        print(f"\nbaseline unlocked candidates: {len(disc)}")
        for t, r, label, fd, bd in disc[:10]:
            print(f"  {label} r={r} fieldDisc={fd} base={bd}")
    # family breakdown (from tags)
    fam_label = Counter()
    for p in ok:
        idx = p['polynomialIndex']
        if idx >= len(lines): continue
        tag = lines[idx].split('#')[1].strip() if '#' in lines[idx] else '?'
        fam = tag.split()[0]
        fam_label[(fam, p['label'])] += 1
    print("\nfamily -> label (top):")
    for (fam, label), n in fam_label.most_common(25):
        print(f"  {fam:10s} {label} x{n}")

if __name__ == "__main__":
    main()
