import json, csv
from collections import Counter

rem = json.load(open('data/remaining_pairs_live.json'))['data']['items']
rems = set((it['t'], it['r']) for it in rem)
base = set()
with open('data/lmfdb_baseline.csv') as f:
    rdr = csv.DictReader(f)
    for row in rdr:
        base.add((int(row['label'][3:]), int(row['r'])))
print('baseline pairs:', len(base), '| remaining pairs:', len(rems))
new_targets = rems.difference(base)
disc_targets = rems.intersection(base)
print('TRUE new-discovery targets (non-baseline):', len(new_targets))
print('baseline pairs needing disc improvement:', len(disc_targets))
print('new targets by r:', dict(Counter(r for t, r in new_targets)))
def bucket(t):
    if t <= 100: return '1-100'
    if t <= 1000: return '101-1000'
    if t <= 5000: return '1001-5000'
    if t <= 10000: return '5001-10000'
    if t <= 20000: return '10001-20000'
    return '20001+'
print('new targets by t-range:', dict(Counter(bucket(t) for t, r in new_targets)))
small = sorted(new_targets)[:30]
for t, r in small:
    print('  new', t, 'r=', r)
json.dump([{'t': t, 'r': r} for t, r in sorted(new_targets)], open('data/new_targets.json', 'w'))
