#!/usr/bin/env python3
"""Recon: how many open pairs are pure 2-groups / small order."""
import json
from collections import defaultdict

nt = json.load(open("data/new_targets.json"))
orders = {}
for l in open("data/group_orders.tsv"):
    t, o, s = l.strip().split("\t")
    orders[int(t)] = (int(o), s)

def is_pow2(o):
    while o > 1 and o % 2 == 0:
        o //= 2
    return o == 1

cnt_2group = set()
cnt_2group_pairs = 0
cnt_small = 0
for x in nt:
    t, r = x["t"], x["r"]
    o, solv = orders.get(t, (0, "?"))
    if o and is_pow2(o):
        cnt_2group.add(t)
        cnt_2group_pairs += 1
    if o and o <= 100000:
        cnt_small += 1
print(f"纯 2-群开放: {len(cnt_2group)} 群 / {cnt_2group_pairs} 对")
print(f"order<=1e5 的开放: {cnt_small} 对")
