#!/bin/bash
# helix2 totally-real batch: 8 seeds x 30 targets, output to explore_batch14.txt
cd /Users/apple/SAIR/igp24 || exit 1
rm -f data/explore_batch14.txt
for s in 141 142 143 144 145 146 147 148; do
  python3 -u scripts/gen_helix2.py "$s" 30 > "data/helix2_s${s}.txt" 2>> /tmp/h2d_err.log
  echo "seed $s: $(wc -l < data/helix2_s${s}.txt)"
done
cat data/helix2_s141.txt data/helix2_s142.txt data/helix2_s143.txt data/helix2_s144.txt \
    data/helix2_s145.txt data/helix2_s146.txt data/helix2_s147.txt data/helix2_s148.txt \
    > data/explore_batch14.txt
echo "TOTAL: $(wc -l < data/explore_batch14.txt)"
awk -F'r=' '{print $2}' data/explore_batch14.txt | sort -n | uniq -c
