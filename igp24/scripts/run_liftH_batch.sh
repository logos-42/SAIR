#!/bin/bash
# liftH batch25: quadray bases scan (8 (D,m) combos x 30 targets)
cd /Users/apple/Downloads/SAIR/igp24 || exit 1
rm -f data/explore_batch25.txt
s=1201
for dm in "-23 4" "-23 5" "-31 3" "-8 3" "-20 3" "-24 3" "-40 3" "-15 4"; do
  set -- $dm
  D=$1; M=$2
  python3 -u scripts/gen_liftH.py "$s" 30 full "$D" "$M" > "data/liftH_${D}_${M}.txt" 2>> /tmp/liftH_err.log
  echo "D=$D m=$M: $(wc -l < data/liftH_${D}_${M}.txt)"
  s=$((s + 1))
done
cat data/liftH_-23_4.txt data/liftH_-23_5.txt data/liftH_-31_3.txt data/liftH_-8_3.txt \
    data/liftH_-20_3.txt data/liftH_-24_3.txt data/liftH_-40_3.txt data/liftH_-15_4.txt \
    > data/explore_batch25.txt
echo "TOTAL: $(wc -l < data/explore_batch25.txt)"
awk -F'r=' '{print $2}' data/explore_batch25.txt | sort -n | uniq -c
