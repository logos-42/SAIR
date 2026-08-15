#!/bin/bash
# cheb batch26: 8 seeds x 60 targets, all (k, base) configs
cd /Users/apple/Downloads/SAIR/igp24 || exit 1
rm -f data/explore_batch26.txt
for s in 1301 1302 1303 1304 1305 1306 1307 1308; do
  python3 -u scripts/gen_cheb.py "$s" 60 0 > "data/cheb_s${s}.txt" 2>> /tmp/cheb_err.log
  echo "seed $s: $(wc -l < data/cheb_s${s}.txt)"
done
cat data/cheb_s1301.txt data/cheb_s1302.txt data/cheb_s1303.txt data/cheb_s1304.txt \
    data/cheb_s1305.txt data/cheb_s1306.txt data/cheb_s1307.txt data/cheb_s1308.txt \
    > data/explore_batch26.txt
echo "TOTAL: $(wc -l < data/explore_batch26.txt)"
awk -F'r=' '{print $2}' data/explore_batch26.txt | sort -n | uniq -c
