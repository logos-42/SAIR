import sys
from cypari2 import Pari
pari = Pari()

# Read first PSL(2,23) polynomial from baseline (ascending powers)
coeffs = [int(x) for x in open('/Users/apple/SAIR/igp24/data/psl223_r0_lmfdb.txt').readline().strip().split(',')]
f = pari.Polrev(coeffs)  # ascending powers -> monic degree 24
print("f =", f)
print("degree:", pari.poldegree(f))
print("irreducible:", pari.polisirreducible(f))
print("nfdisc abs:", abs(pari.nfdisc(f)))
print("real roots (polsturm):", pari.polsturm(f))

# Galois group via galoisinit
import time
t0 = time.time()
g = pari.galoisinit(f)
t1 = time.time()
print("galoisinit time: %.1fs" % (t1 - t0))
print("group order:", g.group.order())
print("num generators:", len(g.group.gen))
# export generators as permutations on roots (1-based in PARI)
perms = []
for gen in g.group.gen:
    # gen is a vector of images
    p = [int(x) for x in gen]
    perms.append(p)
print("generator perms (0-based):", perms)
print("group is solvable?", pari.galoisisabelian(f) if False else "n/a")
# galoisidentify only up to degree 11; use order + transitivity instead
print("frobenius info available; order =", g.group.order())
