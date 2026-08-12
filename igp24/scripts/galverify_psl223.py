from cypari2 import Pari
import time, json, sys
pari = Pari()

coeffs = [int(x) for x in open('/Users/apple/SAIR/igp24/data/psl223_r0_lmfdb.txt').readline().strip().split(',')]
f = pari.Pol(coeffs)  # ascending powers, monic, a24=1
print("f =", f, flush=True)
t0 = time.time()
g = pari.galoisinit(f)
t1 = time.time()
print("galoisinit OK in %.1fs" % (t1 - t0), flush=True)
print("group order:", g.group.order(), flush=True)
gens = []
for gen in g.group.gen:
    p = [int(x) for x in gen]  # 1-based images
    gens.append(p)
print("num gens:", len(gens), flush=True)
print("gens 1-based:", json.dumps(gens), flush=True)
# frob: check orders of Frobenius elements for a few primes (consistency check)
print("done", flush=True)
