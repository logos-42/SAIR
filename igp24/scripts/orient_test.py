from cypari2 import Pari
pari = Pari()

coeffs = [int(x) for x in open('/Users/apple/SAIR/igp24/data/psl223_r0_lmfdb.txt').readline().strip().split(',')]
csv_nfdisc = 36245192713290003306351501713852497190281934538751108553916766514130714624

for name, cl in [("as-given (a24=1)", coeffs),
                 ("reversed (a0=1)", coeffs[::-1])]:
    f = pari.Pol(cl)
    irr = pari.polisirreducible(f)
    nd = abs(pari.nfdisc(f))
    print(name)
    print("  irreducible:", irr)
    print("  nfdisc:", nd)
    print("  nfdisc matches CSV:", nd == csv_nfdisc)
    print("  polsturm:", pari.polsturm(f))
    print()
