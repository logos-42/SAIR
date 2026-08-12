from cypari2 import Pari
pari = Pari()

def analyze(path, nm):
    print(f"=== {nm} ===")
    for i, line in enumerate(open(path)):
        coeffs = [int(x) for x in line.strip().split(',')]
        f = pari.Pol(coeffs[::-1])  # ascending a0..a24, monic
        fp = pari.deriv(f)
        r_f = pari.polsturm(f)
        r_fp = pari.polsturm(fp)
        # real critical points count = r_fp if no repeated roots of f'
        disc_fp = pari.poldisc(fp)
        sq = pari.issquare(disc_fp)
        # interlacing check via signs: number of real roots of f(x)-c max = ?
        # max real roots of f(x)-c over c = 24 - 2*(#complex critical points with real critical value? ) ...
        # simpler: max over c in [-M,M] sampled: skip, use theoretical: if r_fp==23 then can reach 24 only if
        # all critical values real and max(minima)<min(maxima). Check distinctness + count.
        print(f"  poly{i}: r(f)={r_f}  r(f')={r_fp}  disc(f') square?={sq}")
    print()

analyze('/Users/apple/SAIR/igp24/data/psl223_r0_lmfdb.txt', "PSL(2,23) 24T7817")
analyze('/Users/apple/SAIR/igp24/data/pgl223_r0_lmfdb.txt', "PGL(2,23) 24T10255")
