# SAIR

> Tools for the **IGP24** competition hosted by the SAIR Foundation.

This repository contains a Python client for the IGP24 competition: a
candidate generator that produces degree-24, monic, irreducible polynomials,
plus a submission client that verifies them against the SAIR public API.

**中文文档见 [`README.zh-CN.md`](./README.zh-CN.md) · [中文 README](README.zh-CN.md)**

---

## What is IGP24?

IGP24 is a competition about the inverse Galois problem over **Q** for degree
24. The goal is to realize Galois groups as splitting fields of monic
polynomials of degree 24 with integer coefficients. Polynomials are submitted
as **25 comma-separated coefficients** (ascending powers, i.e. constant term
first, leading coefficient last, monic), one polynomial per line, with an
optional trailing `# comment`.

## Repository layout

```
igp24/
├── scripts/
│   ├── gen_explore.py   # candidate polynomial generator (batch 1: exploration)
│   └── submit.py        # submission client for the SAIR public API
└── data/
    ├── lmfdb_baseline.csv  # reference baseline pulled from LMFDB
    ├── labels_progress.json# per-label discovery progress from the API
    ├── remaining_pairs.json# label/r pairs still to be found
    ├── test_batch.txt       # sanity-check batch (docs examples)
    └── explore_batch1.txt   # generated candidate batch
```

## Dependencies

- Python 3.8+
- [`numpy`](https://numpy.org/) — real-root counting via `np.roots`
- [`sympy`](https://www.sympy.org/) — irreducibility checks, factoring,
  resultants

Install with:

```bash
pip install numpy sympy
```

## Usage

### 1. Generate candidates

```bash
python3 igp24/scripts/gen_explore.py > igp24/data/explore_batch1.txt
```

The generator searches several families, all filtered to be **monic, with
non-zero constant term, and irreducible over Q**:

| Family  | Description                                             |
|---------|---------------------------------------------------------|
| `tri`   | Trinomials `x^24 + a·x^k + b`                            |
| `subst` | `g(x^k)` where `g` is a random irreducible deg 24/k poly |
| `rand`  | Random dense monic degree-24, small coefficients        |
| `tr`    | Totally-real search (deg 2,3,4,6,8), all-real-roots test|
| `comp`  | Composita (sum field) of totally-real pairs via resultant |
| `cheb`  | Chebyshev `T24(x) - c`                                  |
| `cyc`   | `x^24 ± a` cyclic candidates                            |

Each output line carries a trailing comment with its family and its real-root
count, e.g. `... # tri a=1 k=4 b=2 r=2`.

### 2. Submit a batch

```bash
export SAIR_API_KEY=<your-key>
python3 igp24/scripts/submit.py igp24/data/test_batch.txt          # quickly check a few
python3 igp24/scripts/submit.py igp24/data/explore_batch1.txt --poll  # submit + wait for verification
```

Key options:

| Flag     | Meaning                                              |
|----------|------------------------------------------------------|
| `--slot` | Optional submission slot number                      |
| `--dry`  | Print the payload size without sending               |
| `--poll` | Keep polling the batch until verification completes  |

The submitted polynomials must have **exactly 25 coefficients**; any line
without 25 is skipped with a warning. `submit.py` reads the key from the
`SAIR_API_KEY` environment variable and never stores it in a file.

## Data notes

- `lmfdb_baseline.csv` enumerates reference discriminants (with per-label
  `r` and polynomial discriminants) used as the scoring baseline.
- `labels_progress.json` / `remaining_pairs.json` are snapshots pulled from the
  public API tracking which `label × r` signatures are still undiscovered.

## License

Released under the [MIT License](./LICENSE).
