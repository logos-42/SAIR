#!/usr/bin/env python3
"""IGP24 submission client.

Reads polynomials (25 comma-separated coeffs, ascending powers, monic) from a
file, one per line; optional trailing # comments allowed. POSTs to the SAIR
public API, prints per-polynomial verified results.

Usage:
  python3 submit.py polys.txt [--slot N] [--dry]
"""
import argparse, json, os, sys, time, urllib.request, urllib.error

API = "https://api.sair.foundation/api/public/v1/competitions/igp24/submissions"
KEY = os.environ.get("SAIR_API_KEY")
if not KEY:
    sys.exit("SAIR_API_KEY env var required")


def parse_polys(path):
    lines = []
    with open(path) as f:
        for ln in f:
            ln = ln.split("#", 1)[0].strip()
            if not ln:
                continue
            parts = [p.strip() for p in ln.split(",")]
            if len(parts) != 25:
                print(f"SKIP (not 25 coeffs): {ln[:60]}", file=sys.stderr)
                continue
            lines.append(",".join(parts))
    return lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--slot", type=int, default=None)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--poll", action="store_true", help="poll batch status until complete")
    args = ap.parse_args()

    polys = parse_polys(args.file)
    if not polys:
        sys.exit("no valid polynomials")
    print(f"parsed {len(polys)} polynomials from {args.file}")
    payload = {"payload": {"polynomials": polys}}
    if args.slot:
        payload["slot"] = args.slot
    body = json.dumps(payload).encode()
    UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
          "(KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36")
    headers = {
        "Authorization": f"Bearer {KEY}",
        "Content-Type": "application/json",
        "User-Agent": UA,
        "Accept": "application/json",
    }
    if args.dry:
        print("DRY RUN, would POST", len(body), "bytes")
        return
    data = http_json(API, body, headers, "POST")
    d = data.get("data", {})
    print("ok:", data.get("ok"), "| submissionId:", d.get("submissionId"),
          "| batchId:", d.get("batchId"), "| status:", d.get("submissionStatus"),
          "| queued:", d.get("queuedCount"), "| rejected:", d.get("rejectedCount"))
    batch_id = d.get("batchId")
    if args.poll and batch_id:
        poll_batch(batch_id, headers)


def poll_batch(batch_id, headers, timeout=1800):
    url = ("https://api.sair.foundation/api/public/v1/competitions/igp24/"
           f"submissions/batch/{batch_id}")
    t0 = time.time()
    while time.time() - t0 < timeout:
        d = http_json(url, None, headers, "GET").get("data", {})
        bs = d.get("batchStatus")
        print(f"[{time.time()-t0:6.0f}s] batch {bs} verified={d.get('verifiedCount')} "
              f"failed={d.get('failedCount')} queued={d.get('queuedCount')}")
        if bs != "queued":
            n_ok = n_new = 0
            for p in d.get("polynomials", []):
                st = p.get("status")
                if st == "ok":
                    n_ok += 1
                    flag = "NEW!" if (p.get("scoreable") and not p.get("inBaseline")) else "    "
                    print(f"  {flag} idx={p.get('polynomialIndex')} {p.get('label')} "
                          f"r={p.get('r')} scoreable={p.get('scoreable')} "
                          f"inBaseline={p.get('inBaseline')} "
                          f"src={p.get('discSource')} fieldDisc={p.get('fieldDiscAbs')} "
                          f"baselineDisc={p.get('baselineDiscAbs')}")
                elif st in ("invalid", "error"):
                    print(f"  FAIL idx={p.get('polynomialIndex')} {st}: {p.get('reason')}")
            print(f"summary: ok={n_ok}")
            return
        time.sleep(10)
    print("TIMEOUT waiting for batch")


def http_json(url, body, headers, method):
    req = urllib.request.Request(url, data=body, method=method, headers=headers)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            msg = e.read().decode()[:500]
            print(f"HTTP {e.code}: {msg}", file=sys.stderr)
            if e.code in (429, 503) and attempt < 3:
                time.sleep(5 * (attempt + 1))
                continue
            sys.exit(1)
    sys.exit(1)


if __name__ == "__main__":
    main()
