#!/usr/bin/env python3
"""
test_checker.py -- fast regression tests for the trusted checker (seconds, not the full proof).

  python3 test_checker.py [path/to/cert16.json]

1. run_chain accepts a 4-step slice of chain B of the committed certificate (base + ghost).
2. Same slice with the ghost frames G/Ginv deleted from one step: must raise CheckFailed
   (the per-step frame switch found by the review of 2026-10-11).
3. Initial data in the outward-rounding slack of steps[0] but outside the exact set c + r:
   must raise CheckFailed (review of 2026-10-11).
4. The outward-rounded printing never moves an endpoint inward.
"""
import copy, json, os, random, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cert_check as drv
import cert_lemmas as cl

here = os.path.dirname(os.path.abspath(__file__))
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "results", "cert16.json")
cert = json.load(open(path))
cl.set_grid(int(cert.get("K", 200)))
p = int(cert["p"])
raw = copy.deepcopy(cert["chainB"])
raw["steps"] = raw["steps"][:4]
ch = drv.parse_chain(raw)
init = [(a[0] + b[0], a[0] + b[1]) for a, b in zip(ch["steps"][0]["c"], ch["steps"][0]["r"])]
init = [((lo + hi) / 2, (lo + hi) / 2) for lo, hi in init]          # a point inside steps[0]

def expect_fail(label, fn):
    try:
        fn()
    except cl.CheckFailed as e:
        print(f"ok   {label}: rejected ({e})")
        return
    raise SystemExit(f"FAIL {label}: accepted")

cl.run_chain(ch, p, [0, 1], "B", False, init)
print("ok   1: 4-step slice of chain B accepted (base + ghost)")

bad = copy.deepcopy(raw)
del bad["steps"][2]["G"], bad["steps"][2]["Ginv"]
expect_fail("2: ghost frames missing on one step", lambda: cl.run_chain(drv.parse_chain(bad), p, [0, 1], "B", False, init))

off = copy.deepcopy(ch)
c0, r0 = off["steps"][0]["c"], off["steps"][0]["r"]
tiny = Fr(1, 2 ** (cl.K + 50))                                     # far below the rounding grid
r0[0] = (r0[0][0], r0[0][1] - tiny)                                 # exact set shrinks off the grid
point = c0[0][0] + r0[0][1] + tiny / 2                              # outside the exact set ...
assert cl.set_box(c0, off["steps"][0]["A"], r0)[0][1] >= point      # ... but inside the rounded box
init3 = list(init); init3[0] = (point, point)
expect_fail("3: initial point in the rounding slack", lambda: cl.run_chain(off, p, [0], "B", False, init3))

rng = random.Random(1)
for _ in range(20000):
    x = Fr(rng.randint(-10**30, 10**30), rng.randint(1, 10**30)) * Fr(10) ** rng.randint(-20, 20)
    sig = rng.randint(1, 17)
    lo, hi = Fr(drv._dec(x, False, sig)), Fr(drv._dec(x, True, sig))
    assert lo <= x <= hi, (x, sig, lo, hi)
print("ok   4: outward printing contains the exact value (20000 random cases)")
print("all checker regression tests passed")
