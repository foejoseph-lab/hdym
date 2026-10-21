# Note on `esc_E.txt` / `esc_E.npz` and `esc_orient_0412.txt` / `.npz` (2026-10-20)

These logs were produced on 19–20 Oct 2026 by `islands/code/offdiag_escape.py` and `ym/csrc/hdym_full_escape.c`
as of commit 1723030 (patch 0016, the first version of the kernel), on Joseph's machine (Windows, MinGW gcc,
OpenMP, 24 threads). Two things changed in the kernel afterwards, and this note is kept instead of rerunning.

1. **Thread-dependent floating point (found 20 Oct, fixed in patches 0019–0020).** The step controller called the
   C library's `pow()`, which MinGW implements in x87 code; its last bit depends on the thread's x87 control word,
   and that word differed between the main thread and the OpenMP workers. Identical initial data therefore gave
   two different escape times on different threads (E = 0.156, θ = 0 row: 123 cells at 2965.95, 5 cells — one
   thread's share of 24 — at 1818.84). About 1/24 of the cells in these logs were integrated by that one thread.
   They are valid trajectories of the same equations differing from their neighbours by a last-bit rounding path,
   i.e. independent draws from the same loiter distribution; every statistic quoted from these logs (medians,
   brackets, leading sector, escape/no-escape by tmax) is unaffected. Cell-level, cross-platform bit
   reproducibility holds only for runs made after patch 0020 (the kernel now pins MXCSR and the x87 word in every
   thread and computes the step factor without the C library; the log of every run prints the FP state found and a
   48-copy determinism CHECK).
2. **Step controller.** With `pow()` removed the controller is the standard err^(−1/5) rule (was PI). Same accuracy
   class and the same scipy cross-check tolerances, but a different step sequence, so escape times from these logs
   are not bit-reproducible with the later kernel even on one thread.

Decision (Joseph, 20 Oct): the physics does not change; these logs stand with this note rather than being rerun.
The orient maps at E = 0.156 and 0.30 are made with the fixed kernel (patch 0020 or later).
