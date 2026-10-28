@echo off
rem Day queue, 22 Oct 2026 -- hole 4 / hole 8: the low-E "order-one ghost, no escape" regime, run four times longer.
rem Last night's esc_E_cross: at E = 0.0585 and 0.110 the ghost reaches |w| = 1 in ~60/64 cells and ~50 of them do not
rem escape by tmax = 3000 (loiters of 16-21 periods, censored).  Same section, tmax 12000, 20 energies over the gap
rem [0.05, 0.13] where the diagonal ghost is Floquet-stable and the off-diagonal bands lead.  Idempotent (--resume, >>);
rem relaunch after any restart.  (islands\results\day_2026-10-22.log is the .bat's own log; not committed.)  ~3-4 h worst case on 24 threads; most low-E rows are much faster.
setlocal
set PY=C:\Users\Owner\AppData\Local\Programs\Python\Python312\python.exe
cd /d %~dp0..\..\..
echo === %date% %time% day queue start >> islands\results\day_2026-10-22.log
"%PY%" -B -u islands\code\offdiag_escape.py --section E --Emin 0.05 --Emax 0.13 --nE 20 --nseed 64 --tmax 12000 --nchk 240 --tag esc_E_lowE_long --resume >> islands\results\esc_E_lowE_long.txt 2>&1
echo === %date% %time% day queue done (exit %errorlevel%) >> islands\results\day_2026-10-22.log
endlocal
