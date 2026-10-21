@echo off
rem Overnight queue, 21 Oct 2026 -- holes 1, 12, 4.  IDEMPOTENT: every run uses --resume and appends to its log, so
rem after a reboot just launch this file again and each run continues from its last completed row (a few minutes
rem lost at most).  Run from the repo root with the Spyder kernel's interpreter, e.g.
rem     C:\GitHub\hdym> rabbit_holes\ym\code\overnight_2026-10-21.bat
rem Optional: register it in Task Scheduler "at log on" for tonight so a restart resumes without you.
setlocal
set PY=C:\Users\Owner\AppData\Local\Programs\Python\Python312\python.exe
cd /d %~dp0..\..\..

echo === %date% %time% queue start >> islands\results\overnight_2026-10-21.log

rem 1. hole 1: phase ensemble at the Krein band, 12 x 256 cells, ~1 h
"%PY%" -B -u islands\code\offdiag_escape.py --section phase --E 0.156 --block 12 --nalpha 12 --nseed 256 --tag esc_phase_0156 --resume >> islands\results\esc_phase_0156.txt 2>&1
echo === %date% %time% phase done (exit %errorlevel%) >> islands\results\overnight_2026-10-21.log

rem 2. hole 12: 64^2 parametric map keeping chk and tcross, ~2.7 h
"%PY%" -B -u islands\code\offdiag_escape.py --section orient --E 0.30 --n 64 --tag esc_orient_030_n64 --resume >> islands\results\esc_orient_030_n64.txt 2>&1
echo === %date% %time% orient done (exit %errorlevel%) >> islands\results\overnight_2026-10-21.log

rem 3. hole 4: section E with crossing times, 20 x 64, ~1 h
"%PY%" -B -u islands\code\offdiag_escape.py --section E --nE 20 --nseed 64 --tag esc_E_cross --resume >> islands\results\esc_E_cross.txt 2>&1
echo === %date% %time% E done (exit %errorlevel%) >> islands\results\overnight_2026-10-21.log

echo === %date% %time% queue end >> islands\results\overnight_2026-10-21.log
endlocal
