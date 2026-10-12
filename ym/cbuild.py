"""
cbuild.py -- run the C compiler for the scripts' on-the-fly kernels, robustly on Windows.

check_call(cmd) replaces subprocess.check_call(cmd) for compiler invocations:
  * the compiler's own bin directory is put FIRST on the PATH of the compiler process.  Other
    software on the PATH can ship older copies of the DLLs gcc's helpers load (zlib1.dll,
    libzstd.dll, libwinpthread-1.dll, libgcc_s_seh-1.dll ...); found first, they make cc1.exe die
    before it prints anything (seen 2026-10-12 from the Spyder 6 console, whose runtime folder
    spyder-runtime\\Library\\bin precedes C:\\msys64\\mingw64\\bin).
  * on failure the compiler's output and exit code are printed (Windows status codes in hex:
    0xC0000135 DLL not found, 0xC0000139 entry point not found), instead of failing silently.
  * on Windows the compiler's bin directory is also added to the DLL search path of this
    process, so the compiled kernel finds libgomp-1.dll and its dependencies when loaded.
Nothing outside the compiler process is changed.
"""
import os
import shutil
import subprocess
import sys

_added = set()


def _bin_dir(exe):
    path = shutil.which(exe)
    return os.path.dirname(os.path.abspath(path)) if path else None


def check_call(cmd):
    env = os.environ.copy()
    bindir = _bin_dir(cmd[0])
    if bindir:
        env['PATH'] = bindir + os.pathsep + env.get('PATH', '')
        if sys.platform == 'win32' and bindir not in _added:
            try:
                os.add_dll_directory(bindir)
                _added.add(bindir)
            except (AttributeError, OSError):
                pass
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        code = r.returncode
        hexcode = f" (0x{code & 0xFFFFFFFF:08X})" if code < 0 or code > 255 else ""
        print(f"compiler failed, exit code {code}{hexcode}: {' '.join(cmd)}", flush=True)
        for stream in (r.stdout, r.stderr):
            if stream.strip():
                print(stream.rstrip(), flush=True)
        if not (r.stdout.strip() or r.stderr.strip()):
            print("(no output: a helper program of the compiler could not start; on Windows "
                  "usually a DLL conflict on the PATH)", flush=True)
        raise subprocess.CalledProcessError(code, cmd, r.stdout, r.stderr)
    if r.stderr.strip():
        print(r.stderr.rstrip(), flush=True)
