 
#!/usr/bin/env python3

import os
import platform
import subprocess
import sys


NTHREADS = "10"
JOB_LIST = "problem_file_list.txt"
WALL_TIME_KEY = ">> total job wall time"


def get_warp3d_exe():
    home = os.environ.get("WARP3D_HOME")
    if not home:
        sys.exit("ERROR: WARP3D_HOME is not defined")

    system = platform.system()
    machine = platform.machine().lower()

    if system == "Darwin":
        exe_name = "warp3d_silicon.exe" if machine in ("arm64", "aarch64") \
                                     else "warp3d_Intel.exe"
        exe = os.path.join(home, "run_macOS", exe_name)

    elif system == "Linux":
        exe = os.path.join(home, "run_linux", "warp3d.exe")

    elif system == "Windows":
        exe = os.path.join(home, "run_windows", "warp3d.exe")

    else:
        sys.exit(f"ERROR: unsupported platform: {system}")

    if not os.path.isfile(exe):
        sys.exit(f"ERROR: executable not found: {exe}")

    print(f"Using WARP3D executable: {exe}")
    return exe


def find_wall_time_line(outname):
    if not os.path.exists(outname):
        return None

    with open(outname, "r", errors="ignore") as f:
        for line in f:
            if WALL_TIME_KEY in line:
                return line.strip()

    return None


def extract_wall_time_seconds(line):
    if not line:
        return None

    fields = line.split()
    try:
        return float(fields[-1])
    except (ValueError, IndexError):
        return None


def run_cleanup():
    try:
        subprocess.run("./clean_dir.bash", shell=True, check=False)
    except Exception as e:
        print(f"cleanup failed: {e}")


def run_case(exe, fname):
    outname = f"{fname}.out"

    if os.path.exists(outname):
        os.remove(outname)

    print(f"\n>>> running: {fname}")

    rc = None
    try:
        with open(fname, "r") as fin, open(outname, "w") as fout:
            result = subprocess.run(
                [exe, NTHREADS],
                stdin=fin,
                stdout=fout,
                stderr=subprocess.STDOUT,
                check=False
            )
        rc = result.returncode

    except Exception as e:
        run_cleanup()
        if os.path.exists(outname):
            os.remove(outname)
        return True, f"exception: {e}", None

    wall_time_line = find_wall_time_line(outname)
    wall_time_seconds = extract_wall_time_seconds(wall_time_line)

    failed = False
    reason = ""

    if rc != 0:
        failed = True
        reason = f"return code {rc}"

    if wall_time_line:
        print(f"{fname}: {wall_time_line}")
    else:
        if not failed:
            failed = True
            reason = "wall time not found"
        print(f"{fname}: wall time not found")

    run_cleanup()

    if os.path.exists(outname):
        os.remove(outname)

    return failed, reason, wall_time_seconds


def main():
    exe = get_warp3d_exe()

    with open(JOB_LIST, "r") as f:
        files = [line.strip() for line in f if line.strip()]

    failures = []
    wall_times = []

    for fname in files:
        failed, reason, wall_time_seconds = run_case(exe, fname)

        if failed:
            failures.append((fname, reason))

        if wall_time_seconds is not None:
            wall_times.append((fname, wall_time_seconds))

    print("\n====================")
    print("Failure Summary")
    print("====================")

    if failures:
        for fname, reason in failures:
            print(f"{fname}: {reason}")
        print(f"\nTotal failures: {len(failures)} / {len(files)}")
    else:
        print("All jobs completed successfully.")

    print("\n====================")
    print("Wall Time Summary")
    print("====================")

    if wall_times:
        total_time = sum(t for _, t in wall_times)
        min_job, min_time = min(wall_times, key=lambda x: x[1])
        max_job, max_time = max(wall_times, key=lambda x: x[1])
        avg_time = total_time / len(wall_times)

        print(f"Jobs with wall times found: {len(wall_times)} / {len(files)}")
        print(f"Total wall time: {total_time:.2f} sec")
        print(f"Average wall time: {avg_time:.2f} sec")
        print(f"Minimum wall time: {min_time:.2f} sec  ({min_job})")
        print(f"Maximum wall time: {max_time:.2f} sec  ({max_job})")
    else:
        print("No wall times found.")


if __name__ == "__main__":
    main()    
    