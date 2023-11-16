#!/usr/bin/env python
import subprocess
import sys

# NB: See docs here: https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#using-cluster-status

jobid = sys.argv[1]

output = str(subprocess.check_output(f"sacct -j {jobid} --format State --noheader | head -1 | awk '{{print $1}}'", shell = True).strip())

# TODO This is probably not comprehensive, what other states?
running_status = ["PENDING", "RUNNING"]

if "COMPLETED" in output:
    print("success")
elif any(r in output for r in running_status):
    print("running")
else:
    print("failed")
