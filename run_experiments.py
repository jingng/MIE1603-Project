import subprocess
import argparse

# try:
#     subprocess.call(["python3", "mip_model.py"])
# except:
#     print("Couldn't find a feasible solution in time limit for mip_model")

try:
    subprocess.call(["python3", "mip_lb.py"])
except:
    print("Couldn't find a feasible solution in time limit for mip_lb")

