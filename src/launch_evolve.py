import os
import sys
import shutil

# Get the path from the command-line argument
if len(sys.argv) != 3:
    print("Usage: python launch_evolve.py <path> <model>")
    sys.exit(1)

path = sys.argv[1]
model = sys.argv[2]

# Create the GTR directory if it doesn't exist
if not os.path.exists(path + model + "/"):
    os.makedirs(path + model + "/")

# Change into the GTR directory
os.chdir(path + model + "/")

shutil.copy(path + "evolve", path + model + "/")

# Launch the evolve command-line tool 100 times
for i in range(100):
    os.system(path + model + "/" + "evolve" + " --n_otu=50 --seq_len=1000 --no_colalias -c 1 --r_seed=" + str(i) + " -m " + model)