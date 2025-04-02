import os
import subprocess

def launch_slurm_job(xml_file):
    # Construct the SLURM job name
    job_name = f"phyml_{xml_file[:-4]}"

    # Write the SLURM script to a file
    slurm_script = f"""#!/bin/bash
#SBATCH -c 1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=long
#SBATCH --mem 500
#SBATCH -t 01-12:00
#SBATCH -J {job_name}
#SBATCH --output=phyml"."_%A_%a".".out
#SBATCH --error=phyml"."_%A_%a".".err

srun /shared/home/sguindon/phyml/src/phyml --xml={xml_file}
"""

    # Write the script to a file
    with open("slurm_script.sh", "w") as f:
        f.write(slurm_script)

    # Make the script executable
    subprocess.run(["chmod", "+x", "slurm_script.sh"])

    # Launch the SLURM job
    subprocess.run(["sbatch", "slurm_script.sh"])

# Get the directory from the user
directory = input("Enter the directory path: ")

# Get the path from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python launch_phyml.py <path to directory>")
    sys.exit(1)

directory = sys.argv[1]

# Change into the working directory
os.chdir(directory)

# List all XML files in the directory
xml_files = [f for f in os.listdir(directory) if f.endswith(".xml")]

# Launch phyml for each XML file using SLURM
for xml_file in xml_files:
    launch_slurm_job(os.path.join(directory, xml_file))