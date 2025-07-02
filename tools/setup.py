#!/usr/bin/env python3

"""
Script and Functions for setting up required Python packages for RACS comparison tools
"""

import subprocess
import sys
import pkg_resources

# Define dependencies needed to be installed
# Analogous to CRAN packages
pckges = ["openpyxl", "plotly", "pandas", "numpy", "matplotlib"]

##############################################################################################
# Functions for checking whether dependencies are installed and install needed packages

def needed_packages(pckges, other_pckgs=None):
    if other_pckgs is None:
        other_pckgs = []

    # Combine package lists
    all_pckges = pckges + other_pckgs

    # Check what packages are already installed
    installed = {pkg.key for pkg in pkg_resources.working_set}
    need_to_install = [pkg for pkg in all_pckges if pkg.lower() not in installed]

    if need_to_install:
        print("Requested packages:")
        print(all_pckges)
        print("Installing:", need_to_install)
        for pkg in need_to_install:
            subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])


def check_version(pckges):
    import platform
    import sys
    print("\nSystem Info:")
    print(f"Python version: {platform.python_version()}")
    print(f"Executable: {sys.executable}\n")

    print("Installed package versions:")
    for pkg in pckges:
        try:
            version = pkg_resources.get_distribution(pkg).version
            print(f"{pkg}: {version}")
        except pkg_resources.DistributionNotFound:
            print(f"{pkg}: NOT INSTALLED")


##############################################################################################

# Check and install required packages
needed_packages(pckges)

# Display versions of the installed/required packages
check_version(pckges)
