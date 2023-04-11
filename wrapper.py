import os
import subprocess

"""whatever = subprocess.run(
    'echo $SHELL', 
    shell=True, executable='/bin/zsh',
    check=True)
"""

subprocess.run("pwd")
#command = "conda init zsh; conda activate /Users/daniel/miniconda/envs/minimap2; minimap2 --help"
#subprocess.run(command, shell=True)
subprocess.run("conda init zsh", shell=True)
#subprocess.run("conda env list", shell=True)
subprocess.run("conda activate minimap2", shell=True)
#subprocess.run(". /Users/daniel/miniconda/envs/minimap2", shell=True)
