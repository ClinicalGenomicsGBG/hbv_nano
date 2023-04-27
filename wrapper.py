import os
import subprocess

subprocess.run("pwd")
subprocess.run("conda env list", shell=True)
subprocess.run("conda activate /Users/daniel/miniconda/envs/minimap2", shell=True)
#command = "conda init zsh; conda activate /Users/daniel/miniconda/envs/minimap2; minimap2 --help"
#subprocess.run(command, shell=True)
#subprocess.run(source /opt/homebrew/Caskroom/miniforge/base)
#subprocess.run("conda init zsh", shell=True)
#subprocess.run("conda init minimap2", shell=True)
#subprocess.run("conda env list", shell=True)
#subprocess.run("conda activate /Users/daniel/miniconda/envs/minimap2", shell=True)
#subprocess.run(". /Users/daniel/miniconda/envs/minimap2", shell=True)
