import subprocess
import glob
import shutil
import os

#file_names = glob.glob("/Users/xschmd/Desktop/test/*")
#folder_names = glob.glob("/Users/xschmd/Desktop/test2/*")



#base = os.path.basename(os.path.normpath(folder_names[1]))
#new = f"{base}.fastq.gz"

#print(len(folder_names))    #how many folders there are
"""
rootdir = "/Users/xschmd/Desktop/test2/"

#gives the name of the files in each folder
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        print(os.path.basename(os.path.normpath(subdir)))


for root, dirs, files in os.walk(rootdir):
    for f in dirs:
        print(f)
        # list all files in directory
"""
new = glob.glob("/Users/xschmd/Desktop/test2/barcode01/*")
print(os.path.basename(os.path.normpath(new[1])))

#this is functional
with open("outfile.fastq.gz",'wb') as wfd:
    for f in new:
        with open(f,'rb') as fd:
            shutil.copyfileobj(fd, wfd)
            print(f)

#OLD CODE (Maybe useful later...)
"""
with open('output_file.fastq.gz','wb') as wfd:
    for f in file_names:
        with open(f,'rb') as fd:
            shutil.copyfileobj(fd, wfd)
"""
#Concatenate all the .fastq.gz files in each barcode folder

"""for i in range(1, 18):
    if i < 10:
        i =f"0{i}"   
    bc_path = f"{path}{i}"    #the path to the barcode folder
    command = f"cat {bc_path}/*.fastq.gz > {bc_path}/bc{i}.fastq.gz"    #the command to concatenate the files
    print(command)
    #subprocess.run(command, shell=True)
"""
