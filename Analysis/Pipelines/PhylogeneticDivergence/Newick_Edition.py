import os

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Phylogenetic_tree/CorePhylo/SH_TEST')

filelist = os.listdir(
    '/Users/u1866168/Documents/OneDrive - University of Warwick/Experiments/Phylogenetic_tree/CorePhylo/SH_TEST')

for file in sorted(filelist):
    if file.endswith('.treefile'):
        with open(file, "r+") as f:
            old = f.read()  # read everything in the file
            f.seek(0)  # rewind
            f.write("[ " + str(file.rsplit('.', 5)[0]) + " ]" + old)  # write the new line before
