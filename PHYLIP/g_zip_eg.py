import os
import tarfile

from contextlib import contextmanager

fname = '/scratch/d85/mc7636/Yeast/D2S_archive/Genomes_for_AFphylogeny_red_40_26_D2S_cp.tz.gz'

@contextmanager
def change_dir(destination):
    try:
        cwd = os.getcwd()
        os.chdir(destination)
        yield
    finally:
        os.chdir(cwd)

dirname = os.path.dirname(fname)
filename = os.path.basename(fname)

with change_dir(dirname):
    tar = tarfile.open(filename, "r:gz")
    tar.extractall()
    tar.close()