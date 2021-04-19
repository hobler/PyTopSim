"""
Cleans the working directory.

Copyright (C) 2021 Gerhard Hobler  (license: GPLv3 or higher)
"""
import os
import shutil
import subprocess

########################################################################
program = 'pytopsim'        # name of the program
########################################################################
# This parameter determines where the program is built:
########################################################################
root_dir = '/home/hobler/work/' + program   # root of program workspace
########################################################################
#command = 'git describe --tags --exact-match || ' + \
#          'git symbolic-ref -q --short HEAD'
#command = 'if test "$(shell git status -z)"; then ' + \
#          'git symbolic-ref --short HEAD; else ' + \
#          'git describe --tags --exact-match || git symbolic-ref --short HEAD; fi'
#p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
p = subprocess.Popen('git symbolic-ref --short HEAD', 
                     shell=True, stdout=subprocess.PIPE)
branch = p.stdout.read()[:-1]   # strip newline


def clean(project):
    """
    Remove directory tree of project in the workspace.
    
    Argument:
        project (string): Name of the project (e.g., 'code', 'scripts', etc.)
    """
    proj_dir = os.path.join(root_dir, branch, project)
    
    if os.path.isdir(proj_dir):
        shutil.rmtree(proj_dir)
        print 'Directory ' + proj_dir + ' has been removed.'
    else:
        print 'Project is already clean.'