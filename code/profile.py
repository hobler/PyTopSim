"""
Usage: python profile.py configfile

Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
"""

from PyTopSim import main
import cProfile


cProfile.run('main()')