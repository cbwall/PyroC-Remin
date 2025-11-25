# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

"""

import sys

file_in = sys.argv[1]
name_in = '.'.join(file_in.split('.')[0:-1])

with open(file_in, 'r') as autobod_in, open(name_in + '.clean.csv', 'w') as autobod_out:
    for line in autobod_in.readlines():
        line = line.rstrip()
        line = line.split()
        if len(line) != 15:
            continue
        else:
            line = ','.join(line)
            print(line, file = autobod_out)