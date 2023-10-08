#!/usr/bin/env python
import sys

pad_length = 3

i = sys.argv[1]
if len(i) > 2:
    print(i)
else:
    print('0'*(pad_length-len(i)) + i)
