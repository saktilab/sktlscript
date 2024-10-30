#!/usr/bin/env python3

import gentools as gt
import sys

print(gt.loadgen(sys.argv[1]).volume/(0.529177**3.0))
