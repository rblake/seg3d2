# run script using exec(open('/Users/ayla/scratch/seg3d2_is/seg3d2_is/scripts/test_matlab.py').read()) or exec(open('../scripts/test_matlab.py').read()) (change to Windows path if necessary)

import os
import sys

dir='/Users/ayla/scratch/seg3d2_is/scripts'
if dir not in sys.path:
  sys.path.append(dir)

import edgequeryutils
import mlabraw3

print(sys.path)
os.chdir(dir)

eq = edgequeryutils.EdgeQueryUtils(dir)
session = mlabraw3.open(dir)
mlabraw3.eval(session, 'UnConstSpecClust_Seg3D', log=True)

#for counter in range(1, eq.max_iterations):
for counter in range(1, 3):
  eq.processEdgeQuery(counter)

mlabraw3.close(session)
