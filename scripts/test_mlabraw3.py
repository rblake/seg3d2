# run script using exec(open('/Users/ayla/scratch/seg3d2_is/scripts/test_mlabraw3.py').read()) or exec(open('../scripts/test_mlabraw3.py').read()) (change to Windows path if necessary)

import os
import sys

MATLAB_PATH = '/usr/sci/OSX/matlab_r2011a/MATLAB_R2011a.app'

dir='/Users/ayla/scratch/seg3d2_is/bin'
if dir not in sys.path:
  sys.path.append(dir)

dir='/Users/ayla/scratch/seg3d2_is/scripts'
if dir not in sys.path:
  sys.path.append(dir)

print(sys.path)
#os.chdir(dir)

import mlabraw3

session = mlabraw3.open(dir, MATLAB_PATH)
mlabraw3.eval(session, 'A=zeros(4, 4)', log=True)
mlabraw3.eval(session, 'A=test2', log=True)
mlabraw3.close(session)
print('session closed')
