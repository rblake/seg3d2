import os
import sys

dir='/Users/ayla/scratch/seg3d2_is/bin'
if dir not in sys.path:
  sys.path.append(dir)

#dir='/Users/ayla/scratch/seg3d2_is/bin/Seg3D2.app/Contents/MacOS/lib'
#if dir not in sys.path:
  #sys.path.append(dir)

dir='/Users/ayla/scratch/seg3d2_is/scripts'
if dir not in sys.path:
  sys.path.append(dir)

print(sys.path)
os.chdir(dir)

import mlabraw3

session = mlabraw3.open(dir)
mlabraw3.eval(session, 'A=zeros(4, 4)', log=True)
mlabraw3.eval(session, 'A=test2', log=True)
mlabraw3.close(session)
print('session closed')
