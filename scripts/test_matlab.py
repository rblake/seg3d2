# run script using exec(open('/Users/ayla/scratch/seg3d2_is/scripts/test_matlab.py').read()) or exec(open('../scripts/test_matlab.py').read()) (change to Windows path if necessary)

import os
import sys
import configparser

dir='/Users/ayla/scratch/seg3d2_is/scripts'
if dir not in sys.path:
  sys.path.append(dir)

config = configparser.ConfigParser()
config.read(os.path.join(dir, 'settings.ini'))
# Matlab config
MATLAB_PATH = config.get('MatlabSettings', 'matlab_path')
shell_working_dir = config.get('MatlabSettings', 'matlab_scripts_dir')

import edgequeryutils
import mlabraw3

print(sys.path)
os.chdir(dir)

eq = edgequeryutils.EdgeQueryUtils(dir)
session = mlabraw3.open(shell_working_dir, MATLAB_PATH)
mlabraw3.eval(session, 'UnConstSpecClust_Seg3D', log=True)

for counter in range(1, eq.maxMatlabIterations):
  eq.processEdgeQuery(counter)
  matlabCommand = "ConstSpecClust_Seg3D(%d)" % counter
  mlabraw3.eval(session, matlabCommand, log=True)


mlabraw3.close(session)
