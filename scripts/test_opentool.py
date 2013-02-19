# run script using exec(open('../scripts/test_opentool.py').read()) (change to Windows path if necessary)

import os
import sys
import configparser

def getPointsFromFile(dir, filename):
  pointsFile = open(os.path.join(dir, filename))
  pointsFromFile = []
  for line in pointsFile:
    floats = [float(points) for points in line.split()]
    pointsFromFile.append(floats)

  pointsFile.close()

  print(len(pointsFromFile))

  # 2D case - add default z location (0)
  for pointList in pointsFromFile:
    if len(pointList) == 2:
      pointList.append(0)

  print("Edge query vertices from Matlab: ", pointsFromFile)
  return pointsFromFile

def writeLabelsToFile(dir, filename, selectedEdge):
  # TODO: can probably find a better way to do this...
  if selectedEdge == 0:
    labels = '1 0'
  elif selectedEdge == 1:
    labels = '0 1'
  else:
    # exception?
    print("Invalid edge")
    return;

  edgeFile = open(os.path.join(dir, filename), 'wt')
  edgeFile.write(labels)
  edgeFile.close()

config = configparser.ConfigParser()
config.read('../scripts/settings.ini')
#config.read('/Users/aylakhan/devel/seg3d2_is/scripts/settings.ini')

seg3DLocation = config.get('Seg3DSettings', 'dir')
matlabDataLocation = config.get('DataSettings', 'data_dir')
seg3DDataLocation = config.get('DataSettings', 'seg3d_data_dir')
pointsFilename = config.get('DataSettings', 'query_file')
labelsFilename = config.get('DataSettings', 'labels_file')

if not os.path.exists(seg3DLocation):
  print("Path %s to Seg3D does not exist." % seg3DLocation)
  sys.exit(0)

if not os.path.exists(seg3DDataLocation):
  print("Path %s to data does not exist." % seg3DDataLocation)
  sys.exit(0)

import seg3d2


# set up Axial view by default
# returns bool
result=set(stateid='view::active_viewer', value=0)

# returns bool
result=set(stateid='viewer0::view_mode', value='Axial')

vertices = getPointsFromFile(seg3DDataLocation, pointsFilename)

# TODO: get actual layerid
activatelayer(layerid='layer_1')
toolid=opentool(tooltype='edgequerytool')
print(toolid)
vertices_stateid=toolid+"::vertices"
set(stateid=vertices_stateid, value=vertices)
