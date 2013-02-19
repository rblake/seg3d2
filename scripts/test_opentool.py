# run script using exec(open('/Users/aylakhan/devel/seg3d2_is/scripts/test_opentool.py').read()) (change to Windows path if necessary)

import os
import sys

##################################################################################
# local functions (needed before we can import packages)
# TODO: all this should be moved to a module and incorporated into Seg3D's python in the IS branch

def checkLocalPath(path):
  if not os.path.exists(path):
    raise ValueError("Path %s does not exist." % path)

def addLocalPath(path):
  try:
    checkLocalPath(path)
    if path not in sys.path:
      sys.path.append(path)

  except ValueError as err:
    print("ValueError", err)
    sys.exit(err)

def stringEscape(s):
  return s.replace("(","\\(").replace(")","\\)").replace(" ","\\ ")

##################################################################################
# configuration

import configparser

# location of both scripts and settings.ini file
scriptsLocation='/Users/aylakhan/devel/seg3d2_is/scripts'

addLocalPath(scriptsLocation)

# import configuration from settings.ini
# (contains file paths)
config = configparser.ConfigParser()
config.read(os.path.join(scriptsLocation, 'settings.ini'))

# where the Seg3D python scripts are located
seg3DLocation = config.get('Seg3DScriptSettings', 'dir')

# volume and mask location
matlabDataLocation = config.get('DataSettings', 'data_dir')
filename = config.get('DataSettings', 'volume_file')
volumeFilename = os.path.join(matlabDataLocation, filename)
volumeFilename = stringEscape(volumeFilename)

filename = config.get('DataSettings', 'mask_file')
maskFilename = os.path.join(matlabDataLocation, filename)

# location of points file read into Seg3D and labels written by Seg3d
seg3DDataLocation = config.get('Seg3DDataSettings', 'seg3d_data_dir')
filename = config.get('Seg3DDataSettings', 'query_file')
pointsFilename = os.path.join(seg3DDataLocation, filename)

filename = config.get('Seg3DDataSettings', 'labels_file')
labelsFilename = os.path.join(seg3DDataLocation, filename)

addLocalPath(seg3DLocation)

try:
  checkLocalPath(seg3DDataLocation)
except ValueError as err:
  print("ValueError", err)
  sys.exit(err)

##################################################################################
# set up data

import seg3d2utils
import seg3d2

# set up view
result = set(stateid='view::layout', value='single')

# sets up Axial view (0) by default
# returns bool
result = set(stateid='view::active_viewer', value=0)

# returns bool
result = set(stateid='viewer0::view_mode', value='Axial')

# TODO: filepaths don't work with strings (even escaped ones) yet...

dataLayerID = seg3d2utils.importMatlabDataLayer(volumeFilename)




maskLayerID = seg3d2utils.importMatlabSingleMaskLayer(maskFilename)
activatelayer(layerid=maskLayerID)

(vertices, z) = seg3d2utils.getPointsFromFile(pointsFilename)

# tool setup
(toolID, targetStateID, verticesStateID, edgeStateID) = seg3d2utils.openEdgeQueryTool()
result = set(stateid=targetStateID, value=maskLayerID)
result = set(stateid=verticesStateID, value=vertices)
# assuming axial view here...
result = set(stateid='viewer0::slice_number', value=z)
