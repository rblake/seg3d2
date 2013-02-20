# run script using exec(open('/Users/aylakhan/devel/seg3d2_is/scripts/test_opentool.py').read()) (change to Windows path if necessary)
#
# sys.path.append('/Users/aylakhan/devel/seg3d2_is/scripts')
# import edgequeryutils
# eq = edgequeryutils.EdgeQueryUtils('/Users/aylakhan/devel/seg3d2_is/scripts')
# eq.processEdgeQuery()

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
# set up configuration

# location of both scripts and settings.ini file
scriptsLocation='/Users/aylakhan/devel/seg3d2_is/scripts'

#addLocalPath(scriptsLocation)

import seg3d2utils

(volumeFilename, maskFilename, pointsFilename, labelsFilename) = seg3d2utils.importConfiguration(scriptsLocation) 

##################################################################################
# set up data
# TODO: filepaths don't work with strings (even escaped ones) yet...

import seg3d2

seg3d2utils.setupGlobalView()
dataLayerID = seg3d2utils.importMatlabDataLayer(volumeFilename)

seg3d2utils.processEdgeQuery(maskFilename, pointsFilename, labelsFilename)
