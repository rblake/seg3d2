import os
import sys
import configparser
import seg3d2
import threading

# TODO: filepaths don't work with strings (even escaped ones) yet...

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

def processFilepath(configureParser, configSection, option, basename):
  fileName = configureParser.get(configSection, option)
  fullPath = os.path.join(basename, fileName)
  return stringEscape(fullPath)

def importConfiguration(scriptsLocation):
  # import configuration from settings.ini
  # (contains file paths)
  config = configparser.ConfigParser()
  config.read(os.path.join(scriptsLocation, 'settings.ini'))

  # where the Seg3D python scripts are located
  seg3DLocation = config.get('Seg3DScriptSettings', 'dir')

  # add custom scripts to PYTHONPATH
  addLocalPath(seg3DLocation)

  # volume and mask location
  matlabDataLocation = config.get('DataSettings', 'data_dir')

  volumeFilename = processFilepath(config, 'DataSettings', 'volume_file', matlabDataLocation)
  maskFilename = processFilepath(config, 'DataSettings', 'mask_file', matlabDataLocation)

  # location of points file read into Seg3D and labels written by Seg3d
  seg3DDataLocation = config.get('Seg3DDataSettings', 'seg3d_data_dir')

  pointsFilename = processFilepath(config, 'Seg3DDataSettings', 'query_file', seg3DDataLocation)
  labelsFilename = processFilepath(config, 'Seg3DDataSettings', 'labels_file', seg3DDataLocation)

  try:
    checkLocalPath(seg3DDataLocation)
  except ValueError as err:
    print("ValueError", err)
    sys.exit(err)

  return (volumeFilename, maskFilename, pointsFilename, labelsFilename)

def getPointsFromFile(filepath):
  pointsFile = open(filepath)
  pointsFromFile = []
  for line in pointsFile:
    floats = [float(points) for points in line.split()]
    pointsFromFile.append(floats)

  pointsFile.close()

  if len(pointsFromFile) < 3:
    # TODO: error class
    raise Exception('ListError', 'Points list should contain 3 points.') 

  # 2D case - add default z location (0)
  for pointList in pointsFromFile:
    if len(pointList) == 2:
      pointList.append(0)

  print("Edge query vertices from Matlab: ", pointsFromFile)
  return pointsFromFile

def writeLabelsToFile(filepath, selectedEdge):
  # TODO: can probably find a better way to do this...
  if selectedEdge == 0:
    labels = '1 0'
  elif selectedEdge == 1:
    labels = '0 1'
  else:
    raise ValueError("Invalid edge ", selectedEdge)

  edgeFile = open(filepath, 'wt')
  edgeFile.write(labels)
  edgeFile.close()


def importMatlabDataLayer(filepath):
  dataLayerID = seg3d2.importlayer(filename=filepath, importer='[Matlab Importer]', mode='data')
  return dataLayerID[0]

def importMatlabSingleMaskLayer(filepath):
  maskLayerID = seg3d2.importlayer(filename=filepath, importer='[Matlab Importer]', mode='single_mask')
  return maskLayerID[0]

# TODO: replace return tuple with named tuple...
def openEdgeQueryTool():
  toolID = seg3d2.opentool(tooltype='edgequerytool')
  targetStateID = toolID + "::target"
  verticesStateID = toolID + "::vertices"
  edgeStateID = toolID + "::edge"
  saveStateID = toolID + "::save"
  return (toolID, targetStateID, verticesStateID, edgeStateID, saveStateID)

def setupGlobalView():
  # set(..) returns bool
  result = seg3d2.set(stateid='view::layout', value='single')
  # sets up Axial view (0) by default
  result = seg3d2.set(stateid='view::active_viewer', value=0)
  result = seg3d2.set(stateid='viewer0::view_mode', value='Axial')

def setupTool(vertices, maskLayerID):
  (toolID, targetStateID, verticesStateID, edgeStateID, saveStateID) = openEdgeQueryTool()
  result = seg3d2.set(stateid=targetStateID, value=maskLayerID)
  result = seg3d2.set(stateid=verticesStateID, value=vertices)
  # assuming axial view here...
  z = vertices[0][2]
  result = seg3d2.set(stateid='viewer0::slice_number', value=z)
  return (toolID, targetStateID, verticesStateID, edgeStateID, saveStateID)

def processEdgeQuery(maskFilename, pointsFilename, labelsFilename):
  try:
    maskLayerID = importMatlabSingleMaskLayer(maskFilename)
    seg3d2.activatelayer(layerid=maskLayerID)

    vertices = getPointsFromFile(pointsFilename)
    (toolID, targetStateID, verticesStateID, edgeStateID, saveStateID) = setupTool(vertices, maskLayerID)

    c = threading.Condition()
    c.acquire()

    # wait for about 2 minutes
    MAX_COUNT = 60
    counter = 1
    with c:
      while not seg3d2.get(stateid=saveStateID):
        if counter > 60:
          break

        counter += 1
        c.wait(2.0)

    saveEdges = seg3d2.get(stateid=saveStateID)

    if saveEdges:
      selectedEdge = seg3d2.get(stateid=edgeStateID)
      writeLabelsToFile(labelsFilename, selectedEdge)
      print("Saved edge labels to ", labelsFilename)
    else:
      print("Timeout before edge could be selected")
  except Exception as err:
    print(err)

  finally:
    seg3d2.closetool(toolid=toolID)
    seg3d2.deletelayers(layers=maskLayerID) 


## class CheckSave:
##   import seg3d2
##   def __init__(self, stateid):
##     self.stateid = stateid

##   def __call__(self):
##     print('called!!!', seg3d2.get(stateid=self.stateid))
##     return seg3d2.get(stateid=self.stateid)
