import os
import seg3d2

# TODO: filepaths don't work with strings (even escaped ones) yet...

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
  # assuming axial view here...
  zlocation = 0
  for pointList in pointsFromFile:
    if len(pointList) == 2:
      pointList.append(0)
    else:
      zlocation = pointList[2]

  print("Edge query vertices from Matlab: ", pointsFromFile)
  return (pointsFromFile, zlocation)

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
  print(filepath)
  dataLayerID = seg3d2.importlayer(filename=filepath, importer='[Matlab Importer]', mode='data')
  return dataLayerID[0]

def importMatlabSingleMaskLayer(filepath):
  maskLayerID = seg3d2.importlayer(filename=filepath, importer='[Matlab Importer]', mode='single_mask')
  return maskLayerID[0]

def openEdgeQueryTool():
  toolID = seg3d2.opentool(tooltype='edgequerytool')
  targetStateID = toolID + "::target"
  verticesStateID = toolID + "::vertices"
  edgeStateID = toolID + "::edge"
  return (toolID, targetStateID, verticesStateID, edgeStateID)
