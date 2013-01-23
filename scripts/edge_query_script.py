# run script using exec(open('../scripts/edge_query_script.py').read())

import os, sys

# must add Seg3D location to python path
seg3d_location='/Users/aylakhan/devel/seg3d2_is/bin'
if not os.path.exists(seg3d_location):
  #print("Path ", seg3d_location, " to Seg3D does not exist.")
  print("Path %s to Seg3D does not exist." % seg3d_location)
  sys.exit(0)

sys.path.append(seg3d_location)

import mlabraw3
import seg3d2

data_dir='/Users/aylakhan/Downloads/SampleData-1'

data_layerid = importseries(filenames=os.path.join(data_dir, 'original_gray.jpg'), importer='[ITK FileSeries Importer]', mode='data')
print(data_layerid)

session = mlabraw3.open()
print('Matlab session opened')

# test
counter = 0
while counter < 1:
# stops looping when points file is empty
#while True:

  # matlab eval call
  # TODO: replace with actual Matlab function
  # TODO: error handling
  print("Call Matlab")
  mlabraw3.eval(session, 'A=zeros(4, 4)', log=True)

  # assuming layer mask will always have the same name...
  # TODO: should check file for changes - checksum?
  mask_layerid = importlayer(filename=os.path.join(data_dir, 'label_mask.mat'), importer='[Matlab Importer]', mode='single_mask')
  print(mask_layerid)

  # set up Axial view by default
  # returns bool
  result = set(stateid='view::active_viewer', value=0)

  # returns bool
  result=set(stateid='viewer0::view_mode', value='Axial')

  # returns bool
  result=set(stateid='viewer0::flip_vertical', value='true')

  verticesFromFile = []
  edges = []

  pointsfile = open(os.path.join(data_dir, 'points.txt'))
  for line in pointsfile:
    floats = [float(points) for points in line.split()]
    verticesFromFile.append(floats)

  print("Edge query vertices from Matlab: ", verticesFromFile)

  if not verticesFromFile:
    print("Done")
    break

  #toolid=opentool tooltype=edgequerytool
  #set stateid=toolid::vertices value=[[],[],[]]
  #
  # axial = 0, coronal = 1, sagittal = 2
  #edgequery(target=layer_0, slice_type=0, slice_number=86, vertices=[[10,20,30],[10,40,50],[10,50,60]], edges=[L1,L2])
  #edgequery(target=mask_layerid, slice_type=0, slice_number=86, vertices=verticesFromFile, edges=[])
  #edgequery(target=mask_layerid, slice_type=0, slice_number=86, vertices=[[10,20,1],[10,40,1],[10,50,1]], edges=[])

  # TODO: should layer be deleted?
  #deletelayers(layers=mask_layerid)

  counter=counter+1

mlabraw3.close(session)
print('Matlab session closed')
