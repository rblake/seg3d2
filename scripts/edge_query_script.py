# run script using exec(open('../scripts/edge_query_script.py').read())

import os, sys

# must add Seg3D location to python path
#
# sys.path.append('')
#
sys.path.append('/Users/aylakhan/devel/seg3d2_is/bin')
import seg3d2

data_layerid = importseries(filenames='/Users/aylakhan/Downloads/SampleData-1/original_gray.jpg', importer='[ITK FileSeries Importer]', mode='data')
print(data_layerid)

mask_layerid = importlayer(filename='/Users/aylakhan/Downloads/SampleData-1/label_mask.mat', importer='[Matlab Importer]', mode='single_mask')
print(mask_layerid)

# set up Axial view by default
# returns bool
result = set(stateid='view::active_viewer', value=0)

# returns bool
result=set(stateid='viewer0::view_mode', value='Axial')

# returns bool
result=set(stateid='viewer0::flip_vertical', value='true')

verticesFromFile = []

pointsfile = open("points.txt")
for line in pointsfile:
  #print(line)
  floats = [float(points) for points in line.split()]
  #print(floats)
  verticesFromFile.append(floats)

print(verticesFromFile)

#toolid=opentool tooltype=edgequerytool
#set stateid=toolid::vertices value=[[],[],[]]
#
# axial = 0, coronal = 1, sagittal = 2
#edgequery(target=layer_0, slice_type=0, slice_number=86, vertices=[[10,20,30],[10,40,50],[10,50,60]], edges=[L1,L2])
edgequery(target=mask_layerid, slice_type=0, slice_number=86, vertices=verticesFromFile)
