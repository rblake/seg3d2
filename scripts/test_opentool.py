# run script using exec(open('/Users/aylakhan/devel/seg3d2_is/scripts/test_opentool.py').read()) or exec(open('../scripts/test_opentool.py').read()) (change to Windows path if necessary)

import edgequeryutils
eq = edgequeryutils.EdgeQueryUtils('/Users/aylakhan/devel/seg3d2_is/scripts')

print(eq.pointsFilename)
