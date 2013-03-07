#!/usr/bin/env python

""" A quick and extremely dirty hack to wrap matlabpipe/matlabcom as if they
were mlabraw.

Author: Dani Valevski <daniva@gmail.com>
License: MIT
"""
import os
import sys

is_win = sys.platform == 'win32'
if is_win:
  from matlabcom import MatlabCom as MatlabConnection
  from matlabcom import MatlabError as error
else:
  from matlabpipe3 import MatlabPipe as MatlabConnection
  from matlabpipe3 import MatlabError as error

#try:
#import settings
#except:
#class settings:
#MATLAB_PATH = '/usr/sci/OSX/matlab_r2011a/MATLAB_R2011a.app'


def open(matlab_scripts_dir, matlab='guess'):
  if is_win:
    ret = MatlabConnection()
    ret.open()
  else:
    if matlab != 'guess':
      matlab_path = matlab + '/bin/matlab'
    else:
      matlab_path = 'guess'

    try:
      ret = MatlabConnection(matlab_path, timeout=1000)
      ret.open(print_matlab_welcome=False, working_directory=matlab_scripts_dir)
    except:
      print('Could not open matlab, is it in %s?' % matlab_path)
  return ret
  
def close(matlab):
  matlab.close()

def eval(matlab, exp, log=False):
  if log or is_win:
    matlab.eval(exp)
  else:
    matlab.eval(exp, print_expression=False, on_new_output=None)
  return ''

#def get(matlab, var_name):
#  return matlab.get(var_name)
#
#def put(matlab, var_name, val):
#  matlab.put({var_name : val})
