""" Helper functions related to treating histogram data that is stored in NumPy
    arrays.
"""

import logging as log
import numpy as np

def proj_2from3D(x, y, i_coord):
  """ Project the 3D histogram {x,y} into a 2D histogram by integrating over the
      coordinate of index i_coord.
  """
  # Get the other coordinate indices and unique values
  other_cs = np.array([i for i in [0,1,2] if not i == i_coord])
  xo = x[:,other_cs]
  unique_xo = np.unique(xo, axis=0)
  
  # Build a mask that sorts the entries into the unique coordinates
  xo_mask = np.all(xo == unique_xo[:,np.newaxis], axis=2)
  
  # Use the mask to perform the projection
  return unique_xo, np.sum(xo_mask * y, axis=1)

def test():
  x = np.array([
    [1.0,0.05,3.3],
    [2.0,0.05,3.4],
    [3.0,0.05,3.3],
    [1.0,0.07,3.3],
  ])
  y = np.array([
    2.2, 4.4, 1.1, 0.1
  ])
  
  log.error(proj_2from3D(x, y, 0))
  log.error(proj_2from3D(x, y, 1))
  log.error(proj_2from3D(x, y, 2))