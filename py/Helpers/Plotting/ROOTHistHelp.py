""" Helper function related to the ROOT histograms.
"""

import numpy as np

def TH1_to_arrays(th1):
  """ Transform a ROOT TH1 into numpy x and y arrays.
  """
  # Bins excluding underflow / overflow
  bins = np.arange(th1.GetNbinsX()) + 1
  x = []
  y = []
  for bin in bins:
    x.append(th1.GetBinCenter(int(bin)))
    y.append(th1.GetBinContent(int(bin)))
  return np.array(x), np.array(y)