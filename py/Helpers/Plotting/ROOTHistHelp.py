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
  
def TH3_to_arrays(th3):
  """ Transform a ROOT TH3 into numpy x and y arrays.
  """
  # Find bin centers and values, skip overflow/underflow bins
  x=[]
  y=[]
  for x_bin in range(0, th3.GetNbinsX()+2):
    for y_bin in range(0, th3.GetNbinsY()+2):
      for z_bin in range(0, th3.GetNbinsZ()+2):
        bin = th3.GetBin(x_bin, y_bin, z_bin)

        # Skip overflow and underflow bins
        if th3.IsBinUnderflow(bin) or th3.IsBinOverflow(bin): 
          continue

        # Collect bin information and fill into arrays
        _x = [th3.GetXaxis().GetBinCenter(x_bin),
              th3.GetYaxis().GetBinCenter(y_bin),
              th3.GetZaxis().GetBinCenter(z_bin)]
        if not x:
          x = [_x] # Make sure to get correct array structure
        else:
          x.append(_x)  

        y.append(th3.GetBinContent(bin))

  return np.array(x), np.array(y)