""" Helper functions related to file name conventions.
"""

import re

def find_chirality(file_path):
  """ Find the chirality in the given file path.
  """
  eM_chi = re.search(r"e(L|R)",file_path).group(0)
  eP_chi = re.search(r"p(L|R)",file_path).group(0)
  return "{}{}".format(eM_chi,eP_chi)