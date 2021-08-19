""" Helpers for the naming of parameters etc.
"""

def chirality_str(chirality):
  """ Translate chiralities from eL/pL style to latex style. 
  """
  eM_chi = None
  if "eL" in chirality:
    eM_chi = "e_{{L}}^{{-}}"
  elif "eR" in chirality:
    eM_chi = "e_{{R}}^{{-}}"
    
  eP_chi = None
  if "pL" in chirality:
    eP_chi = "e_{{L}}^{{+}}"
  elif "pR" in chirality:
    eP_chi = "e_{{R}}^{{+}}"
    
  return eM_chi, eP_chi
  