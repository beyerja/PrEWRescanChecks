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
  
def sign_str(num, spelled=False):
  """ Return the sign of the number as string, either just the sign or spelled.
  """
  if spelled:
    if num >= 0:
      return "plus"
    else:
      return "minus"
  else:
    if num >= 0:
      return "+"
    else:
      return "-"
  
def process_str(chirality, mu_charge):
  """ Return the latex math string that describes the semilept. WW process with
      the given chirality and muon charge.
  """
  eM_chi, eP_chi = chirality_str(chirality)
  charge_str = sign_str(mu_charge)
  return "{}{}\\rightarrow\\mu^{{{}}}\\nu qq".format(
          eM_chi, eP_chi, charge_str)

