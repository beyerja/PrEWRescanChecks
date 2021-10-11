""" Helpers for the naming of parameters etc.
"""

def chirality_str(chirality):
  """ Translate chiralities from eL/pL style to latex style. 
  """
  eM_chi = None
  if chirality is None:
    eM_chi = "e^{{-}}"
  elif "eL" in chirality:
    eM_chi = "e_{{L}}^{{-}}"
  elif "eR" in chirality:
    eM_chi = "e_{{R}}^{{-}}"
    
  eP_chi = None
  if chirality is None:
    eP_chi = "e^{{+}}"
  elif "pL" in chirality:
    eP_chi = "e_{{L}}^{{+}}"
  elif "pR" in chirality:
    eP_chi = "e_{{R}}^{{+}}"
    
  return eM_chi, eP_chi
  
def sign_str(num, spelled=False):
  """ Return the sign of the number as string, either just the sign or spelled.
  """
  if num is None:
    return ""
  elif spelled:
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

def fermion_str(fermion):
  """ Return the latex string to a given fermion.
  """
  if fermion == "mu": 
    return "\\mu"
  else:
    raise Exception("Unknown fermion {}".format(fermion))
    
def difermion_mass_str(label):
  """ Latex Labelling of a difermion mass range.
  """
  if label == "return_to_Z":
    return "m_{{Z}}"
  elif label == "high_Q2":
    return "250\\,GeV"
  else:
    raise Exception("Unknown mass label {}".format(label))

def difermion_process_str(fermion, chirality, label, Z_direction=None):
  """ Return the latex math string that describes the difermion produciton with
      the given chirality in the given mass range, optional with Z flight 
      direction.
  """
  f_str = fermion_str(fermion)
  eM_chi, eP_chi = chirality_str(chirality)
  m_str = difermion_mass_str(label)
  
  if Z_direction:
    m_str = "({},{})".format(m_str, Z_direction)
  else:
    m_str = "({})".format(m_str)
    
  return "{}{}\\rightarrow{}{} {}".format(
          eM_chi, eP_chi, f_str, f_str, m_str)

def observable_str(obs_name, process):
  """ Return the LaTeX string for a given observable.
  """
  if process == "WW":
    if obs_name == "costh_Wminus_star":
      return "\\cos\\theta_{W^-}"
    elif obs_name == "costh_l_star":
      return "\\cos\\theta_l^*"
    elif obs_name == "phi_l_star":
      return "\\phi_l^*"
    else:
      raise Exception("Unknown WW observable {}".format(obs_name))
  else:
    raise Exception("Unknown process {}".format(process))
