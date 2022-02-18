""" Helper functions related to ROOT RDataFrame.
"""

import Plotting.Naming as PN

import logging as log

def select_mu(rdf, charge):
  """ Return an RDataFrame that only consideres events with muon of the given 
      charge.
  """
  if charge is None:
    return rdf
  else:
    select_str = "(decay_to_mu == 1) && (l_charge == {}1)".format(PN.sign_str(charge))
    return rdf.Filter(select_str)

def skip_0weight(rdf):
  """ Return an RDataFrame that ignores events with 0-weights
  """
  return rdf.Filter("rescan_weights.weight1 > 0.01")
  
def select_mumu(rdf, m_min=None, m_max=None, Z_direction=None):
  """ Return an RDataFrame that only considers di-muon events in the given mass
      frame.
  """
  cut_str = "(f_pdg == 13)"
  
  if not (m_min is None) and not (m_max is None):
    cut_str += " && (m_ff > {}) && (m_ff < {})".format(m_min, m_max)
  else:
    log.debug("Not applying m_ff cut, need non-None cut values.")
  
  if Z_direction == "FZ":
    cut_str += " && (pz_ff > 0)"
  elif Z_direction == "BZ":
    cut_str += " && (pz_ff < 0)"
  else:
    log.debug("Not applying Z direction cut, got {}".format(Z_direction))
  
  return rdf.Filter(cut_str)