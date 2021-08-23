""" Helper functions related to ROOT RDataFrame.
"""

import Plotting.Naming as PN

def select_mu(rdf, charge):
  """ Return an RDataFrame that only consideres events with muon of the given 
      charge.
  """
  select_str = "(decay_to_mu == 1) && (l_charge == {}1)".format(PN.sign_str(charge))
  return rdf.Filter(select_str)

def skip_0weight(rdf):
  """ Return an RDataFrame that ignores events with 0-weights
  """
  return rdf.Filter("rescan_weights.weight1 > 0.01")