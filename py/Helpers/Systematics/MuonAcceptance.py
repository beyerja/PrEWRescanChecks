""" Code related to muon acceptance cuts.
"""


def get_costh_cut(cut_val, center_shift, width_shift, costh_branch):
  """ Get the string that describes the cos(theta) cut for the given cut values.
      Four inputs are needed:
        The central cut value (same for both sides +-).
        The change that of the acceptance center (center_shift).
        The change that of the acceptance width (width_shift).
        The name of the cos(Theta) branch.
  """
  pos_cut =   abs(cut_val) + center_shift + width_shift/2.0
  neg_cut = - abs(cut_val) + center_shift - width_shift/2.0
  cut_str = "({} > {}) && ({} < {})".format(costh_branch, neg_cut, costh_branch,
                                            pos_cut)
  return cut_str
  
def get_ndim_costh_cut(cut_val, center_shift, width_shift, costh_branches):
  """ Applying the same cos(theta) cut to multiple branches (/particles) at the
      same time. 
  """
  cut_str = ""
  for costh_branch in costh_branches:
    if cut_str != "":
      cut_str += " && "  
    cut_str += get_costh_cut(cut_val, center_shift, width_shift, costh_branch)
  return cut_str