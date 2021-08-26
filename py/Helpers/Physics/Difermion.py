""" Code related to the difermion parametrisation.
"""

def dif_param_LR(cos_th_name, xs0, Ae, Af, ef, k0, dk):
  """ Difermion parametrisation in the difermion rest frame for the left-handed
      electron right-handed positron initial state.
  """
  return f"3./8. * {xs0} * (1.0 + {Ae})/2.0 * ( (1. + ({k0} + {dk})/2.0) + ({ef} + 2.0 * {Af}) * {cos_th_name} + (1.0 - 3.0 * ({k0} + {dk})/2.0) * {cos_th_name}*{cos_th_name} )"

def dif_param_RL(cos_th_name, xs0, Ae, Af, ef, k0, dk):
  """ Difermion parametrisation in the difermion rest frame for the right-handed
      electron left-handed positron initial state.
  """
  return f"3./8. * {xs0} * (1.0 - {Ae})/2.0 * ( (1. + ({k0} - {dk})/2.0) + ({ef} - 2.0 * {Af}) * {cos_th_name} + (1.0 - 3.0 * ({k0} - {dk})/2.0) * {cos_th_name}*{cos_th_name} )"

      
def get_dif_param(chirality):
  """ Get the difermion parametrisation for a given chirality.
  """
  if (chirality == "eLpR"):
    return dif_param_LR
  elif (chirality == "eRpL"):
    return dif_param_RL
  else:
    raise Exception("Unknown chirality " + chirality)
    
class DifParamSet:
  """ Class summarizing the set of difermion parameters.
  """
  def __init__(self, xs0, Ae, Af, ef, k0, dk):
    self.xs0 = xs0
    self.Ae = Ae
    self.Af = Af
    self.ef = ef
    self.k0 = k0
    self.dk = dk
  
  def arr(self):
    return [self.xs0, self.Ae, self.Af, self.ef, self.k0, self.dk]
    
  def add_dev(self, d_xs0, d_Ae, d_Af, d_ef, d_k0, d_dk):
    """ Return the same parameter set with some deviations added.
    """
    return DifParamSet( self.xs0 + d_xs0, self.Ae + d_Ae, self.Af + d_Af, 
                        self.ef + d_ef, self.k0 + d_k0, self.dk + d_dk )
    
def SM_parset(fermion, mass_label):
  """ Return the Standard Model set of difermion parameters for a given mass 
      range.
      The total chiral cross section is set to 1 in all of them (does not matter
      for reweighting tests, can be set from dataset).
  """
  difpar_set = None
  if fermion == "mu":
    if mass_label == "return_to_Z":
      difpar_set = DifParamSet(1.0, 0.21360014, 0.20281099, 0.01580906, 0.07471141, 0.00059199)
    elif mass_label == "high_Q2":
      difpar_set = DifParamSet(1.0, 0.11251847, 0.03217479, 1.42594481, 0.00033356, 0.00031470)
  
  if difpar_set:
    return difpar_set
  else:
    raise Exception("No difpar set for {} {}".format(fermion, mass_label))

def weight(chirality, cos_th_name, param_set_ini, param_set_dev):
  """ Calculate the weight for an event of given chirality and at given 
      cos(theta) when the difermion parameters are changed from the initial set 
      to a set with deviations.
  """
  dif_param = get_dif_param(chirality)
  p_ini = param_set_ini.arr()
  p_dev = param_set_dev.arr()
  return "({})/({})".format( dif_param(cos_th_name, *p_dev),
                             dif_param(cos_th_name, *p_ini) )
  
def add_weight_branch(rdf, w_name, cos_th_name, param_set_ini, param_set_dev, 
                      chirality):
  """ Add a difermion weight branch with the given name to the RDataFrame.
      The weight represents a change from the initial parameter set to the 
      deviated one.
      Needs the name of the cos(theta) branch to make the weight dependend on 
      it.
  """
  w_fct = "return {};".format(
            weight(chirality, cos_th_name, param_set_ini, param_set_dev))
  return rdf.Define(w_name, w_fct)
  