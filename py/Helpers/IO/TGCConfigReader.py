import numpy as np

class TGCConfigReader:
  """ Reader that can interpret the TGC configuration used in the 
      PrEWSampleProduction.
  """
  
  def __init__(self, config_path, point_path):
    self.scale = None
    self.read_config(config_path)
    
    self.dev_points = None
    self.read_points(point_path)
    
  def read_config(self, config_path):
    """ Run the config file (trivial bash notation works in python) and store
        the variables.
    """
    var_dict = {}
    exec(open(config_path).read(), globals(), var_dict)
    self.scale = var_dict["scale"]
    
  def read_points(self, point_path):
    """ Read and interpret the file containing the TGC deviation points.
    """
    dev_points = []
    for line in open(point_path):
      dev_points.append(np.array(line.strip().split(" "), dtype=float))
    self.dev_points = np.array(dev_points)
      
  def point_index(self, point):
    """ Given a TGC deviation point (g1Z ka la) return it's index in the point 
        array.
    """
    matches = np.where(np.all(self.dev_points == np.array(point), axis=1))
    if len(matches) == 0:
      raise Exception("No point {} found.".format(point))
    elif len(matches) > 1:
      raise Exception("More than one point {} found.".format(point))
    return matches[0][0]