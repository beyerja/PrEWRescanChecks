""" Helper functions related to MatPlotLib.
"""

def get_hist_color(hist):
  """ Get the color of a histogram from the object that was returned by its 
      constructor.
  """
  color_tuple = hist[2][0].get_facecolor()[:-1] # Ignore alpha=0
  return (color_tuple[0], color_tuple[1], color_tuple[2], 1.0)

def get_plot_color(plot):
  """ Get the color of a plot from the object that was returned by its 
      constructor.
  """
  return plot[0].get_color()