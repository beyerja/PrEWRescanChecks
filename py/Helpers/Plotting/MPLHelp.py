""" Helper functions related to MatPlotLib.
"""

import matplotlib.pyplot as plt

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
  
def reverse_color_cycle():
  """ Reverse the current color cycle.
  """
  plt.rcParams['axes.prop_cycle'] = plt.cycler('color', reversed(plt.rcParams['axes.prop_cycle'].by_key()['color']))
  
def skip_n_colors(ax, n):
  """ Skip n colors in the matplotlib color cycle.
  """
  for _ in range(n):
    ax.plot([],[])
    ax.bar([],[])