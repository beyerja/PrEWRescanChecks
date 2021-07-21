import matplotlib.pyplot as plt
import numpy as np
import ROOT

# Use the the local helper modules
import sys
sys.path.append("../Helpers") 
import IO.FilenameHelp as IFH
import IO.SysHelpers as ISH
import IO.TGCConfigReader as ITCR
import Plotting.DefaultFormat as PDF
import Plotting.MPLHelp as PMH
import Plotting.ROOTHistHelp as PRHH

def add_rel_dev_plot(ax, tgc_name, x, y_dict, n_bins, hist_range, 
                     SM_xsection, lumi):
  """ Add the plots showing the relative deviations from the SM.
  """
  y_SM = y_dict["SM"]
  y_SM_0replaced = y_SM + np.where(y_SM > 0, y_SM, 1)
  
  d_dict = {}
  for point, y in y_dict.items():
    d_dict[point] = np.where(y_SM > 0, (y - y_SM)/ y_SM_0replaced, 0)
  
  d_P = d_dict["{} +".format(tgc_name)]
  d_M = d_dict["{} -".format(tgc_name)]
  
  dmin_P = np.amin(d_P, axis=1)
  dmax_P = np.amax(d_P, axis=1)
  dmean_P = np.mean(d_P, axis=1)
  
  dmin_M = np.amin(d_M, axis=1)
  dmax_M = np.amax(d_M, axis=1)
  dmean_M = np.mean(d_M, axis=1)
  
  scale=100. # in percent
  
  h_P_min = ax.hist(x, weights=dmin_P * scale, alpha=0.9, lw=1, range=hist_range, bins=n_bins, histtype=u'step')
  color_P = PMH.get_hist_color(h_P_min)
  h_P_max = ax.hist(x, weights=dmax_P * scale, color=color_P, alpha=0.9, lw=1, range=hist_range, bins=n_bins, histtype=u'step')
  ax.bar(x=h_P_max[1][:-1], height=h_P_max[0]-h_P_min[0], bottom=h_P_min[0], width=np.diff(h_P_max[1]), align='edge', linewidth=0, color=color_P, alpha=0.25, zorder=-1)
  h_P_mean = ax.hist(x, weights=dmean_P * scale, color=color_P, ls="dashed", alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  
  h_M_min = ax.hist(x, weights=dmin_M * scale, alpha=0.9, lw=1, range=hist_range, bins=n_bins, histtype=u'step')
  color_M = PMH.get_hist_color(h_M_min)
  h_M_max = ax.hist(x, weights=dmax_M * scale, color=color_M, alpha=0.9, lw=1, range=hist_range, bins=n_bins, histtype=u'step')
  ax.bar(x=h_M_max[1][:-1], height=h_M_max[0]-h_M_min[0], bottom=h_M_min[0], width=np.diff(h_M_max[1]), align='edge', linewidth=0, color=color_M, alpha=0.25, zorder=-1)
  h_M_mean = ax.hist(x, weights=dmean_M * scale, color=color_M, ls="dotted", alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  
  p_SM = ax.hist([], ls="solid", color="black", alpha=0.5, lw=1, range=hist_range, bins=n_bins, histtype=u'step')

def add_dev_sig_plot(ax, tgc_name, x, y_dict, n_bins, hist_range, 
                     SM_xsection, lumi):
  """ Add the plots showing the significance of each deviation.
  """
  y_SM = y_dict["SM"]
  
  n_MC = np.sum(y_SM)
  normed = lambda y: y * SM_xsection * lumi / n_MC
  
  # significance here: 
  # how many sigma if we see signal? => (signal - background) / sqrt(signal)
  s_dict = {}
  for point, y in y_dict.items():
    y_0replaced = y + np.where(y > 0, y, 1)
    s_dict[point] = np.where(y_SM>0, np.abs(normed(y) - normed(y_SM)) / np.sqrt(normed(y_0replaced)), 0)

  tgc_name_P = "{} +".format(tgc_name)
  tgc_name_M = "{} -".format(tgc_name)
  s_P = s_dict[tgc_name_P]
  s_M = s_dict[tgc_name_M]
  
  smin_P = np.amin(s_P, axis=1)
  smax_P = np.amax(s_P, axis=1)
  smean_P = np.mean(s_P, axis=1)
  
  smin_M = np.amin(s_M, axis=1)
  smax_M = np.amax(s_M, axis=1)
  smean_M = np.mean(s_M, axis=1)
  
  scale=100. # in percent
  
  title_dict = { "g1z +" : r"$\Delta g_1^{Z} = +\delta$",
                 "g1z -" : r"$\Delta g_1^{Z} = -\delta$",
                 "ka +" : r"$\Delta \kappa_{\gamma} = +\delta$",
                 "ka -" : r"$\Delta \kappa_{\gamma} = -\delta$",
                 "la +" : r"$\Delta \lambda_{\gamma} = +\delta$",
                 "la -" : r"$\Delta \lambda_{\gamma} = -\delta$" }
  
  h_P_min = ax.hist(x, weights=smin_P, alpha=0.9, lw=1, range=hist_range, bins=n_bins, histtype=u'step')
  color_P = PMH.get_hist_color(h_P_min)
  h_P_max = ax.hist(x, weights=smax_P, color=color_P, alpha=0.9, lw=1, range=hist_range, bins=n_bins, histtype=u'step')
  ax.bar(x=h_P_max[1][:-1], height=h_P_max[0]-h_P_min[0], bottom=h_P_min[0], width=np.diff(h_P_max[1]), align='edge', linewidth=0, color=color_P, alpha=0.25, label=title_dict[tgc_name_P], zorder=-1)
  h_P_mean = ax.hist(x, weights=smean_P, color=color_P, ls="dashed", alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  
  h_M_min = ax.hist(x, weights=smin_M, alpha=0.9, lw=1, range=hist_range, bins=n_bins, histtype=u'step')
  color_M = PMH.get_hist_color(h_M_min)
  h_M_max = ax.hist(x, weights=smax_M, color=color_M, alpha=0.9, lw=1, range=hist_range, bins=n_bins, histtype=u'step')
  ax.bar(x=h_M_max[1][:-1], height=h_M_max[0]-h_M_min[0], bottom=h_M_min[0], width=np.diff(h_M_max[1]), align='edge', linewidth=0, color=color_M, alpha=0.25, label=title_dict[tgc_name_M], zorder=-1)
  h_M_mean = ax.hist(x, weights=smean_M, color=color_M, ls="dotted", alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  
  line_dummy = ax.plot([], [], color='black', ls="-.", label="bin-to-bin mean")
  
def create_reweighting_plot(x, y_dict, output_base, 
                            obs_name, obs_range, lumi, SM_xsection, dev_scale, 
                            output_formats=["pdf","png"]):
  """ Create the 3D reweighting plot for one observable.
  """
  for tgc_name in ["g1z","ka","la"]:
    fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True, figsize=(10,10))
    ax_up, ax_down = axs
    
    n_bins = obs_range[0]
    hist_range=(obs_range[1],obs_range[2])
    
    # Upper plot: Relative deviation wrt. SM
    add_rel_dev_plot(ax_up, tgc_name, x, y_dict, n_bins, hist_range, SM_xsection, lumi)
    
    ax_up.set_title("Bin-to-bin sensitivity range in 3D")
    ax_up.set_ylabel(r"$\frac{\# weighted - \# SM}{\# SM} [\%]$", fontsize=30)
    
    # Lower plot: Significance of deviation (in case it would be measured)
    add_dev_sig_plot(ax_down, tgc_name, x, y_dict, n_bins, hist_range, SM_xsection, lumi)  
    
    ax_down.set_ylabel(r"$\frac{\# weighted - \# SM}{\sqrt{\# weighted}}$", fontsize=30)
    ax_down.legend(ncol=1, fontsize=16, title=r"$\delta={}$, $L={}$ab$^{{-1}}$".format(dev_scale, lumi/1000), title_fontsize=16)
    ax_down.set_xlim(hist_range)
    ax_down.set_ylim(0, ax_down.get_ylim()[1])
    ax_down.set_xlabel(obs_name)

    for format in output_formats:
      format_dir = "{}/{}".format(output_base,format)
      ISH.create_dir(format_dir)
      fig.savefig("{}/{}_{}_3DSensitivity.{}".format(format_dir, obs_name, tgc_name, format))

    plt.close(fig)


def check_reweighting(root_file, tgc_config_path, tgc_point_path, output_dir,
                      observables, lumi, tree_name = "WWObservables"):
  """ Check what reweighting does to the observable distributions for the given 
      events in the ROOT file.
  """
  rdf = ROOT.RDataFrame(tree_name, root_file).Filter("rescan_weights.weight1 > 0.01")
  tcr = ITCR.TGCConfigReader(tgc_config_path, tgc_point_path)
  
  chirality = IFH.find_chirality(root_file)
  output_base = "{}/3DSensitivity/{}/".format(output_dir,chirality)
  
  index_dict = {
    "g1z +": tcr.point_index([1,0,0]),
    "g1z -": tcr.point_index([-1,0,0]),
    "ka +": tcr.point_index([0,1,0]),
    "ka -": tcr.point_index([0,-1,0]),
    "la +": tcr.point_index([0,0,1]),
    "la -": tcr.point_index([0,0,-1])
  }
  dev_scale = tcr.scale
  
  # Dict that stores histograms for each observables and each dev point
  hist_dict = {}
  
  obs_names = ["costh_Wminus_star", "costh_l_star", "phi_l_star"]
  obs_1 = observables[obs_names[0]]
  obs_2 = observables[obs_names[1]]
  obs_3 = observables[obs_names[2]]
  h_SM_def = ("SM", "SM", obs_1[0], obs_1[1], obs_1[2],
                          obs_2[0], obs_2[1], obs_2[2],
                          obs_3[0], obs_3[1], obs_3[2])
              
  # Book all the histograms as lazy actions of the RDataFrame
  hist_dict["SM"] = rdf.Histo3D(h_SM_def, obs_names[0], obs_names[1], 
                                          obs_names[2])
  # Histograms with TGC deviation
  for point, index in index_dict.items():
    i_str = str(index+1) # weight naming starts with index 1
    h_def = ("w"+i_str, "w"+i_str, obs_1[0], obs_1[1], obs_1[2],
                                   obs_2[0], obs_2[1], obs_2[2],
                                   obs_3[0], obs_3[1], obs_3[2])
    weight_name = "rescan_weights.weight"+i_str
    hist_dict[point] = rdf.Histo3D(h_def, 
                                   obs_names[0], obs_names[1], obs_names[2], 
                                   weight_name)
    
  # Get the SM cross section, needed for correct uncertainty weighting
  SM_xsection = rdf.Mean("cross_section")
      
  # Get histogram data as numpy arrays 
  xy_dict = {}
  for point, hist in hist_dict.items():
    xy_dict[point] = PRHH.TH3_to_arrays(hist) 
  x = xy_dict["SM"][0]
  
  # # Make the plots for each observable
  for o in range(len(obs_names)):
    obs_name = obs_names[o]
    obs_range = observables[obs_name]
    
    # Find unique x coordinate values
    x_obs = x[:,o]
    unique_x_obs = np.unique(x[:,o].round(decimals=4))
    
    # Determine the y values found for each unique x value
    y_obs_dict = {}
    for point, xy in xy_dict.items():
      y = xy[1]
      y_obs_dict[point] = np.array([y[np.where(np.abs(x_obs - _x)<1.e-4)[0]] for _x in unique_x_obs])
      
    create_reweighting_plot(unique_x_obs, y_obs_dict, output_base, 
                            obs_name, obs_range, 
                            lumi, SM_xsection.GetValue(), dev_scale)
      
def main():
  """ Create check plots that investigate how the weights in the given ROOT file
      affect the distribution.
  """
  ROOT.EnableImplicitMT() # Enable multithreading in RDataFrame
  ROOT.gROOT.SetBatch(True) # Don't show graphics at runtime
  PDF.set_default_mpl_format()
  
  config_dir = "/afs/desy.de/group/flc/pool/beyerjac/TGCAnalysis/SampleProduction/MCProduction/PrEWSampleProduction/scripts/config"
  tgc_config_path = config_dir + "/tgc.config"
  tgc_point_path = config_dir + "/tgc_dev_points_g1z_ka_la.config"
  
  output_dir = "../../output"
  
  observables = { "costh_Wminus_star" : (30, -1., 1.),
                  "costh_l_star" :      (30, -1., 1.),
                  "phi_l_star":         (30, -np.pi, np.pi) }
  
  lumi = 1000 # 2ab^-1 total ~> 1ab^-1 for each of the two allowed chiral states
  
  RL_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eR_pL.root"
  LR_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eL_pR.root"
  
  check_reweighting(RL_path, tgc_config_path, tgc_point_path, output_dir, observables, lumi)
  check_reweighting(LR_path, tgc_config_path, tgc_point_path, output_dir, observables, lumi)
  

if __name__ == "__main__":
  main()