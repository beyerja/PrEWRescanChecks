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

def add_abs_dev_plot(ax, x, xy_dict, n_bins, hist_range):
  """ Add the plots showing the absolute distributions.
  """
  h_SM = ax.hist(x, weights=xy_dict["SM"][1], ls="-", label="SM", range=hist_range, bins=n_bins, histtype=u'step')
  SM_color = PMH.get_hist_color(h_SM)
  
  h_g1z_P = ax.hist(x, weights=xy_dict["g1z +"][1], ls="dashed", label=r"$\Delta g_1^{Z} = +\delta$", range=hist_range, bins=n_bins, histtype=u'step')
  g1z_color = PMH.get_hist_color(h_g1z_P)
  h_g1z_M = ax.hist(x, weights=xy_dict["g1z -"][1], ls="dotted" , label=r"$\Delta g_1^{Z} = -\delta$", color=g1z_color, range=hist_range, bins=n_bins, histtype=u'step')
  
  h_ka_P = ax.hist(x, weights=xy_dict["ka +"][1], ls="dashed", label=r"$\Delta \kappa_{\gamma} = +\delta$", range=hist_range, bins=n_bins, histtype=u'step')
  ka_color = PMH.get_hist_color(h_ka_P)
  h_ka_M = ax.hist(x, weights=xy_dict["ka -"][1], ls="dotted" , label=r"$\Delta \kappa_{\gamma} = -\delta$", color=ka_color, range=hist_range, bins=n_bins, histtype=u'step')
  
  h_la_P = ax.hist(x, weights=xy_dict["la +"][1], ls="dashed", label=r"$\Delta \lambda_{\gamma} = +\delta$", range=hist_range, bins=n_bins, histtype=u'step')
  la_color = PMH.get_hist_color(h_la_P)
  h_la_M = ax.hist(x, weights=xy_dict["la -"][1], ls="dotted" , label=r"$\Delta \lambda_{\gamma} = -\delta$", color=la_color, range=hist_range, bins=n_bins, histtype=u'step')
  
def add_rel_dev_plot(ax, x, xy_dict, SM_xsection, lumi):
  """ Add the plots showing the relative deviations from the SM.
  """
  y_SM = xy_dict["SM"][1]
  d_g1z_P = (xy_dict["g1z +"][1] - y_SM) / y_SM
  d_g1z_M = (xy_dict["g1z -"][1] - y_SM) / y_SM
  d_ka_P  = (xy_dict["ka +"][1] - y_SM) / y_SM
  d_ka_M  = (xy_dict["ka -"][1] - y_SM) / y_SM
  d_la_P  = (xy_dict["la +"][1] - y_SM) / y_SM
  d_la_M  = (xy_dict["la -"][1] - y_SM) / y_SM
  
  scale=100. # in percent
  p_SM = ax.plot(x, np.zeros(len(x)), ls="dashed")
  SM_line = ax.axline((np.amin(x), 0), (np.amax(x), 0), ls='--')
  SM_color = PMH.get_plot_color(p_SM)
  rel_SM_unc = 1/np.sqrt(y_SM / np.sum(y_SM) * SM_xsection * lumi) 
  ax.fill_between(x, -rel_SM_unc*scale, +rel_SM_unc*scale, color=SM_color, alpha=0.25, label=r"SM unc. for {}ab$^{{-1}}$".format(lumi/1000.))
  
  p_g1z_P = ax.plot(x, d_g1z_P * scale, ls="dashed")
  g1z_color = PMH.get_plot_color(p_g1z_P)
  p_g1z_M = ax.plot(x, d_g1z_M * scale, ls="dotted", color=g1z_color)
  
  p_ka_P = ax.plot(x, d_ka_P * scale, ls="dashed")
  ka_color = PMH.get_plot_color(p_ka_P)
  p_ka_M = ax.plot(x, d_ka_M * scale, ls="dotted", color=ka_color)
  
  p_la_P = ax.plot(x, d_la_P * scale, ls="dashed")
  la_color = PMH.get_plot_color(p_la_P)
  p_la_M = ax.plot(x, d_la_M * scale, ls="dotted", color=la_color)

def create_reweighting_plot(xy_dict, output_base, obs_name, obs_range, 
                            lumi, SM_xsection, dev_scale, 
                            output_formats=["pdf","png"]):
  """ Create the reweighting plot for one observable.
  """
  fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True, figsize=(10,10))
  ax_up, ax_down = axs
  
  x = xy_dict["SM"][0]
  n_bins = obs_range[0]
  hist_range=(obs_range[1],obs_range[2])
  
  add_abs_dev_plot(ax_up, x, xy_dict, n_bins, hist_range)
  
  ax_up.set_ylabel("#MC Events")
  ax_up.set_yscale('log')
  ax_up.legend(ncol=3, fontsize=13, loc='lower center', title=r"$\delta={}$".format(dev_scale), title_fontsize=13)
  ax_up.set_xlim(hist_range)
  
  # Calculate the relative ratios
  add_rel_dev_plot(ax_down, x, xy_dict, SM_xsection, lumi)
  
  ax_down.legend(fontsize=17)
  
  ax_down.set_ylabel(r"$\frac{\# weighted - \# SM}{\# SM} [\%]$")
  ax_down.set_xlabel(obs_name)

  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    fig.savefig("{}/{}_ReweightCheck.{}".format(format_dir, obs_name, format))

  plt.close(fig)


def check_reweighting(root_file, tgc_config_path, tgc_point_path, output_dir,
                      observables, lumi, tree_name = "WWObservables"):
  """ Check what reweighting does to the observable distributions for the given 
      events in the ROOT file.
  """
  rdf = ROOT.RDataFrame(tree_name, root_file).Filter("rescan_weights.weight1 > 0.01")
  tcr = ITCR.TGCConfigReader(tgc_config_path, tgc_point_path)
  
  chirality = IFH.find_chirality(root_file)
  output_base = "{}/ReweightCheck/{}/".format(output_dir,chirality)
  
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
  
  # Book all the histograms as lazy actions of the RDataFrame
  for obs_name, obs_range in observables.items():
    hist_dict[obs_name] = {}
    
    # Standard model histograms for comparison
    h_SM_def = (obs_name+"SM", obs_name+"SM", 
                obs_range[0], obs_range[1], obs_range[2])
    hist_dict[obs_name]["SM"] = rdf.Histo1D(h_SM_def, obs_name)
    
    # Histograms with TGC deviation
    for point, index in index_dict.items():
      i_str = str(index+1) # weight naming starts with index 1
      h_def = (obs_name+"_w"+i_str, obs_name+"_w"+i_str, 
               obs_range[0], obs_range[1], obs_range[2])
      weight_name = "rescan_weights.weight"+i_str
      hist_dict[obs_name][point] = rdf.Histo1D(h_def, obs_name, weight_name)
    
  # Get the SM cross section, needed for correct uncertainty weighting
  SM_xsection = rdf.Mean("cross_section")
      
  # Make the plots for each observable
  for obs_name, obs_range in observables.items():
    # First transform hists to arrays for easier matplotlib plotting
    xy_dict = {}
    for point, hist in hist_dict[obs_name].items():
      xy_dict[point] = PRHH.TH1_to_arrays(hist) 
    create_reweighting_plot(xy_dict, output_base, obs_name, obs_range, lumi, 
                            SM_xsection.GetValue(), dev_scale)
      
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
                  "phi_l_star":         (30, 0., np.pi) }
  
  lumi = 2000 # 2ab^-1 to calculate the expected SM uncertainty
  
  RL_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eR_pL.root"
  LR_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eL_pR.root"
  
  check_reweighting(RL_path, tgc_config_path, tgc_point_path, output_dir, observables, lumi)
  check_reweighting(LR_path, tgc_config_path, tgc_point_path, output_dir, observables, lumi)
  

if __name__ == "__main__":
  main()