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
import Plotting.Naming as PN
import Plotting.RDFHelp as PRH
import Plotting.ROOTHistHelp as PRHH

def add_rel_dev_plot(ax, x, xy_dict, n_bins, hist_range, lumi):
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
  p_SM = ax.hist([], ls="solid", color="black", alpha=0.9, lw=1, range=hist_range, bins=n_bins, histtype=u'step')
  
  p_g1z_P = ax.hist(x, weights=d_g1z_P * scale, ls="dashed", alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  g1z_color = PMH.get_hist_color(p_g1z_P)
  p_g1z_M = ax.hist(x, weights=d_g1z_M * scale, ls="dotted", color=g1z_color, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  
  p_ka_P = ax.hist(x, weights=d_ka_P * scale, ls="dashed", alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ka_color = PMH.get_hist_color(p_ka_P)
  p_ka_M = ax.hist(x, weights=d_ka_M * scale, ls="dotted", color=ka_color, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  
  p_la_P = ax.hist(x, weights=d_la_P * scale, ls="dashed", alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  la_color = PMH.get_hist_color(p_la_P)
  p_la_M = ax.hist(x, weights=d_la_M * scale, ls="dotted", color=la_color, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')

def add_dev_sig_plot(ax, x, xy_dict, n_bins, hist_range, MC_norm):
  """ Add the plots showing the significance of each deviation.
  """
  y_SM = xy_dict["SM"][1]
  n_MC = np.sum(y_SM)
  normed = lambda y: y * MC_norm
  
  # significance here: 
  # how many sigma if we see signal? => (signal - background) / sqrt(signal)
  s_g1z_P = np.abs(normed(xy_dict["g1z +"][1]) - normed(y_SM)) / np.sqrt(normed(xy_dict["g1z +"][1]))
  s_g1z_M = np.abs(normed(xy_dict["g1z -"][1]) - normed(y_SM)) / np.sqrt(normed(xy_dict["g1z -"][1]))
  s_ka_P  = np.abs(normed(xy_dict["ka +"][1]) - normed(y_SM)) / np.sqrt(normed(xy_dict["ka +"][1]))
  s_ka_M  = np.abs(normed(xy_dict["ka -"][1]) - normed(y_SM)) / np.sqrt(normed(xy_dict["ka -"][1]))
  s_la_P  = np.abs(normed(xy_dict["la +"][1]) - normed(y_SM)) / np.sqrt(normed(xy_dict["la +"][1]))
  s_la_M  = np.abs(normed(xy_dict["la -"][1]) - normed(y_SM)) / np.sqrt(normed(xy_dict["la -"][1]))
  
  p_g1z_P = ax.hist(x, weights=s_g1z_P, label=r"$\Delta g_1^{Z} = +\delta$", ls="dashed", alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  g1z_color = PMH.get_hist_color(p_g1z_P)
  p_g1z_M = ax.hist(x, weights=s_g1z_M, label=r"$\Delta g_1^{Z} = -\delta$", ls="dotted", color=g1z_color, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  
  p_ka_P = ax.hist(x, weights=s_ka_P, label=r"$\Delta \kappa_{\gamma} = +\delta$", ls="dashed", alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ka_color = PMH.get_hist_color(p_ka_P)
  p_ka_M = ax.hist(x, weights=s_ka_M, label=r"$\Delta \kappa_{\gamma} = -\delta$", ls="dotted", color=ka_color, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  
  p_la_P = ax.hist(x, weights=s_la_P, label=r"$\Delta \lambda_{\gamma} = +\delta$", ls="dashed", alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  la_color = PMH.get_hist_color(p_la_P)
  p_la_M = ax.hist(x, weights=s_la_M, label=r"$\Delta \lambda_{\gamma} = -\delta$", ls="dotted", color=la_color, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  
  
def create_reweighting_plot(xy_dict, output_base, obs_name, obs_range, 
                            lumi, MC_norm, dev_scale, 
                            output_formats=["pdf","png"]):
  """ Create the reweighting plot for one observable.
  """
  fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True, figsize=(10,10))
  ax_up, ax_down = axs
  
  x = xy_dict["SM"][0]
  n_bins = obs_range[0]
  hist_range=(obs_range[1],obs_range[2])
  
  # Upper plot: Relative deviation wrt. SM
  add_rel_dev_plot(ax_up, x, xy_dict, n_bins, hist_range, lumi)
  
  ax_up.set_ylabel(r"$\frac{\# weighted - \# SM}{\# SM} [\%]$", fontsize=30)
  
  # Lower plot: Significance of deviation (in case it would be measured)
  add_dev_sig_plot(ax_down, x, xy_dict, n_bins, hist_range, MC_norm)  
  
  ax_down.set_ylabel(r"$\frac{\# weighted - \# SM}{\sqrt{\# weighted}}$", fontsize=30)
  ax_down.legend(ncol=3, fontsize=16, title=r"$\delta={}$, $L={}$ab$^{{-1}}$".format(dev_scale, lumi/1000), title_fontsize=16)
  ax_down.set_xlim(hist_range)
  ax_down.set_ylim(0, ax_down.get_ylim()[1])
  ax_down.set_xlabel(obs_name)

  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    fig.savefig("{}/{}_ReweightCheck.{}".format(format_dir, obs_name, format))

  plt.close(fig)


def check_reweighting(root_file, tgc_config_path, tgc_point_path, output_dir,
                      observables, mu_charge, lumi, 
                      tree_name = "WWObservables"):
  """ Check what reweighting does to the observable distributions for the given 
      events in the ROOT file.
  """
  rdf = PRH.skip_0weight(ROOT.RDataFrame(tree_name, root_file))
  mu_rdf = PRH.select_mu(rdf, mu_charge)
  tcr = ITCR.TGCConfigReader(tgc_config_path, tgc_point_path)
  
  chirality = IFH.find_chirality(root_file)
  output_base = "{}/ReweightCheck/mu{}/{}/".format(
                  output_dir, PN.sign_str(mu_charge,spelled=True), chirality)
  
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
    hist_dict[obs_name]["SM"] = mu_rdf.Histo1D(h_SM_def, obs_name)
    
    # Histograms with TGC deviation
    for point, index in index_dict.items():
      i_str = str(index+1) # weight naming starts with index 1
      h_def = (obs_name+"_w"+i_str, obs_name+"_w"+i_str, 
               obs_range[0], obs_range[1], obs_range[2])
      weight_name = "rescan_weights.weight"+i_str
      hist_dict[obs_name][point] = mu_rdf.Histo1D(h_def, obs_name, weight_name)
    
  # Get the SM cross section, needed for correct uncertainty weighting
  SM_xsection = mu_rdf.Mean("cross_section")
      
  # Calculate the normalisation of the MC events
  n_MC = rdf.Count() # Use the original rdf with all MC events (except 0-weight)
  MC_norm = SM_xsection.GetValue() * lumi / n_MC.GetValue()
      
  # Make the plots for each observable
  for obs_name, obs_range in observables.items():
    # First transform hists to arrays for easier matplotlib plotting
    xy_dict = {}
    for point, hist in hist_dict[obs_name].items():
      xy_dict[point] = PRHH.TH1_to_arrays(hist) 
    create_reweighting_plot(xy_dict, output_base, obs_name, obs_range, lumi, 
                            MC_norm, dev_scale)
      
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
  
  check_reweighting(RL_path, tgc_config_path, tgc_point_path, output_dir, observables, +1, lumi)
  check_reweighting(RL_path, tgc_config_path, tgc_point_path, output_dir, observables, -1, lumi)
  check_reweighting(LR_path, tgc_config_path, tgc_point_path, output_dir, observables, +1, lumi)
  check_reweighting(LR_path, tgc_config_path, tgc_point_path, output_dir, observables, -1, lumi)
  

if __name__ == "__main__":
  main()