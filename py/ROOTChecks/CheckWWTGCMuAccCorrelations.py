import functools
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
import Plotting.ROOTHistHelp as PRHH
import Systematics.MuonAcceptance as SMA

TGC_name_map = {
  'g1z' : r"g_{{1}}^{{Z}}",
  'ka' : r"\kappa_{{\gamma}}",
  'la' : r"\lambda_{{\gamma}}"
}

def ratio(a,b):
  """ Calculate ratio between two arrays, if denominator point is 0 set ratio=1.
  """
  return np.where(np.abs(b)>0.01,a/b,a)

def TGC_lin_coef(y_TGC, y_SM, dev_scale):
  """ Calculate the normalised gradient wrt. to the TGC at the SM point.
      y_TGC represents the distribution for a change of dev_scale in one TGC.
      Assumes only linear dependence on the TGC.
  """
  return (ratio(y_TGC,y_SM)-1.0) / dev_scale


def add_derivative_plot(ax, TGC_name, xy_dict, n_bins, hist_range, dev_scale):
  """ Add the plots showing the derivative of the cross section in each bin wrt.
      the TGC.
  """
  n_MC = np.sum(xy_dict["no cut"]["SM"][1])
  TGC_grad = functools.partial(TGC_lin_coef, dev_scale=dev_scale)
  
  x = xy_dict["no cut"]["SM"][0]
  d_P_nocut = TGC_grad(xy_dict["no cut"][TGC_name][1], xy_dict["no cut"]["SM"][1])
  d_P_cut = TGC_grad(xy_dict["cut"][TGC_name][1], xy_dict["cut"]["SM"][1])
  d_P_dc = TGC_grad(xy_dict["c shift"][TGC_name][1], xy_dict["c shift"]["SM"][1])
  d_P_dw = TGC_grad(xy_dict["w shift"][TGC_name][1], xy_dict["w shift"]["SM"][1])
  
  ax.hist(x, weights=d_P_nocut, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="no cut")
  ax.hist(x, weights=d_P_cut, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="7deg-cut")
  ax.hist(x, weights=d_P_dc, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$\Delta c += \delta_{cut}$")
  ax.hist(x, weights=d_P_dw, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$\Delta w += \delta_{cut}$")

def add_rel_change_plot(ax, TGC_name, xy_dict, n_bins, hist_range, dev_scale):
  """ Add the plots showing the relative change from the SM for different cut 
      values.
  """
  n_MC = np.sum(xy_dict["no cut"]["SM"][1])
  TGC_grad = functools.partial(TGC_lin_coef, dev_scale=dev_scale)

  x = xy_dict["no cut"]["SM"][0]
  grad_P_nocut = TGC_grad(xy_dict["no cut"][TGC_name][1], xy_dict["no cut"]["SM"][1])
  grad_P_cut = TGC_grad(xy_dict["cut"][TGC_name][1], xy_dict["cut"]["SM"][1])
  grad_P_dc = TGC_grad(xy_dict["c shift"][TGC_name][1], xy_dict["c shift"]["SM"][1])
  grad_P_dw = TGC_grad(xy_dict["w shift"][TGC_name][1], xy_dict["w shift"]["SM"][1])

  rel_P_cut = ratio(grad_P_cut, grad_P_nocut) - 1.0
  rel_P_dc = ratio(grad_P_dc, grad_P_nocut) - 1.0
  rel_P_dw = ratio(grad_P_dw, grad_P_nocut) - 1.0

  scale=100.0
  ax.hist([], alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=rel_P_cut*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=rel_P_dc*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=rel_P_dw*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')

  
def create_reweighting_plot(xy_dict, output_base, obs_name, obs_range, 
                            dev_scale, cut_dev, output_formats=["pdf","png"]):
  """ Create the reweighting plot for one observable.
  """
  eM_chi, eP_chi = PN.chirality_str(IFH.find_chirality(output_base))
  
  for TGC_name in ["g1z", "ka", "la"]:
    fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True, figsize=(10,10))
    ax_up, ax_down = axs
    
    n_bins = obs_range[0]
    hist_range=(obs_range[1],obs_range[2])
    
    # Upper plot
    add_derivative_plot(ax_up, TGC_name, xy_dict, n_bins, hist_range, dev_scale)
    
    ax_up.set_ylabel(r"$\frac{{1}}{{\sigma^{{bin}}}}\frac{{d\sigma^{{bin}}}}{{d {}}} [\frac{{fb}}{{1}}]$".format(TGC_name_map[TGC_name]), fontsize=30)
    ax_up.set_xlim(hist_range)

    ax_up.legend(fontsize=16, title="${}{}$, $\delta_{{cut}}={}$".format(eM_chi, eP_chi, dev_scale, cut_dev), title_fontsize=16)
    
    # Lower plot
    add_rel_change_plot(ax_down, TGC_name, xy_dict, n_bins, hist_range, dev_scale)  
    ax_down.set_ylabel(r"$\frac{y_{cut}-y_{nocut}}{y_{nocut}}$ [%]")
    ax_down.set_xlabel(obs_name)
    

    for format in output_formats:
      format_dir = "{}/{}".format(output_base,format)
      ISH.create_dir(format_dir)
      fig.savefig("{}/{}_{}_TGCMuAccCorr.{}".format(format_dir, obs_name, 
                                                    TGC_name, format))

    plt.close(fig)

def check_reweighting(root_file, tgc_config_path, tgc_point_path, output_dir,
                      observables, tree_name = "WWObservables"):
  """ Check what reweighting does to the observable distributions for the given 
      events in the ROOT file.
  """
  rdf = ROOT.RDataFrame(tree_name, root_file).Filter("rescan_weights.weight1 > 0.01")
  
  # Muon Acceptance related setup
  cut_val = 0.9925 # 7deg
  cut_dev = 5.0e-4
  costh_branch="costh_l"
  
  mu_acc_rdf_dict = {
    "no cut": rdf,
    "cut": rdf.Filter(SMA.get_costh_cut(cut_val, 0, 0, costh_branch)),
    "c shift": rdf.Filter(SMA.get_costh_cut(cut_val, cut_dev, 0, costh_branch)),
    "w shift": rdf.Filter(SMA.get_costh_cut(cut_val, 0, cut_dev, costh_branch))
  }
  
  # TGC related setup
  tcr = ITCR.TGCConfigReader(tgc_config_path, tgc_point_path)
  
  chirality = IFH.find_chirality(root_file)
  output_base = "{}/TGCMuAccCorrelations/{}/".format(output_dir,chirality)
  
  index_dict = {
    "g1z": tcr.point_index([1,0,0]),
    "ka": tcr.point_index([0,1,0]),
    "la": tcr.point_index([0,0,1]),
  }
  dev_scale = tcr.scale
  
  # Dict that stores histograms for each observables and each dev point
  hist_dict = {}
  
  # Book all the histograms as lazy actions of the RDataFrame
  for obs_name, obs_range in observables.items():
    hist_dict[obs_name] = {}
    
    for cut_name, cut_rdf in mu_acc_rdf_dict.items():
      hist_dict[obs_name][cut_name] = {}
      
      # Standard model histograms for comparison
      h_SM_def = (obs_name+"SM", obs_name+"SM", 
                  obs_range[0], obs_range[1], obs_range[2])
      hist_dict[obs_name][cut_name]["SM"] = cut_rdf.Histo1D(h_SM_def, obs_name)
      
      # Histograms with TGC deviation
      for point, index in index_dict.items():
        i_str = str(index+1) # weight naming starts with index 1
        h_def = (obs_name+"_w"+i_str, obs_name+"_w"+i_str, 
                 obs_range[0], obs_range[1], obs_range[2])
        weight_name = "rescan_weights.weight"+i_str
        hist_dict[obs_name][cut_name][point] = cut_rdf.Histo1D(h_def, obs_name, 
                                                               weight_name)
    
  # Make the plots for each observable
  for obs_name, obs_range in observables.items():
    # First transform hists to arrays for easier matplotlib plotting
    xy_dict = {}
    for cut_name, cut_dict in hist_dict[obs_name].items():
      xy_dict[cut_name] = {}
      for point, hist in cut_dict.items():
        xy_dict[cut_name][point] = PRHH.TH1_to_arrays(hist) 
    create_reweighting_plot(xy_dict, output_base, obs_name, obs_range, 
                            dev_scale, cut_dev)
      
def main():
  """ Create check plots that investigate how the weights in the given ROOT file
      affect the distribution.
  """
  ROOT.EnableImplicitMT(7) # Enable multithreading in RDataFrame
  ROOT.gROOT.SetBatch(True) # Don't show graphics at runtime
  PDF.set_default_mpl_format()
  
  config_dir = "/afs/desy.de/group/flc/pool/beyerjac/TGCAnalysis/SampleProduction/MCProduction/PrEWSampleProduction/scripts/config"
  tgc_config_path = config_dir + "/tgc.config"
  tgc_point_path = config_dir + "/tgc_dev_points_g1z_ka_la.config"
  
  output_dir = "../../output"
  
  observables = { "costh_Wminus_star" : (30, -1., 1.),
                  "costh_l_star" :      (30, -1., 1.),
                  "phi_l_star":         (30, -np.pi, np.pi) }
  
  RL_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eR_pL.root"
  LR_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eL_pR.root"
  
  check_reweighting(RL_path, tgc_config_path, tgc_point_path, output_dir, observables)
  check_reweighting(LR_path, tgc_config_path, tgc_point_path, output_dir, observables)
  

if __name__ == "__main__":
  main()