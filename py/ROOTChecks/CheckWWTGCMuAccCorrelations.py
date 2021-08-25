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
import Plotting.RDFHelp as PRH
import Plotting.ROOTHistHelp as PRHH
import Systematics.MuonAcceptance as SMA

TGC_name_map = {
  'g1z' : r"g_{{1}}^{{Z}}",
  'ka' : r"\kappa_{{\gamma}}",
  'la' : r"\lambda_{{\gamma}}"
}

cut_name_map = {
  'c_shift': r"\Delta c",
  'w_shift': r"\Delta w"
}

def ratio(a,b):
  """ Calculate ratio between two arrays, if denominator point is 0 set ratio=1.
  """
  return np.where(np.abs(b)>0.01,a/b,1.0)

def lin_coef(y_delta, y_zero, delta):
  """ Calculate the normalised gradient:
        1/y_zero * (y_delta-y_zero)/delta
  """
  return (ratio(y_delta,y_zero)-1.0) / delta


def add_TGC_derivative_plot(ax, TGC_name, xy_dict, n_bins, hist_range, dev_scale):
  """ Add the plots showing the derivative of the cross section in each bin wrt.
      the TGC.
  """
  grad = functools.partial(lin_coef, delta=dev_scale)
  
  x = xy_dict["no cut"]["SM"][0]
  d_P_nocut = grad(xy_dict["no cut"][TGC_name][1], xy_dict["no cut"]["SM"][1])
  d_P_cut = grad(xy_dict["cut"][TGC_name][1], xy_dict["cut"]["SM"][1])
  d_P_dc = grad(xy_dict["c_shift"][TGC_name][1], xy_dict["c_shift"]["SM"][1])
  d_P_dw = grad(xy_dict["w_shift"][TGC_name][1], xy_dict["w_shift"]["SM"][1])
  
  ax.hist(x, weights=d_P_nocut, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="no cut")
  ax.hist(x, weights=d_P_cut, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="7deg-cut")
  ax.hist(x, weights=d_P_dc, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$\Delta c $+=$ \delta_{cut}$")
  ax.hist(x, weights=d_P_dw, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$\Delta w $+=$ \delta_{cut}$")

def add_TGC_diff_plot(ax, TGC_name, xy_dict, n_bins, hist_range, dev_scale):
  """ Add the plots showing the relative change from the SM for different cut 
      values.
  """
  grad = functools.partial(lin_coef, delta=dev_scale)

  x = xy_dict["no cut"]["SM"][0]
  grad_P_nocut = grad(xy_dict["no cut"][TGC_name][1], xy_dict["no cut"]["SM"][1])
  grad_P_cut = grad(xy_dict["cut"][TGC_name][1], xy_dict["cut"]["SM"][1])
  grad_P_dc = grad(xy_dict["c_shift"][TGC_name][1], xy_dict["c_shift"]["SM"][1])
  grad_P_dw = grad(xy_dict["w_shift"][TGC_name][1], xy_dict["w_shift"]["SM"][1])

  diff_P_cut = grad_P_cut - grad_P_nocut
  diff_P_dc = grad_P_dc - grad_P_nocut
  diff_P_dw = grad_P_dw - grad_P_nocut

  scale=1.0
  ax.hist([], alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_P_cut*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_P_dc*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_P_dw*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')

def add_cut_derivative_plot(ax, cut_name, xy_dict, n_bins, hist_range, dev_scale):
  """ Add the plots showing the derivative of the cross section in each bin wrt.
      the muon acceptance parameter shift.
  """
  grad = functools.partial(lin_coef, delta=dev_scale)
  
  x = xy_dict["no cut"]["SM"][0]
  d_P_SM = grad(xy_dict[cut_name]["SM"][1], xy_dict["cut"]["SM"][1])
  d_P_g1z = grad(xy_dict[cut_name]["g1z"][1], xy_dict["cut"]["g1z"][1])
  d_P_ka = grad(xy_dict[cut_name]["ka"][1], xy_dict["cut"]["ka"][1])
  d_P_la = grad(xy_dict[cut_name]["la"][1], xy_dict["cut"]["la"][1])
  
  ax.hist(x, weights=d_P_SM, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="SM")
  ax.hist(x, weights=d_P_g1z, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$g_{1}^{Z}$+=$ \delta_{TGC}$")
  ax.hist(x, weights=d_P_ka, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$\kappa_{\gamma}$+=$ \delta_{TGC}$")
  ax.hist(x, weights=d_P_la, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$\lambda_{\gamma}$+=$ \delta_{TGC}$")

def add_cut_diff_plot(ax, cut_name, xy_dict, n_bins, hist_range, dev_scale):
  """ Add the plots showing the relative change from the base cut for different 
      TGC values.
  """
  grad = functools.partial(lin_coef, delta=dev_scale)

  x = xy_dict["no cut"]["SM"][0]
  grad_P_SM = grad(xy_dict[cut_name]["SM"][1], xy_dict["cut"]["SM"][1])
  grad_P_g1z = grad(xy_dict[cut_name]["g1z"][1], xy_dict["cut"]["g1z"][1])
  grad_P_ka = grad(xy_dict[cut_name]["ka"][1], xy_dict["cut"]["ka"][1])
  grad_P_la = grad(xy_dict[cut_name]["la"][1], xy_dict["cut"]["la"][1])

  diff_P_g1z = grad_P_g1z - grad_P_SM
  diff_P_ka = grad_P_ka - grad_P_SM
  diff_P_la = grad_P_la - grad_P_SM

  scale=1.0
  ax.hist([], alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_P_g1z*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_P_ka*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_P_la*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')

  
def create_reweighting_plots(xy_dict, output_base, obs_name, obs_range, 
                             dev_scale, cut_dev, process_str,
                             output_formats=["pdf","png"]):
  """ Create the reweighting plot for one observable.
  """
  for TGC_name in ["g1z", "ka", "la"]:
    fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True, figsize=(10,10))
    ax_up, ax_down = axs
    
    n_bins = obs_range[0]
    hist_range=(obs_range[1],obs_range[2])
    
    # Upper plot
    add_TGC_derivative_plot(ax_up, TGC_name, xy_dict, n_bins, hist_range, dev_scale)
    
    ax_up.set_ylabel(r"$y=\frac{{1}}{{\sigma^{{bin}}}}\frac{{d\sigma^{{bin}}}}{{d {}}}$".format(TGC_name_map[TGC_name]), fontsize=30)
    ax_up.set_xlim(hist_range)

    ax_up.legend(fontsize=16, title="${}$, $\delta_{{cut}}={}$".format(process_str, cut_dev), title_fontsize=16)
    
    # Lower plot
    add_TGC_diff_plot(ax_down, TGC_name, xy_dict, n_bins, hist_range, dev_scale)  
    ax_down.set_ylabel(r"$y_{cut}-y_{nocut}$")
    ax_down.set_xlabel(obs_name)
    
    for format in output_formats:
      format_dir = "{}/{}".format(output_base,format)
      ISH.create_dir(format_dir)
      fig.savefig("{}/{}_{}_TGCMuAccCorr.{}".format(format_dir, obs_name, 
                                                    TGC_name, format))

    plt.close(fig)
    
  for cut_name in ["c_shift", "w_shift"]:
    fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True, figsize=(10,10))
    ax_up, ax_down = axs
    
    n_bins = obs_range[0]
    hist_range=(obs_range[1],obs_range[2])
    
    # Upper plot
    add_cut_derivative_plot(ax_up, cut_name, xy_dict, n_bins, hist_range, dev_scale)
    
    ax_up.set_ylabel(r"$y=\frac{{1}}{{\sigma^{{bin}}}}\frac{{d\sigma^{{bin}}}}{{d {}}}$".format(cut_name_map[cut_name]), fontsize=30)
    ax_up.set_xlim(hist_range)

    ax_up.legend(fontsize=16, title="${}$, $\delta_{{TGC}}={}$".format(process_str, dev_scale), title_fontsize=16)
    
    # Lower plot
    add_cut_diff_plot(ax_down, cut_name, xy_dict, n_bins, hist_range, dev_scale)  
    ax_down.set_ylabel(r"$y_{TGC}-y_{SM}$")
    ax_down.set_xlabel(obs_name)
    
    for format in output_formats:
      format_dir = "{}/{}".format(output_base,format)
      ISH.create_dir(format_dir)
      fig.savefig("{}/{}_{}_TGCMuAccCorr.{}".format(format_dir, obs_name, 
                                                    cut_name, format))

    plt.close(fig)

def check_reweighting(root_file, tgc_config_path, tgc_point_path, output_dir,
                      observables, mu_charge, tree_name = "WWObservables"):
  """ Check what reweighting does to the observable distributions for the given 
      events in the ROOT file.
  """
  rdf = PRH.skip_0weight(ROOT.RDataFrame(tree_name, root_file))
  mu_rdf = PRH.select_mu(rdf, mu_charge)
  
  # Muon Acceptance related setup
  cut_val = 0.9925 # 7deg
  cut_dev = 5.0e-4
  costh_branch="costh_l"
  
  mu_acc_rdf_dict = {
    "no cut": mu_rdf,
    "cut": mu_rdf.Filter(SMA.get_costh_cut(cut_val, 0, 0, costh_branch)),
    "c_shift": mu_rdf.Filter(SMA.get_costh_cut(cut_val, cut_dev, 0, costh_branch)),
    "w_shift": mu_rdf.Filter(SMA.get_costh_cut(cut_val, 0, cut_dev, costh_branch))
  }
  
  # TGC related setup
  tcr = ITCR.TGCConfigReader(tgc_config_path, tgc_point_path)
  
  chirality = IFH.find_chirality(root_file)
  output_base = "{}/TGCMuAccCorrelations/mu{}/{}/".format(
                  output_dir, PN.sign_str(mu_charge,spelled=True), chirality)
  process_str = PN.process_str(chirality, mu_charge)
  
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
    create_reweighting_plots(xy_dict, output_base, obs_name, obs_range, 
                             dev_scale, cut_dev, process_str)
      
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
  
  check_reweighting(RL_path, tgc_config_path, tgc_point_path, output_dir, observables, +1)
  check_reweighting(RL_path, tgc_config_path, tgc_point_path, output_dir, observables, -1)
  check_reweighting(LR_path, tgc_config_path, tgc_point_path, output_dir, observables, +1)
  check_reweighting(LR_path, tgc_config_path, tgc_point_path, output_dir, observables, -1)
  

if __name__ == "__main__":
  main()