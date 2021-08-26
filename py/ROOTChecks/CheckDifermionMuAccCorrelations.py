import functools
import matplotlib.pyplot as plt
import numpy as np
import ROOT

# Use the the local helper modules
import sys
sys.path.append("../Helpers") 
import IO.FilenameHelp as IFH
import IO.SysHelpers as ISH
import Physics.Difermion as PD
import Plotting.DefaultFormat as PDF
import Plotting.Naming as PN
import Plotting.RDFHelp as PRH
import Plotting.ROOTHistHelp as PRHH
import Systematics.MuonAcceptance as SMA

dif_par_name_map = {
  "d_Ae": "A_{{e}}",
  "d_Af": "A_{{f}}",
  "d_ef": "\\epsilon_{{f}}",
  "d_k0": "k_{{0}}",
  "d_dk": "\\Delta k"
}

cut_name_map = {
  'c_shift': r"\Delta c",
  'w_shift': r"\Delta w"
}

def dev_name_index(dev_name):
  return ["d_xs0", "d_Ae", "d_Af", "d_ef", "d_k0", "d_dk"].index(dev_name)

def get_dev_delta(dev_name, SM_difpar_set, dev_difpar_set):
  return (np.array(dev_difpar_set.arr()) - np.array(SM_difpar_set.arr()))[dev_name_index(dev_name)]

def ratio(a,b,default=1.0):
  """ Calculate ratio between two arrays, if denominator point is 0 set ratio=1.
  """
  return np.where(np.abs(b)>0.01,a/b,default)

def lin_coef(y_delta, y_zero, delta):
  """ Calculate the normalised gradient:
        1/y_zero * (y_delta-y_zero)/delta
  """
  return (ratio(y_delta,y_zero)-1.0) / delta


def add_dif_derivative_plot(ax, dev_name, xy_dict, n_bins, hist_range, 
                           SM_difpar_set, difpar_set):
  """ Add the plots showing the derivative of the cross section in each bin wrt.
      the difermion parameters.
  """
  dev_scale = get_dev_delta(dev_name, SM_difpar_set, difpar_set)
  grad = functools.partial(lin_coef, delta=dev_scale)

  x = xy_dict["no cut"]["SM"][0]
  grad_nocut = grad(xy_dict["no cut"][dev_name][1], xy_dict["no cut"]["SM"][1])
  grad_cut = grad(xy_dict["cut"][dev_name][1], xy_dict["cut"]["SM"][1])
  grad_dc = grad(xy_dict["c_shift"][dev_name][1], xy_dict["c_shift"]["SM"][1])
  grad_dw = grad(xy_dict["w_shift"][dev_name][1], xy_dict["w_shift"]["SM"][1])

  ax.hist(x, weights=grad_nocut, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="no cut")
  ax.hist(x, weights=grad_cut, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="7deg-cut")
  ax.hist(x, weights=grad_dc, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$\Delta c $+=$ \delta_{c}$")
  ax.hist(x, weights=grad_dw, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$\Delta w $+=$ \delta_{w}$")

def add_dif_diff_plot(ax, dev_name, xy_dict, n_bins, hist_range, 
                      SM_difpar_set, difpar_set):
  """ Add the plots showing the relative change from the SM for different cut 
      values.
  """
  dev_scale = get_dev_delta(dev_name, SM_difpar_set, difpar_set)
  grad = functools.partial(lin_coef, delta=dev_scale)

  x = xy_dict["no cut"]["SM"][0]
  grad_nocut = grad(xy_dict["no cut"][dev_name][1], xy_dict["no cut"]["SM"][1])
  grad_cut = grad(xy_dict["cut"][dev_name][1], xy_dict["cut"]["SM"][1])
  grad_dc = grad(xy_dict["c_shift"][dev_name][1], xy_dict["c_shift"]["SM"][1])
  grad_dw = grad(xy_dict["w_shift"][dev_name][1], xy_dict["w_shift"]["SM"][1])

  diff_P_cut = grad_cut - grad_nocut
  diff_P_dc = grad_dc - grad_nocut
  diff_P_dw = grad_dw - grad_nocut

  scale=1.0
  ax.hist([], alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_P_cut*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_P_dc*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_P_dw*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')

def add_cut_derivative_plot(ax, cut_name, xy_dict, n_bins, hist_range, 
                            dev_scale, SM_difpar_set, dev_difpar_sets):
  """ Add the plots showing the derivative of the cross section in each bin wrt.
      the muon acceptance parameter shift.
  """
  grad = functools.partial(lin_coef, delta=dev_scale)

  x = xy_dict["no cut"]["SM"][0]
  grad_SM = grad(xy_dict[cut_name]["SM"][1], xy_dict["cut"]["SM"][1])
  grad_d_Ae = grad(xy_dict[cut_name]["d_Ae"][1], xy_dict["cut"]["d_Ae"][1])
  grad_d_Af = grad(xy_dict[cut_name]["d_Af"][1], xy_dict["cut"]["d_Af"][1])
  grad_d_ef = grad(xy_dict[cut_name]["d_ef"][1], xy_dict["cut"]["d_ef"][1])
  grad_d_k0 = grad(xy_dict[cut_name]["d_k0"][1], xy_dict["cut"]["d_k0"][1])
  grad_d_dk = grad(xy_dict[cut_name]["d_dk"][1], xy_dict["cut"]["d_dk"][1])

  delta_Ae = get_dev_delta("d_Ae", SM_difpar_set, dev_difpar_sets["d_Ae"])
  delta_Af = get_dev_delta("d_Af", SM_difpar_set, dev_difpar_sets["d_Af"])
  delta_ef = get_dev_delta("d_ef", SM_difpar_set, dev_difpar_sets["d_ef"])
  delta_k0 = get_dev_delta("d_k0", SM_difpar_set, dev_difpar_sets["d_k0"])
  delta_dk = get_dev_delta("d_dk", SM_difpar_set, dev_difpar_sets["d_dk"])

  ax.hist(x, weights=grad_SM, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="SM")
  ax.hist(x, weights=grad_d_Ae, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$A_{{e}}$+={:.1e}".format(delta_Ae))
  ax.hist(x, weights=grad_d_Af, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$A_{{f}}$+={:.1e}".format(delta_Af))
  ax.hist(x, weights=grad_d_ef, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$\epsilon_{{f}}$+={:.1e}".format(delta_ef))
  ax.hist(x, weights=grad_d_k0, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$k_{{0}}$+={:.1e}".format(delta_k0))
  ax.hist(x, weights=grad_d_dk, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step', label="$\Delta k$+={:.1e}".format(delta_dk))

def add_cut_diff_plot(ax, cut_name, xy_dict, n_bins, hist_range, dev_scale):
  """ Add the plots showing the relative change from the base cut for different 
      difermion parameter values.
  """
  grad = functools.partial(lin_coef, delta=dev_scale)

  grad = functools.partial(lin_coef, delta=dev_scale)

  x = xy_dict["no cut"]["SM"][0]
  grad_SM = grad(xy_dict[cut_name]["SM"][1], xy_dict["cut"]["SM"][1])
  grad_d_Ae = grad(xy_dict[cut_name]["d_Ae"][1], xy_dict["cut"]["d_Ae"][1])
  grad_d_Af = grad(xy_dict[cut_name]["d_Af"][1], xy_dict["cut"]["d_Af"][1])
  grad_d_ef = grad(xy_dict[cut_name]["d_ef"][1], xy_dict["cut"]["d_ef"][1])
  grad_d_k0 = grad(xy_dict[cut_name]["d_k0"][1], xy_dict["cut"]["d_k0"][1])
  grad_d_dk = grad(xy_dict[cut_name]["d_dk"][1], xy_dict["cut"]["d_dk"][1])

  diff_d_Ae = grad_d_Ae - grad_SM
  diff_d_Af = grad_d_Af - grad_SM
  diff_d_ef = grad_d_ef - grad_SM
  diff_d_k0 = grad_d_k0 - grad_SM
  diff_d_dk = grad_d_dk - grad_SM

  scale=1.0
  ax.hist([], alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_d_Ae*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_d_Af*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_d_ef*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_d_k0*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')
  ax.hist(x, weights=diff_d_dk*scale, alpha=0.9, lw=3, range=hist_range, bins=n_bins, histtype=u'step')


def create_reweighting_plots(xy_dict, output_base, obs_name, obs_range, 
                             SM_difpar_set, dev_difpar_sets, c_dev, w_dev, 
                             process_str, output_formats=["pdf","png"]):
  """ Create the reweighting plot for one observable.
  """
  for dev_name, difpar_set in dev_difpar_sets.items():
    fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True, figsize=(10,10))
    ax_up, ax_down = axs

    n_bins = obs_range[0]
    hist_range=(obs_range[1],obs_range[2])

    # Upper plot
    add_dif_derivative_plot(ax_up, dev_name, xy_dict, n_bins, hist_range, SM_difpar_set, difpar_set)

    ax_up.set_ylabel(r"$y=\frac{{1}}{{\sigma^{{bin}}}}\frac{{d\sigma^{{bin}}}}{{d {}}}$".format(dif_par_name_map[dev_name]), fontsize=30)
    ax_up.set_xlim(hist_range)

    ax_up.legend(fontsize=16, title="${}$,\n$\delta_{{c}}={}$, $\delta_{{w}}={}$".format(process_str, c_dev, w_dev), title_fontsize=16)

    # Lower plot
    add_dif_diff_plot(ax_down, dev_name, xy_dict, n_bins, hist_range, SM_difpar_set, difpar_set)
    ax_down.set_ylabel(r"$y_{cut}-y_{nocut}$")
    ax_down.set_xlabel(obs_name)

    for format in output_formats:
      format_dir = "{}/{}".format(output_base,format)
      ISH.create_dir(format_dir)
      fig.savefig("{}/{}_{}_DifMuAccCorr.{}".format(format_dir, obs_name, 
                                                    dev_name, format))

    plt.close(fig)

  for cut_name in ["c_shift", "w_shift"]:
    fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True, figsize=(10,10))
    ax_up, ax_down = axs
  
    n_bins = obs_range[0]
    hist_range=(obs_range[1],obs_range[2])
  
    # Upper plot
    dev_scale = c_dev if (cut_name == "c_shift") else w_dev
    add_cut_derivative_plot(ax_up, cut_name, xy_dict, n_bins, hist_range, dev_scale, SM_difpar_set, dev_difpar_sets)
  
    ax_up.set_ylabel(r"$y=\frac{{1}}{{\sigma^{{bin}}}}\frac{{d\sigma^{{bin}}}}{{d {}}}$".format(cut_name_map[cut_name]), fontsize=30)
    ax_up.set_xlim(hist_range)
    # 
    ax_up.legend(fontsize=16, title="${}$".format(process_str), title_fontsize=16)
    # 
    # # Lower plot
    add_cut_diff_plot(ax_down, cut_name, xy_dict, n_bins, hist_range, dev_scale)  
    ax_down.set_ylabel(r"$y_{deviation}-y_{SM}$")
    ax_down.set_xlabel(obs_name)
  
    for format in output_formats:
      format_dir = "{}/{}".format(output_base,format)
      ISH.create_dir(format_dir)
      fig.savefig("{}/{}_{}_DifMuAccCorr.{}".format(format_dir, obs_name, 
                                                    cut_name, format))
  
    plt.close(fig)

def check_reweighting(root_file, output_dir, observables, label, mass_range, 
                      Z_direction, tree_name = "DifermionObservables"):
  """ Check what reweighting does to the observable distributions for the given 
      events in the ROOT file.
  """
  # Initial RDataFrame with all events
  rdf = ROOT.RDataFrame(tree_name, root_file)
  
  chirality = IFH.find_chirality(root_file)
  output_base = "{}/DifermionMuAccCorrelations/{}/{}/".format(
                  output_dir, 
                  label if not Z_direction else label + "_" + Z_direction, 
                  chirality)
  process_str = PN.difermion_process_str("mu", chirality, label, Z_direction)
  
  # Define the parameter deviations to check
  SM_difpar_set = PD.SM_parset("mu", label)
  dev_difpar_sets = {
    "d_Ae": SM_difpar_set.add_dev(0, 8.e-4, 0, 0, 0, 0),
    "d_Af": SM_difpar_set.add_dev(0, 0, 7.e-4, 0, 0, 0),
    "d_ef": SM_difpar_set.add_dev(0, 0, 0, 1.e-3, 0, 0),
    "d_k0": SM_difpar_set.add_dev(0, 0, 0, 0, 1.5e-3, 0),
    "d_dk": SM_difpar_set.add_dev(0, 0, 0, 0, 0, 1.7e-3)
  }
  
  # Add the weight branches for the different models
  for dev_name, dev_difpar_set in dev_difpar_sets.items():
    # Updates the rdf to include additional functional branches
    rdf = PD.add_weight_branch(rdf, "w_"+dev_name, "costh_f_star", 
                               SM_difpar_set, dev_difpar_set, chirality)
  
  mu_rdf = PRH.select_mumu(rdf, mass_range[0], mass_range[1], Z_direction)

  # Muon Acceptance related setup
  cut_val = 0.9925 # 7deg
  c_dev = 2.0e-5
  w_dev = 4.0e-5
  costh_cut_branches = ["costh_f","costh_fbar"]
  
  mu_acc_rdf_dict = {
    "no cut": mu_rdf,
    "cut": mu_rdf.Filter(SMA.get_ndim_costh_cut(cut_val, 0, 0, costh_cut_branches)),
    "c_shift": mu_rdf.Filter(SMA.get_ndim_costh_cut(cut_val, c_dev, 0, costh_cut_branches)),
    "w_shift": mu_rdf.Filter(SMA.get_ndim_costh_cut(cut_val, 0, w_dev, costh_cut_branches))
  }
  
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
  
      # Histograms with difermion parameter deviations
      for dev_name in dev_difpar_sets:
        h_def = (obs_name+dev_name, obs_name+dev_name, 
                 obs_range[0], obs_range[1], obs_range[2])
        weight_name = "w_"+dev_name
        hist_dict[obs_name][cut_name][dev_name] = cut_rdf.Histo1D(
                                                  h_def, obs_name, weight_name)
  
  # Make the plots for each observable
  for obs_name, obs_range in observables.items():
    # First transform hists to arrays for easier matplotlib plotting
    xy_dict = {}
    for cut_name, cut_dict in hist_dict[obs_name].items():
      xy_dict[cut_name] = {}
      for dev_name, hist in cut_dict.items():
        xy_dict[cut_name][dev_name] = PRHH.TH1_to_arrays(hist) 
    create_reweighting_plots(xy_dict, output_base, obs_name, obs_range, 
                             SM_difpar_set, dev_difpar_sets, c_dev, w_dev, 
                             process_str)
      
def main():
  """ Create check plots that investigate how the weights in the given ROOT file
      affect the distribution.
  """
  ROOT.EnableImplicitMT(7) # Enable multithreading in RDataFrame
  ROOT.gROOT.SetBatch(True) # Don't show graphics at runtime
  PDF.set_default_mpl_format()
  
  output_dir = "../../output"
  
  observables = { "costh_f_star": (20, -1.0, 1.0) }
  
  RL_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/2f_Z_l_eR_pL.root"
  LR_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/2f_Z_l_eL_pR.root"

  energy = 250
  m_ranges = {
    "return_to_Z": [81, 101],
    "high_Q2": [180, 1.1*energy] }
  FB_split = {
    "return_to_Z": True,
    "high_Q2": False }
  
  for m_label, m_range in m_ranges.items():
    if FB_split[m_label]:
      for Z_direction in ["FZ", "BZ"]:
        check_reweighting(RL_path, output_dir, observables, m_label, m_range, Z_direction)
        check_reweighting(LR_path, output_dir, observables, m_label, m_range, Z_direction)
    else:
      check_reweighting(RL_path, output_dir, observables, m_label, m_range, None)
      check_reweighting(LR_path, output_dir, observables, m_label, m_range, None)
    
if __name__ == "__main__":
  main()