import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import ROOT

# Use the the local helper modules
import sys
sys.path.append("../Helpers") 
import IO.FilenameHelp as IFH
import IO.SysHelpers as ISH
import Plotting.DefaultFormat as PDF
import Plotting.Naming as PN
import Plotting.RDFHelp as PRH
import Plotting.ROOTHistHelp as PRHH

def create_comparison_plot(x, y_cut, y_nocut, output_base, 
                           obs_name, obs_range, output_formats=["pdf","png"]):
  """ Create the comparison plot the given observable, given the values with and
      without the cut.
  """
  fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True, figsize=(8,10))
  ax_up, ax_down = axs
  
  hist_range=(obs_range[1],obs_range[2])
  ax_up.hist(x, weights=y_nocut, range=hist_range, bins=obs_range[0], histtype=u'step', label="no cut")
  ax_up.hist(x, weights=y_cut, range=hist_range, bins=obs_range[0], histtype=u'step', label="0-weight cut")
  ax_up.set_ylabel("# MC events")
  ax_up.set_yscale('log')
  ax_up.legend()
  ax_up.set_xlim(hist_range)
  
  # Calculate the ratio and the error on the ratio
  delta_y = y_cut-y_nocut
  r = delta_y/y_nocut
  r_err = np.sqrt(y_cut * np.abs(delta_y)) / y_nocut**(3/2)
  r_mean = np.mean(r)
  
  scale=100. # in percent
  ax_down.plot(x, r * scale)
  ax_down.fill_between(x, (r-r_err)*scale, (r+r_err)*scale, alpha=0.5)
  ax_down.axline((np.amin(x), r_mean*scale), (np.amax(x), r_mean*scale), ls='--', color='black')
  ax_down.set_ylabel(r"$\frac{\# cut - \# nocut}{\# nocut} [\%]$")
  ax_down.set_xlabel(obs_name)

  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    fig.savefig("{}/{}_ZeroWeightCut.{}".format(format_dir, obs_name, format))

  plt.close(fig)

def check_0weight_cut(root_file, output_dir, observables, mu_charge,
                      tree_name = "WWObservables"):
  """ Create all 0-weight cut comparisons for the given observables on this ROOT
      file.
  """
  rdf = PRH.select_mu(ROOT.RDataFrame(tree_name, root_file), mu_charge)
  rdf_no_0_w1 = PRH.skip_0weight(rdf)
  
  chirality = IFH.find_chirality(root_file)
  output_base = "{}/ZeroWeightCut/mu{}/{}/".format(
                  output_dir, PN.sign_str(mu_charge,spelled=True), chirality)
  
  # First "lazy" definition of histograms, cummulating the tasks that are to be 
  # performed on the rdf
  obs_hists = {}
  for obs_name, obs_range in observables.items():
    h_nocut_def = ("nocut_"+obs_name, "nocut_"+obs_name, 
                   obs_range[0], obs_range[1], obs_range[2])
    h_cut_def = ("cut_"+obs_name, "cut_"+obs_name, 
                 obs_range[0], obs_range[1], obs_range[2])
    obs_hists[obs_name] = [ rdf.Histo1D(h_nocut_def, obs_name),
                            rdf_no_0_w1.Histo1D(h_cut_def, obs_name) ]
    
  # Now start plotting
  for obs_name in observables:
    x_nocut, y_nocut = PRHH.TH1_to_arrays(obs_hists[obs_name][0])
    x_cut, y_cut = PRHH.TH1_to_arrays(obs_hists[obs_name][1])
  
    create_comparison_plot(x_nocut, y_cut, y_nocut, output_base, 
                           obs_name, observables[obs_name])
  
def main():
  """ Check how cutting away the events where the reweighting failed (and gave 0
      weights) affects the distributions.
  """
  ROOT.EnableImplicitMT() # Enable multithreading in RDataFrame
  ROOT.gROOT.SetBatch(True) # Don't show graphics at runtime
  PDF.set_default_mpl_format()
  
  output_dir = "../../output"
  
  observables = { "costh_Wminus_star" : (30, -1., 1.),
                  "costh_l_star" :      (30, -1., 1.),
                  "phi_l_star":         (30, -np.pi, np.pi) }
  
  LR_file = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eL_pR.root"
  RL_file = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eR_pL.root"

  check_0weight_cut(RL_file, output_dir, observables, +1)
  check_0weight_cut(RL_file, output_dir, observables, -1)
  check_0weight_cut(LR_file, output_dir, observables, +1)
  check_0weight_cut(LR_file, output_dir, observables, -1)

if __name__ == "__main__":
  main()