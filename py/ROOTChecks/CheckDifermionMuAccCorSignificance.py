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

def ratio(a,b,default=0.0):
  """ Calculate ratio between two arrays, if denominator point is 0 set default.
  """
  return np.where(np.abs(b)>0.01,a/b,default)

def chisq_delta(y_SM_cut, y_TGC_shift, MC_norm):
  """ Calculate the chi-squared that describes the significance of the given 
      cobined shift in the muon acceptance and difermion parameters.
      MC_norm : normalisation of the Monte Carlo events (needs to include lumi)
  """
  return MC_norm * np.sum( 
    ( ratio(y_TGC_shift - y_SM_cut, np.sqrt(y_SM_cut)) )**2 )

def chisq_cor(y_SM_cut, y_SM_shift, y_TGC_cut, y_TGC_shift, MC_norm):
  """ Calculate the chi-squared that describes the significance of the deviation
      of the factorized approach for the difermion parametrisation & muon 
      acceptance from the actual combined application of both on the MC.
      MC_norm : normalisation of the Monte Carlo events (needs to include lumi)
  """
  return MC_norm * np.sum( 
    ( ratio( (ratio(y_TGC_cut, y_SM_cut) * y_SM_shift  - y_TGC_shift)**2,
             y_SM_cut) ) )

def add_cor_sig_plot(ax, y_dict, MC_norm):
  """ Add a scatter plot comparing the chi-squared that results from treating
      the difermion parameter deviations and the muon acceptance variations as 
      independent (factorizing), to the chi-squared that a true shift in both of 
      them causes.
  """
  chisq_delta_arr = []
  chisq_cor_arr = []
  
  y_SM_cut = y_dict["cut"]["SM"]
  for dev_name in y_dict["no cut"].keys():
    if dev_name == "SM":
      continue # Only check shifts
      
    y_dev_cut = y_dict["cut"][dev_name]
    for cut in y_dict.keys():
      if cut in ["no cut", "cut"]:
        continue # Only check shifts
      
      y_SM_shift = y_dict[cut]["SM"]
      y_dev_shift = y_dict[cut][dev_name]
      
      chisq_delta_arr.append(chisq_delta(y_SM_cut, y_dev_shift, MC_norm))
      chisq_cor_arr.append(chisq_cor(y_SM_cut, y_SM_shift, y_dev_cut, 
                                     y_dev_shift, MC_norm))
                                     
  chisq_delta_arr = np.array(chisq_delta_arr)
  chisq_cor_arr = np.array(chisq_cor_arr)

  label = "Combined 2f par. &\n $\mu$ acc. shifts"
  scatter = ax.scatter(chisq_delta_arr, chisq_cor_arr, label=label)
  
def create_cor_sig_plot(y_dict, output_base, SM_difpar_set, dev_difpar_sets, 
                        c_cut_dev, w_cut_dev, MC_norm, process_str, 
                        output_formats=["pdf","png"]):
  """ Create the correlation significance plot for one observable.
  """
  fig, ax = plt.subplots(tight_layout=True, figsize=(8,6))
  
  add_cor_sig_plot(ax, y_dict, MC_norm)
  
  ax.set_xlabel(r"$\chi_{combined\,shift}^{2}$")
  ax.set_ylabel(r"$\chi_{factorization\,error}^{2}$")
  
  ax.set_xscale('log')
  ax.set_ylim(0, ax.get_ylim()[1])
  
  # legend_title = "${}$\n$\delta_{{TGC}} = {},$\n$\delta_{{c}}^{{\mu-acc.}} = {},$\n$\delta_{{w}}^{{\mu-acc.}} = {}$".format(process_str, dev_scale, c_cut_dev, w_cut_dev)
  ax.legend(loc='upper left')#title=legend_title, )

  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    fig.savefig("{}/CorSig.{}".format(format_dir, format))

  plt.close(fig)

def check_cor_sig(root_file, output_dir, observables, label, mass_range, 
                  Z_direction, lumi, tree_name = "DifermionObservables"):
  """ Check if the muon acceptance and the difermion parameter deviations can be
      factorized for the given events in the ROOT file.
  """
  # Initial RDataFrame with all events
  rdf = ROOT.RDataFrame(tree_name, root_file)
  
  chirality = IFH.find_chirality(root_file)
  output_base = "{}/DifermionMuAccCorSignificance/{}/{}/".format(
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
  c_cut_dev = 2.0e-5
  w_cut_dev = 4.0e-5
  costh_cut_branches = ["costh_f","costh_fbar"]
  
  mu_acc_rdf_dict = {
    "no cut": mu_rdf,
    "cut": mu_rdf.Filter(SMA.get_ndim_costh_cut(cut_val, 0, 0, costh_cut_branches)),
    "c shift +": mu_rdf.Filter(SMA.get_ndim_costh_cut(cut_val, +c_cut_dev, 0, costh_cut_branches)),
    "c shift -": mu_rdf.Filter(SMA.get_ndim_costh_cut(cut_val, -c_cut_dev, 0, costh_cut_branches)),
    "w shift +": mu_rdf.Filter(SMA.get_ndim_costh_cut(cut_val, 0, +w_cut_dev, costh_cut_branches)),
    "w shift -": mu_rdf.Filter(SMA.get_ndim_costh_cut(cut_val, 0, -w_cut_dev, costh_cut_branches)),
  }
  
  # Dict that stores histograms for each observables and each dev point
  hist_dict = {}
  
  obs = observables["costh_f_star"]
  
  # Book all the histograms as lazy actions of the RDataFrame
  for cut_name, cut_rdf in mu_acc_rdf_dict.items():
    hist_dict[cut_name] = {}
    
    # Standard model histograms for comparison
    h_SM_name = "SM_{}".format(cut_name)
    h_SM_def = (h_SM_name, h_SM_name, obs[0], obs[1], obs[2])
                
    # Book all the histograms as lazy actions of the RDataFrame
    hist_dict[cut_name]["SM"] = cut_rdf.Histo1D(
      h_SM_def, "costh_f_star")
    
    # Histograms with TGC deviation
    for dev_name in dev_difpar_sets:
      h_def = (dev_name, dev_name, obs[0], obs[1], obs[2])
      weight_name = "w_"+dev_name
      hist_dict[cut_name][dev_name] = cut_rdf.Histo1D(h_def, "costh_f_star", 
                                                   weight_name)
  
  # Book finding the cross section
  SM_xsection = mu_rdf.Mean("cross_section")
  
  # Calculate the normalisation of the MC events
  n_MC = rdf.Count() # Use the original rdf with all MC events (except 0-weight)
  MC_norm = SM_xsection.GetValue() * lumi / n_MC.GetValue()
  
  # First transform hists to arrays for easier matplotlib plotting
  y_dict = {}
  for cut_name, cut_dict in hist_dict.items():
    y_dict[cut_name] = {}
    for dev_name, hist in cut_dict.items():
      y_dict[cut_name][dev_name] = PRHH.TH1_to_arrays(hist)[1]
  
  create_cor_sig_plot(y_dict, output_base, SM_difpar_set, dev_difpar_sets, 
                      c_cut_dev, w_cut_dev, MC_norm, process_str)
      
def main():
  ROOT.EnableImplicitMT(7) # Enable multithreading in RDataFrame
  ROOT.gROOT.SetBatch(True) # Don't show graphics at runtime
  PDF.set_default_mpl_format()
  
  output_dir = "../../output"
  
  lumi = 1000. # 2ab^-1 total ~> 1ab^-1 for each of the two allowed chiral states
  
  # Important here to use the realistic 3D binning
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
        check_cor_sig(RL_path, output_dir, observables, m_label, m_range, Z_direction, lumi)
        check_cor_sig(LR_path, output_dir, observables, m_label, m_range, Z_direction, lumi)
    else:
      check_cor_sig(RL_path, output_dir, observables, m_label, m_range, None, lumi)
      check_cor_sig(LR_path, output_dir, observables, m_label, m_range, None, lumi)
  

if __name__ == "__main__":
  main()