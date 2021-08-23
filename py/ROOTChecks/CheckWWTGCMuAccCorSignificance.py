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
import Plotting.Naming as PN
import Plotting.ROOTHistHelp as PRHH
import Systematics.MuonAcceptance as SMA

def ratio(a,b,default=0.0):
  """ Calculate ratio between two arrays, if denominator point is 0 set default.
  """
  return np.where(np.abs(b)>0.01,a/b,default)

def chisq_delta(y_SM_cut, y_TGC_shift, MC_norm):
  """ Calculate the chi-squared that describes the significance of the given 
      cobined shift in the muon acceptance and TGCs.
      MC_norm : normalisation of the Monte Carlo events (needs to include lumi)
  """
  return MC_norm * np.sum( 
    ( ratio(y_TGC_shift - y_SM_cut, np.sqrt(y_SM_cut)) )**2 )

def chisq_cor(y_SM_cut, y_SM_shift, y_TGC_cut, y_TGC_shift, MC_norm):
  """ Calculate the chi-squared that describes the significance of the deviation
      of the factorized approach for the TGC & muon acceptance from the actual
      combined application of both on the MC.
      MC_norm : normalisation of the Monte Carlo events (needs to include lumi)
  """
  return MC_norm * np.sum( 
    ( ratio( ratio(y_SM_shift * y_TGC_cut, y_SM_cut)  - y_TGC_shift,
             np.sqrt(y_SM_cut)) )**2 )

def add_cor_sig_plot(ax, y_dict, MC_norm):
  """ Add a scatter plot comparing the chi-squared that results from treating
      the TGC deviations and the muon acceptance variations as independent 
      (factorizing), to the chi-squared that a true shift in both of them 
      causes.
  """
  chisq_delta_arr = []
  chisq_cor_arr = []
  
  y_SM_cut = y_dict["cut"]["SM"]
  for TGC_name in y_dict["no cut"].keys():
    if TGC_name == "SM":
      continue # Only check shifts
      
    y_TGC_cut = y_dict["cut"][TGC_name]
    for cut in y_dict.keys():
      if cut in ["no cut", "cut"]:
        continue # Only check shifts
      
      y_SM_shift = y_dict[cut]["SM"]
      y_TGC_shift = y_dict[cut][TGC_name]
      
      chisq_delta_arr.append(chisq_delta(y_SM_cut, y_TGC_shift, MC_norm))
      chisq_cor_arr.append(chisq_cor(y_SM_cut, y_SM_shift, y_TGC_cut, 
                                     y_TGC_shift, MC_norm))
                                     
  chisq_delta_arr = np.array(chisq_delta_arr)
  chisq_cor_arr = np.array(chisq_cor_arr)

  label = "Combined TGC &\n $\mu$ acc. shifts"
  scatter = ax.scatter(chisq_delta_arr, chisq_cor_arr, label=label)
  
def create_cor_sig_plot(y_dict, output_base, dev_scale, c_cut_dev, w_cut_dev, 
                        MC_norm, output_formats=["pdf","png"]):
  """ Create the correlation significance plot for one observable.
  """
  eM_chi, eP_chi = PN.chirality_str(IFH.find_chirality(output_base))
  fig, ax = plt.subplots(tight_layout=True, figsize=(8,6))
  
  add_cor_sig_plot(ax, y_dict, MC_norm)
  
  ax.set_xlabel(r"$\chi_{combined\,shift}^{2}$")
  ax.set_ylabel(r"$\chi_{factorization\,error}^{2}$")
  
  ax.set_xscale('log')
  ax.set_ylim(0, ax.get_ylim()[1])
  
  legend_title = "${}{}$\n$\delta_{{TGC}} = {},$\n$\delta_{{c}}^{{\mu-acc.}} = {},$\n$\delta_{{w}}^{{\mu-acc.}} = {}$".format(eM_chi, eP_chi, dev_scale, c_cut_dev, w_cut_dev)
  ax.legend(title=legend_title, loc='upper left')

  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    fig.savefig("{}/CorSig.{}".format(format_dir, format))

  plt.close(fig)

def check_cor_sig(root_file, tgc_config_path, tgc_point_path, output_dir,
                  observables, lumi, tree_name = "WWObservables"):
  """ Check if the muon acceptance and the TGC deviations can be factorized
      for the given events in the ROOT file.
  """
  rdf = ROOT.RDataFrame(tree_name, root_file).Filter("rescan_weights.weight1 > 0.01")
  
  # Muon Acceptance related setup
  cut_val = 0.9925 # 7deg
  c_cut_dev = 5.0e-4 
  w_cut_dev = 1.0e-4
  costh_branch="costh_l"
  
  mu_acc_rdf_dict = {
    "no cut": rdf,
    "cut": rdf.Filter(SMA.get_costh_cut(cut_val, 0, 0, costh_branch)),
    "c shift +": rdf.Filter(SMA.get_costh_cut(cut_val, +c_cut_dev, 0, costh_branch)),
    "c shift -": rdf.Filter(SMA.get_costh_cut(cut_val, -c_cut_dev, 0, costh_branch)),
    "w shift +": rdf.Filter(SMA.get_costh_cut(cut_val, 0, +w_cut_dev, costh_branch)),
    "w shift -": rdf.Filter(SMA.get_costh_cut(cut_val, 0, -w_cut_dev, costh_branch)),
  }
  
  # TGC related setup
  tcr = ITCR.TGCConfigReader(tgc_config_path, tgc_point_path)
  
  chirality = IFH.find_chirality(root_file)
  output_base = "{}/TGCMuAccCorSignificance/{}/".format(output_dir,chirality)
  
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
  
  # Book all the histograms as lazy actions of the RDataFrame
  for cut_name, cut_rdf in mu_acc_rdf_dict.items():
    hist_dict[cut_name] = {}
    
    # Standard model histograms for comparison
    h_SM_name = "SM_{}".format(cut_name)
    h_SM_def = (h_SM_name, h_SM_name, obs_1[0], obs_1[1], obs_1[2],
                                      obs_2[0], obs_2[1], obs_2[2],
                                      obs_3[0], obs_3[1], obs_3[2])
                
    # Book all the histograms as lazy actions of the RDataFrame
    hist_dict[cut_name]["SM"] = cut_rdf.Histo3D(
      h_SM_def, obs_names[0], obs_names[1], obs_names[2])
    
    # Histograms with TGC deviation
    for point, index in index_dict.items():
      i_str = str(index+1) # weight naming starts with index 1
      h_name = "w{}_{}".format(i_str, cut_name) 
      h_def = (h_name, h_name, obs_1[0], obs_1[1], obs_1[2],
                               obs_2[0], obs_2[1], obs_2[2],
                               obs_3[0], obs_3[1], obs_3[2])
      weight_name = "rescan_weights.weight"+i_str
      hist_dict[cut_name][point] = cut_rdf.Histo3D(
        h_def, obs_names[0], obs_names[1], obs_names[2], weight_name)
  
  # Book finding the cross section
  SM_xsection = rdf.Mean("cross_section")
  
  # Make the plots for each observable
  y_dict = {}
  # First transform hists to arrays for easier matplotlib plotting
  y_dict = {}
  for cut_name, cut_dict in hist_dict.items():
    y_dict[cut_name] = {}
    for point, hist in cut_dict.items():
      y_dict[cut_name][point] = PRHH.TH3_to_arrays(hist)[1]
  
  # Calculate the normalisation of the MC events
  n_MC = np.sum(y_dict["no cut"]["SM"])
  MC_norm = SM_xsection.GetValue() * lumi / n_MC
  
  create_cor_sig_plot(y_dict, output_base, dev_scale, c_cut_dev, w_cut_dev, 
                      MC_norm)
      
def main():
  ROOT.EnableImplicitMT(7) # Enable multithreading in RDataFrame
  ROOT.gROOT.SetBatch(True) # Don't show graphics at runtime
  PDF.set_default_mpl_format()
  
  config_dir = "/afs/desy.de/group/flc/pool/beyerjac/TGCAnalysis/SampleProduction/MCProduction/PrEWSampleProduction/scripts/config"
  tgc_config_path = config_dir + "/tgc.config"
  tgc_point_path = config_dir + "/tgc_dev_points_g1z_ka_la.config"
  
  output_dir = "../../output"
  
  lumi = 1000. # 2ab^-1 total ~> 1ab^-1 for each of the two allowed chiral states
  
  # Important here to use the realistic 3D binning
  observables = { "costh_Wminus_star" : (20, -1.0, 1.0),
                  "costh_l_star" :      (10, -1., 1.),
                  "phi_l_star":         (10, -np.pi, np.pi) }
                  
  RL_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eR_pL.root"
  LR_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eL_pR.root"
  
  check_cor_sig(RL_path, tgc_config_path, tgc_point_path, output_dir, observables, lumi)
  check_cor_sig(LR_path, tgc_config_path, tgc_point_path, output_dir, observables, lumi)
  

if __name__ == "__main__":
  main()