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
  
def create_xs0_plot(xs0_dict, output_base, dev_scale, mu_charge, 
                    output_formats=["pdf","png"]):
  """ Check how the total chiral cross section changes when reweighting.
  """
  fig = plt.figure(tight_layout=True, figsize=(9,6.5))
  ax = plt.gca()
  
  x = [0.5,1.5,2.5]
  ax.set_xticks(x, minor=True)
  x_ticks = [ "$g_1^{Z}$", "$\kappa_{\gamma}$", "$\lambda_{\gamma}$" ]
  plt.xticks(x, x_ticks, size='large')
  
  y_SM = xs0_dict["SM"]
  
  y_TGC_P = np.array([xs0_dict[point] for point in ["g1z +", "ka +", "la +"]])
  y_TGC_M = np.array([xs0_dict[point] for point in ["g1z -", "ka -", "la -"]])
  
  ax.plot([0,3], [1.,1.])
  ax.errorbar(x, y_TGC_P / y_SM, xerr=0.5, marker="*", ls='none', label="$+\delta$")
  ax.errorbar(x, y_TGC_M / y_SM, xerr=0.5, marker="o", ls='none', label="$-\delta$")
    
  ax.set_ylabel(r"$\sigma_0^{TGC} / \sigma_0^{SM}$", fontsize=30)
  
  process_str = PN.process_str(None,mu_charge)
  ax.legend(fontsize=20, loc='center right', title="${}$\n$\sigma_0^{{SM}}={}$\n$\delta={}$".format(process_str,np.round(y_SM,1),dev_scale), title_fontsize=20)
  
  plt.tight_layout()
  
  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    fig.savefig("{}/WWTGCChiralCheck_xs0.{}".format(format_dir, format))
  
  plt.close(fig)
  
def create_asymm_plot(asymm_dict, output_base, dev_scale, mu_charge, 
                      output_formats=["pdf","png"]):
  """ Check how the asymmetry changes when reweighting.
  """
  fig = plt.figure(tight_layout=True, figsize=(9,6.5))
  ax = plt.gca()
  ax.set_xlim(0, 3)
  
  x = [0.5,1.5,2.5]
  ax.set_xticks(x, minor=True)
  x_ticks = [ "$g_1^{Z}$", "$\kappa_{\gamma}$", "$\lambda_{\gamma}$" ]
  plt.xticks(x, x_ticks, size='large')
  
  y_SM = asymm_dict["SM"]
  
  y_TGC_P = np.array([asymm_dict[point] for point in ["g1z +", "ka +", "la +"]])
  y_TGC_M = np.array([asymm_dict[point] for point in ["g1z -", "ka -", "la -"]])
  
  ax.plot([0,3], [0.,0.])
  ax.errorbar(x, y_TGC_P - y_SM, xerr=0.5, marker="*", ls='none', label="$+\delta$")
  ax.errorbar(x, y_TGC_M - y_SM, xerr=0.5, marker="o", ls='none', label="$-\delta$")
    
  ax.set_ylabel(r"$A_{LR}^{TGC} - A_{LR}^{SM}$", fontsize=30)
  
  process_str = PN.process_str(None,mu_charge)
  ax.legend(fontsize=20, loc='center right', title="${}$\n$A_{{LR}}^{{SM}}={}$\n$\delta={}$".format(process_str,np.round(y_SM,5),dev_scale), title_fontsize=20)
  
  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    fig.savefig("{}/WWTGCChiralCheck_asymm.{}".format(format_dir, format))
  
  plt.close(fig)

def get_point_counts(rdf, tgc_config_reader, points_dict):
  """ Get the total weight sum the weights given in the weight dict.
  """
  counts_dict = {}
  for name, point in points_dict.items():
    index_str = str(tgc_config_reader.point_index(point) + 1)
    weight_name = "rescan_weights.weight" + index_str
    counts_dict[name] = rdf.Sum(weight_name)
  return counts_dict
  
def normalize_counts(counts_dict, norm):
  """ Normalise all counts in the dictionary by the given norm.
  """
  for name, count in counts_dict.items():
    counts_dict[name] = count.GetValue() * norm
  return counts_dict
    
def get_chiral_dicts(counts_LR, counts_RL):
  """ Transform the dictionaries with the chiral cross sections into two 
      dictionaries containing chiral asymmetry and total chiral cross section.
  """
  xs0_dict = {}
  asymm_dict = {}
  for key in counts_LR:
    xs_LR = counts_LR[key]
    xs_RL = counts_RL[key]
    xs0_dict[key] = xs_LR + xs_RL
    asymm_dict[key] = (xs_LR - xs_RL) / (xs_LR + xs_RL)
  return xs0_dict, asymm_dict
  
def check_chiral_dep(root_file_LR, root_file_RL, tgc_config_path, 
                     tgc_point_path, output_dir, mu_charge, 
                     tree_name = "WWObservables"):
  """ Check how reweighting changes the asymmetry and the total chiral cross 
      section for the given events in the ROOT file.
  """
  rdf_LR = PRH.skip_0weight(ROOT.RDataFrame(tree_name, root_file_LR))
  rdf_RL = PRH.skip_0weight(ROOT.RDataFrame(tree_name, root_file_RL))
  
  output_base = "{}/WWTGCChiralCheck/mu{}/".format(
                  output_dir, PN.sign_str(mu_charge,spelled=True))
  
  tcr = ITCR.TGCConfigReader(tgc_config_path, tgc_point_path)
  dev_scale = tcr.scale
  points_dict = {
    "g1z +": [1,0,0],
    "g1z -": [-1,0,0],
    "ka +":  [0,1,0],
    "ka -":  [0,-1,0],
    "la +":  [0,0,1],
    "la -":  [0,0,-1]
  }
  
  mu_rdf_LR = PRH.select_mu(rdf_LR, mu_charge)
  mu_rdf_RL = PRH.select_mu(rdf_RL, mu_charge)
  
  counts_LR = get_point_counts(mu_rdf_LR, tcr, points_dict)
  counts_RL = get_point_counts(mu_rdf_RL, tcr, points_dict)
  
  counts_LR["SM"] = mu_rdf_LR.Count() 
  counts_RL["SM"] = mu_rdf_RL.Count() 
  
  # Get the SM cross section, needed for correct uncertainty weighting
  SM_xsection_LR = mu_rdf_LR.Mean("cross_section")
  SM_xsection_RL = mu_rdf_RL.Mean("cross_section")
      
  # Calculate the normalisation of the MC events
  n_MC_LR = rdf_LR.Count() # Use the original rdf with all MC events (except 0-weight)
  n_MC_RL = rdf_RL.Count() # Use the original rdf with all MC events (except 0-weight)
  MC_norm_LR = SM_xsection_LR.GetValue() / n_MC_LR.GetValue()
  MC_norm_RL = SM_xsection_RL.GetValue() / n_MC_RL.GetValue()
      
  # Transform the counts into cross sections
  counts_LR = normalize_counts(counts_LR, MC_norm_LR)
  counts_RL = normalize_counts(counts_RL, MC_norm_RL)
  
  # From the chiral cross sections get the asymmetry and total chiral cross sec.
  xs0_dict, asymm_dict = get_chiral_dicts(counts_LR, counts_RL)
      
  create_xs0_plot(xs0_dict, output_base, dev_scale, mu_charge)
  create_asymm_plot(asymm_dict, output_base, dev_scale, mu_charge)
      
def main():
  """ Create check plots that investigate how the TGCs affect the chiral 
      asymmetry and the total chiral cross section.
  """
  ROOT.EnableImplicitMT() # Enable multithreading in RDataFrame
  ROOT.gROOT.SetBatch(True) # Don't show graphics at runtime
  PDF.set_default_mpl_format()
  
  config_dir = "/afs/desy.de/group/flc/pool/beyerjac/TGCAnalysis/SampleProduction/MCProduction/PrEWSampleProduction/scripts/config"
  tgc_config_path = config_dir + "/tgc.config"
  tgc_point_path = config_dir + "/tgc_dev_points_g1z_ka_la.config"
  
  output_dir = "../../output"
  
  RL_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eR_pL.root"
  LR_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eL_pR.root"
  
  check_chiral_dep(LR_path, RL_path, tgc_config_path, tgc_point_path, 
                   output_dir, None)
  check_chiral_dep(LR_path, RL_path, tgc_config_path, tgc_point_path, 
                   output_dir, +1)
  check_chiral_dep(LR_path, RL_path, tgc_config_path, tgc_point_path, 
                   output_dir, -1)

if __name__ == "__main__":
  main()