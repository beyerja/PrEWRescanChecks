import matplotlib.pyplot as plt
import numpy as np
import ROOT

# Use the the local helper modules
import sys
sys.path.append("../Helpers") 
import IO.SysHelpers as ISH
import IO.TGCConfigReader as ITCR
import NumPyHelp.Histograms as NHPH
import Plotting.DefaultFormat as PDF
import Plotting.MPLHelp as PMH
import Plotting.Naming as PN
import Plotting.RDFHelp as PRH
import Plotting.ROOTHistHelp as PRHH

def ratio(a,b,default=0.0):
  """ Calculate ratio between two arrays, if denominator point is 0 set default.
  """
  return np.where(np.abs(b)>0.01,a/b,default)

def norm_chisq(SM_array, dev_array, MC_norm):
  """ Calculate the chi-squared from the total normalisation information alone.
      (Assumes SM measured)
  """
  SM_count = np.sum(SM_array) * MC_norm
  dev_count = np.sum(dev_array) * MC_norm
  return (SM_count - dev_count)**2 / dev_count 
  
  
def total_chisq(SM_array, dev_array, MC_norm):
  """ Calculate the chi-squared from the all available information.
      (Assumes SM measured)
  """
  SM_normed = SM_array * MC_norm
  dev_normed = dev_array * MC_norm
  return np.sum( ratio((SM_normed - dev_normed)**2,dev_normed) )
  
def shape_chisq(chisq_total, chisq_norm):
  """ Calculate the chi-squared that comes from the shape information, assuming
      that the total uncertainty comes from combining independent chi-squared's
      of shape and normalisation.
  """
  return chisq_total - chisq_norm
  
def third_coord_chisq(x, coord_index, SM_array, dev_array, MC_norm):
  """ Calculate indirectly the importance of one coordinate for the chi-squared:
      - Calculate the total chi-sqaured from the 3D distribution
      - Calculate the chi-squared from the 2D distribution that integrates over
        the relevant coordinate
      => The "chi-squared" of the coordinate is the difference between the two.
  """
  # Project away the relevant coordinate
  _, SM_proj  = NHPH.proj_2from3D(x, SM_array, coord_index)
  _, dev_proj = NHPH.proj_2from3D(x, dev_array, coord_index)
  
  chisq_total = total_chisq(SM_array, dev_array, MC_norm) 
  chisq_proj  = total_chisq(SM_proj, dev_proj, MC_norm) 
  return chisq_total - chisq_proj
  
def draw_grid(ax,x,ls,color="black",alpha=1.0,lw=0.5):
  """ Helper function to simple axis grid lines at the given points
  """
  # Get the current limits to ensure they won't be changed
  xlims = ax.get_xlim()
  ylims = ax.get_ylim()

  # Draw all requested grid lines
  for _x in x:
    ax.plot([x,x],ylims,ls=ls,color=color,alpha=alpha,lw=lw)
    
  # Reset the limits
  ax.set_xlim(xlims)
  ax.set_ylim(ylims)
  
def draw_x_labelling(ax, dy):
  """ Draw the x axis labelling.
  """
  ax.set_xticklabels([])

  dev_labels = [ "$\Delta g_1^{Z}$", "$\Delta \kappa_{\gamma}$", "$\Delta \lambda_{\gamma}$"]
  for i_TGC in range(3):
    ax.text(i_TGC*2+1, -3*dy, dev_labels[i_TGC], ha='center', va='bottom')
    
  for xi_M in np.arange(0, 6, step=2)+0.5:
    ax.text(xi_M, -2*dy, "$-\delta$", ha='center', va='bottom')
  for xi_P in np.arange(1, 6, step=2)+0.5:
    ax.text(xi_P, -2*dy, "$+\delta$", ha='center', va='bottom')
  
  for xi_L in np.arange(6)+0.25:
    ax.text(xi_L, -dy, "LR", ha='center', va='bottom')
  for xi_R in np.arange(6)+0.75:
    ax.text(xi_R, -dy, "RL", ha='center', va='bottom')
  
  
def create_sig_plot(LR_chisqs, RL_chisqs, output_base, dev_scale, process_str, 
                    lumi, output_formats=["pdf","png"]):
  """ Create a plot that showcases where the TGC sensitivities come from.
  """
  fig = plt.figure(tight_layout=True, figsize=(10,8.5))
  ax = plt.gca()
  
  x = np.arange(6)+0.5
  x_L = x-0.25 
  x_R = x+0.25
  ax.set_xlim(0, x_R[-1]+0.25)
  ax.set_xticks(x, minor=True)
  order = ["g1z -", "g1z +", "ka -", "ka +", "la -", "la +"]
         
  # Collect the y points
  chisq_L_total = np.array([LR_chisqs[point]["total"] for point in order])
  chisq_L_shape = np.array([LR_chisqs[point]["shape"] for point in order])
  chisq_L_norm = np.array([LR_chisqs[point]["norm"] for point in order])
  
  chisq_R_total = np.array([RL_chisqs[point]["total"] for point in order])
  chisq_R_shape = np.array([RL_chisqs[point]["shape"] for point in order])
  chisq_R_norm = np.array([RL_chisqs[point]["norm"] for point in order])

  # Common style elements
  ms = 9.0 # Marker size
  ls = 'none' # Line style -> No connecting line
  m_L = "v" # eLpR marker
  m_R = "^" # eRpL marker
  width = 0.3 # x errorbar length 
  
  # Create the plots
  ax.bar(x_L, chisq_L_shape, width=width, label="shape")
  ax.bar(x_L, chisq_L_norm, width=width, bottom=chisq_L_shape, label="normalisation")
  
  plt.gca().set_prop_cycle(None) # Reset the color cycle
  
  ax.bar(x_R, chisq_R_shape, width=width)
  ax.bar(x_R, chisq_R_norm, width=width, bottom=chisq_R_shape)
  
  # Draw some axis grids
  draw_grid(ax,[1,3,5],"dotted")
  draw_grid(ax,[2,4],"solid")
  
  ymax = ax.get_ylim()[1]
  ax.set_ylim(0, ymax)
  ax.set_ylabel("$\\chi^2$", fontsize=30)
  
  dy = ymax/15.
  draw_x_labelling(ax, dy)
  
  ax.legend(ncol=2, fontsize=20, title="${}$, $L={}\,$ab$^{{-1}}$, $\delta={}$".format(process_str,np.round(lumi/1000,1),dev_scale), 
            title_fontsize=20, bbox_to_anchor=(0.0, 1.03), loc='lower left')
  
  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    fig.savefig("{}/WWTGCSensitivityCheck.{}".format(format_dir, format))
  
  plt.close(fig)

def create_shape_contr_plot(LR_chisqs, RL_chisqs, output_base, process_str, 
                            output_formats=["pdf","png"]):
  """ Create a plot that shows the estimated relative contribution of each 
      coordinate to the shape deviation significance. 
  """
  fig = plt.figure(tight_layout=True, figsize=(10,8.5))
  ax = plt.gca()
  
  x = np.arange(6)+0.5
  x_L = x-0.25 
  x_R = x+0.25
  ax.set_xlim(0, x_R[-1]+0.25)
  ax.set_xticks(x, minor=True)
  order = ["g1z -", "g1z +", "ka -", "ka +", "la -", "la +"]
         
  # Collect the y points
  chisq_L_c1 = np.array([LR_chisqs[point]["costh_Wminus_star"] for point in order])
  chisq_L_c2 = np.array([LR_chisqs[point]["costh_l_star"] for point in order])
  chisq_L_c3 = np.array([LR_chisqs[point]["phi_l_star"] for point in order])
  chisq_L_norm = chisq_L_c1 + chisq_L_c2 + chisq_L_c3
  for chisq_L_arr in [chisq_L_c1, chisq_L_c2, chisq_L_c3]:
    chisq_L_arr *= 100./chisq_L_norm
  
  chisq_R_c1 = np.array([RL_chisqs[point]["costh_Wminus_star"] for point in order])
  chisq_R_c2 = np.array([RL_chisqs[point]["costh_l_star"] for point in order])
  chisq_R_c3 = np.array([RL_chisqs[point]["phi_l_star"] for point in order])
  chisq_R_norm = chisq_R_c1 + chisq_R_c2 + chisq_R_c3
  for chisq_R_arr in [chisq_R_c1, chisq_R_c2, chisq_R_c3]:
    chisq_R_arr *= 100./chisq_R_norm
    
  # Labelling of the three coordinates
  c1_label = "${}$".format(PN.observable_str("costh_Wminus_star", "WW"))
  c2_label = "${}$".format(PN.observable_str("costh_l_star", "WW"))
  c3_label = "${}$".format(PN.observable_str("phi_l_star", "WW"))

  # Common style elements
  ms = 9.0 # Marker size
  ls = 'none' # Line style -> No connecting line
  m_L = "v" # eLpR marker
  m_R = "^" # eRpL marker
  width = 0.25 # x errorbar length 
  
  # Create the plots
  PMH.skip_n_colors(ax, 2)
  ax.bar(x_L, chisq_L_c1, width=width, label=c1_label)
  ax.bar(x_L, chisq_L_c2, width=width, bottom=chisq_L_c1, label=c2_label)
  ax.bar(x_L, chisq_L_c3, width=width, bottom=chisq_L_c1+chisq_L_c2, label=c3_label)
  
  plt.gca().set_prop_cycle(None) # Reset the color cycle
  
  PMH.skip_n_colors(ax, 2)
  ax.bar(x_R, chisq_R_c1, width=width)
  ax.bar(x_R, chisq_R_c2, width=width, bottom=chisq_R_c1)
  ax.bar(x_R, chisq_R_c3, width=width, bottom=chisq_R_c1+chisq_R_c2)
  
  dy = ax.get_ylim()[1]/15.
  draw_x_labelling(ax, dy)
  
  # Draw some axis grids
  draw_grid(ax,[1,3,5],"dotted")
  draw_grid(ax,[2,4],"solid")
  
  ax.set_ylim(0, 100)
  ax.set_ylabel("$\\frac{\chi^2_{w/\,\, coord.} - \chi^2_{w/o\,\, coord.}}{norm.}$ [%]", fontsize=30)
  
  ax.legend(ncol=3, fontsize=20, title_fontsize=20, bbox_to_anchor=(0.0, 1.03), 
            loc='lower left', title="Rel. importance of coord. in 3D distr., ${}$".format(process_str))
  
  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    fig.savefig("{}/WWTGCShapeContributions.{}".format(format_dir, format))
  
  plt.close(fig)

def get_chisqs(root_file, tgc_config_path, tgc_point_path,
               observables, mu_charge, lumi, tree_name = "WWObservables"):
  """ Calculate the chi-squared for the TGC deviations given on the distribution
      in the given ROOT file, assuming the given luminosity.
  """
  # Get the RDataFrame (skipping failed rescan events)
  rdf = PRH.skip_0weight(ROOT.RDataFrame(tree_name, root_file))
  mu_rdf = PRH.select_mu(rdf, mu_charge) # Events with given mu charge
  
  # Find the indices to the TGC points are being tested
  tcr = ITCR.TGCConfigReader(tgc_config_path, tgc_point_path)
  index_dict = {
    "g1z +": tcr.point_index([1,0,0]),
    "g1z -": tcr.point_index([-1,0,0]),
    "ka +": tcr.point_index([0,1,0]),
    "ka -": tcr.point_index([0,-1,0]),
    "la +": tcr.point_index([0,0,1]),
    "la -": tcr.point_index([0,0,-1])
  }
  
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
  hist_SM = mu_rdf.Histo3D(h_SM_def, obs_names[0], obs_names[1], obs_names[2])
  
  # Histograms with TGC deviation
  for point, index in index_dict.items():
    i_str = str(index+1) # weight naming starts with index 1
    h_def = ("w"+i_str, "w"+i_str, obs_1[0], obs_1[1], obs_1[2],
                                   obs_2[0], obs_2[1], obs_2[2],
                                   obs_3[0], obs_3[1], obs_3[2])
    weight_name = "rescan_weights.weight"+i_str
    hist_dict[point] = mu_rdf.Histo3D(h_def, 
                                      obs_names[0], obs_names[1], obs_names[2], 
                                      weight_name)
    
  # Get the SM cross section, needed for correct uncertainty weighting
  SM_xsection = mu_rdf.Mean("cross_section")
  
  # Calculate the normalisation of the MC events
  n_MC = rdf.Count() # Use the original rdf with all MC events (except 0-weight)
  MC_norm = SM_xsection.GetValue() * lumi / n_MC.GetValue()
  
  # Find the chi-squared's from normalisation, shape, and both
  chisq_dict = {}
    
  # Calculate the chi-squared's
  x, y_SM = PRHH.TH3_to_arrays(hist_SM)
  for point, hist in hist_dict.items():
    y_dev = PRHH.TH3_to_arrays(hist)[1]
  
    chisq_total = total_chisq(y_SM, y_dev, MC_norm)
    chisq_norm = norm_chisq(y_SM, y_dev, MC_norm)
    chisq_shape = shape_chisq(chisq_total, chisq_norm)
    chisq_dict[point] = {
      "total": chisq_total, "norm": chisq_norm, "shape": chisq_shape }
      
    for i_c in range(len(obs_names)):
      chisq_dict[point][obs_names[i_c]] = third_coord_chisq(x, i_c, y_SM, y_dev, 
                                                            MC_norm)
    
  return chisq_dict
  

def check_sensitivities(LR_file, RL_file, tgc_config_path, tgc_point_path, 
                        output_dir, observables, mu_charge, lumi):
  """ Find the chi-squared's caused by different available informations and 
      trigger the plotting.
  """
  
  # Find some meta info needed for proper output/labelling
  output_base = "{}/TGCSensitivityCheck/mu{}/".format(
                  output_dir, PN.sign_str(mu_charge,spelled=True))
  process_str = PN.process_str(None, mu_charge)
  dev_scale = ITCR.TGCConfigReader(tgc_config_path, tgc_point_path).scale

  # Calculate the chi-squared's 
  chisqs_LR = get_chisqs(LR_file, tgc_config_path, tgc_point_path, observables, mu_charge, lumi)
  chisqs_RL = get_chisqs(RL_file, tgc_config_path, tgc_point_path, observables, mu_charge, lumi)
  
  # Trigger the plotting
  create_sig_plot(chisqs_LR, chisqs_RL, output_base, dev_scale, process_str, lumi)
  create_shape_contr_plot(chisqs_LR, chisqs_RL, output_base, process_str)
      
def main():
  """ Check where the sensitivity to the different TGCs comes from:
      from shape or from normalisation
      (Assumption: Those are the two available kinds of info for 
                   a given chiral initial state)
  """
  ROOT.EnableImplicitMT() # Enable multithreading in RDataFrame
  ROOT.gROOT.SetBatch(True) # Don't show graphics at runtime
  PDF.set_default_mpl_format()
  
  config_dir = "/afs/desy.de/group/flc/pool/beyerjac/TGCAnalysis/SampleProduction/MCProduction/PrEWSampleProduction/scripts/config"
  tgc_config_path = config_dir + "/tgc.config"
  tgc_point_path = config_dir + "/tgc_dev_points_g1z_ka_la.config"
  
  output_dir = "../../output"
  
  observables = { "costh_Wminus_star" : (20, -1., 1.),
                  "costh_l_star" :      (10, -1., 1.),
                  "phi_l_star":         (10, -np.pi, np.pi) }
  
  lumi = 1000.0 # 2ab^-1 total ~> 1ab^-1 for each of the two allowed chiral states
  
  RL_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eR_pL.root"
  LR_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eL_pR.root"
  
  check_sensitivities(LR_path, RL_path, tgc_config_path, tgc_point_path, output_dir, observables, +1, lumi)
  check_sensitivities(LR_path, RL_path, tgc_config_path, tgc_point_path, output_dir, observables, -1, lumi)
  check_sensitivities(LR_path, RL_path, tgc_config_path, tgc_point_path, output_dir, observables, None, lumi)
  

if __name__ == "__main__":
  main()