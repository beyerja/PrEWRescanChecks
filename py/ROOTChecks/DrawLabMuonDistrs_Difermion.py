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

def create_plot(x, y, h_range, m_label, Z_dir, chirality, 
                output_base, output_formats=["pdf","png"]):
  """ Create the plot for one specific mumu distribution.
  """
  fig = plt.figure(figsize=(7.5,6), tight_layout=True)
  ax = plt.gca()
  
  process_str = PN.difermion_process_str("mu", chirality, m_label, Z_dir)
  ax.set_title("${}$".format(process_str))
  
  z_max = 55 if m_label == "return_to_Z" else 320

  hist = ax.hist2d(
    x=x[:,0], y=x[:,1], weights=y,
    bins=(h_range[0],h_range[3]), 
    range=[[h_range[1], h_range[2]],[h_range[4], h_range[5]]],
    vmin=0, vmax=z_max )
  
  fig.colorbar(hist[3], ax=ax, label="$d\\sigma [$fb$]$")
  ax.set_xlabel("$\\cos \\theta_{\\mu^-}$")
  ax.set_ylabel("$\\cos \\theta_{\\mu^+}$")
  
  
  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    label = m_label if Z_dir=="" else m_label + "_" + Z_dir
    fig.savefig("{}/{}_{}.{}".format(
      format_dir, label, chirality, format) )

def book_TH2D(rdf, h_range, m_ranges, FB_split):
  """ Book all the necessary TH2D's, so that only one loop over the RDF is 
      required.
  """
  obs_names = ("costh_f", "costh_fbar")
  
  hist_dict = {}
  for m_label, m_range in m_ranges.items():
    hist_dict[m_label] = {}
    if FB_split[m_label]:
      for Z_direction in ["FZ", "BZ"]:
          mu_rdf = PRH.select_mumu(rdf, m_range[0], m_range[1], Z_direction)
          h_name = m_label + Z_direction
          h_def = (h_name, h_name, *h_range)
          hist_dict[m_label][Z_direction] = mu_rdf.Histo2D(h_def, *obs_names)
    else:
      mu_rdf = PRH.select_mumu(rdf, m_range[0], m_range[1], None)
      h_name = m_label
      h_def = (h_name, h_name, *h_range)
      hist_dict[m_label][""] = mu_rdf.Histo2D(h_def, *obs_names)
      
  return hist_dict

def check_labmu_distr(root_file, output_dir, 
                      tree_name="DifermionObservables"):
  rdf = ROOT.RDataFrame(tree_name, root_file)
  
  chirality = IFH.find_chirality(root_file)
  output_base = "{}/MuLabDistr_Difermion/".format(output_dir)
  
  obs_name = "costh_l"
  nbins, xmin, xmax = 30, -1., 1.
  h_range = (nbins, xmin, xmax, nbins, xmin, xmax)
  
  energy = 250
  m_ranges = {
    "return_to_Z": [81, 101],
    "high_Q2": [180, 1.1*energy] }
  FB_split = {
    "return_to_Z": True,
    "high_Q2": False }

  # Book all histograms
  hist_dict = book_TH2D(rdf, h_range, m_ranges, FB_split)
  
  # Get the SM cross section, needed for correct uncertainty weighting
  SM_xsection = rdf.Mean("cross_section")
      
  # Calculate the normalisation of the MC events
  n_MC = rdf.Count() 
  MC_norm = SM_xsection.GetValue() / n_MC.GetValue()
  
  # First transform hists to arrays and plot
  for m_label, subdict in hist_dict.items():
    for Z_dir, hist in subdict.items():
      x, y = PRHH.TH2_to_arrays(hist)
      y *= MC_norm
      create_plot(x, y, h_range, m_label, Z_dir, chirality, output_base)
      
def main():
  ROOT.EnableImplicitMT() # Enable multithreading in RDataFrame
  ROOT.gROOT.SetBatch(True) # Don't show graphics at runtime
  PDF.set_default_mpl_format()
  
  output_dir = "../../output"
  
  RL_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/2f_Z_l_eR_pL.root"
  LR_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/2f_Z_l_eL_pR.root"
  
  check_labmu_distr(RL_path, output_dir)
  check_labmu_distr(LR_path, output_dir)
  
if __name__ == "__main__":
  main()