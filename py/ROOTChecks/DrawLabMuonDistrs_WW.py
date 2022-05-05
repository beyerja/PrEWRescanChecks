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

def get_hist(root_file, mu_charge, tree_name="WWObservables"):
  rdf = PRH.skip_0weight(ROOT.RDataFrame(tree_name, root_file))
  mu_rdf = PRH.select_mu(rdf, mu_charge)
  
  obs_name = "costh_l"
  nbins, xmin, xmax = 30, -1., 1.
  
  # Book the histogram
  h_def = (obs_name, obs_name, nbins, xmin, xmax)
  rhist = mu_rdf.Histo1D(h_def, obs_name)
  
  # Get the SM cross section
  SM_xsection = rdf.Mean("cross_section")
      
  # Calculate the normalisation of the MC events
  n_MC = rdf.Count()
  MC_norm = SM_xsection.GetValue() / n_MC.GetValue()
  
  x,y = PRHH.TH1_to_arrays(rhist) 
  y *= MC_norm
  return x, y

def check_labmu_distr(root_file_LR, root_file_RL, output_dir, mu_charge, 
                      output_formats=["pdf","png"]):
  output_base = "{}/MuLabDistr_WW/".format(output_dir)
  process_str = PN.process_str(None, mu_charge)
  
  obs_name = "costh_l"
  nbins, xmin, xmax = 30, -1., 1.
  
  x_LR, y_LR = get_hist(root_file_LR, mu_charge)
  x_RL, y_RL = get_hist(root_file_RL, mu_charge)
      
  fig = plt.figure(figsize=(7.5,6), tight_layout=True)
  ax = plt.gca()
  
  ax.hist(x_LR, weights=y_LR, 
          ls="-", lw=3, range=(xmin,xmax), bins=nbins, histtype=u'step',
          label="$e_{{L}}^{{-}}e_{{R}}^{{+}}$")
  ax.hist(x_RL, weights=y_RL, 
          ls="--", lw=3, range=(xmin,xmax), bins=nbins, histtype=u'step',
          label="$e_{{R}}^{{-}}e_{{L}}^{{+}}$")
  
  ax.legend(title="${}$".format(process_str), loc="upper center")
  ax.set_xlabel("$\\cos \\theta_{\\mu}$")
  ax.set_ylabel("$d\\sigma [$fb$]$")
  ax.set_xlim((xmin,xmax))
  ax.set_ylim((0,ax.get_ylim()[1]))
  
  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    fig.savefig("{}/{}.{}".format(
      format_dir, PN.sign_str(mu_charge,spelled=True), format) )
      
def main():
  ROOT.EnableImplicitMT() # Enable multithreading in RDataFrame
  ROOT.gROOT.SetBatch(True) # Don't show graphics at runtime
  PDF.set_default_mpl_format()
  
  output_dir = "../../output"
  
  RL_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eR_pL.root"
  LR_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/4f_WW_sl_eL_pR.root"
  
  check_labmu_distr(LR_path, RL_path, output_dir, +1)
  check_labmu_distr(LR_path, RL_path, output_dir, -1)

if __name__ == "__main__":
  main()