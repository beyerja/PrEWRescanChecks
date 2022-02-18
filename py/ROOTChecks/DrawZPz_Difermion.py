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

def check_ZPz_distr(root_file, output_dir, tree_name="DifermionObservables", 
                       output_formats=["pdf","png"]):
  rdf = ROOT.RDataFrame(tree_name, root_file)
  mu_rdf = PRH.select_mumu(rdf, m_min=81, m_max=101)
  
  chirality = IFH.find_chirality(root_file)
  output_base = "{}/ZPz_Difermion/".format(output_dir)
  
  obs_name = "pz_ff"
  energy = 250.0
  nbins, xmin, xmax = 1000, -energy/2, energy/2
  h_def = (obs_name, obs_name, nbins, xmin, xmax)
  rhist = mu_rdf.Histo1D(h_def, obs_name)
  
  # Get the SM cross section, needed for correct uncertainty weighting
  SM_xsection = rdf.Mean("cross_section")
      
  # Calculate the normalisation of the MC events
  n_MC = rdf.Count() 
  MC_norm = SM_xsection.GetValue() / n_MC.GetValue()
  
  x,y = PRHH.TH1_to_arrays(rhist) 
  y *= MC_norm
      
  fig = plt.figure(figsize=(7.5,5.5), tight_layout=True)
  ax = plt.gca()
  
  colors =  plt.rcParams['axes.prop_cycle'].by_key()['color']
  
  h = ax.hist(x, weights=y, range=(xmin,xmax), bins=nbins,
              ls="-", lw=3, histtype=u'step', color=colors[0])
  print(h[1][np.argmax(h[0])])
  
  y_max = 150
  
  y_arrow = 0.8 * y_max
  arrow_kwargs = {"zorder": 2, "width": 2, "length_includes_head": True,
                  "head_width": 7, "head_length": 3, 
                  "fc": colors[1], "ec": "black"}
  ax.arrow(0, y_arrow, -20, 0, **arrow_kwargs)
  ax.arrow(0, y_arrow, 20, 0, **arrow_kwargs)
  
  y_line = y_max
  line_kwargs = {"color": colors[1], "ls": "--", "lw": 2}
  ax.plot([0,0], [0,y_max], **line_kwargs)
  
  text_kwargs = {"color": colors[1], "fontsize":25, "va": "center"}
  ax.text(-25, y_arrow, "$p_z^{\\mu\\mu}<0$", ha="right", **text_kwargs)
  ax.text(25, y_arrow, "$p_z^{\\mu\\mu}>0$", ha="left", **text_kwargs)
  
  process_str = PN.difermion_process_str("mu", chirality, label="return_to_Z")
  ax.legend(title="${}$".format(process_str),loc="upper left", 
            title_fontsize=20, frameon=False, borderpad=0)
  ax.set_xlabel("$p_z^{\\mu\\mu}\\,[$GeV$]$")
  ax.set_ylabel("$d\\sigma [$fb$^{{-1}}]$")
  ax.set_xlim((xmin,xmax))
  ax.set_ylim((0,y_max))
  
  peak_pz = np.round(h[1][np.argmax(h[0])])
  x = [-peak_pz,0,peak_pz]
  ax.set_xticks(x)
  ax.locator_params(axis='y', nbins=4)
  
  for format in output_formats:
    format_dir = "{}/{}".format(output_base,format)
    ISH.create_dir(format_dir)
    fig.savefig("{}/ZPz_{}.{}".format(format_dir, chirality, format) )
      
def main():
  ROOT.EnableImplicitMT() # Enable multithreading in RDataFrame
  ROOT.gROOT.SetBatch(True) # Don't show graphics at runtime
  PDF.set_default_mpl_format()
  
  output_dir = "../../output"
  
  RL_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction"+\
            "/NewMCProduction/2f_Z_l/2f_Z_l_eR_pL.root"
  LR_path = "/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction"+\
            "/NewMCProduction/2f_Z_l/2f_Z_l_eL_pR.root"
  
  check_ZPz_distr(RL_path, output_dir)
  check_ZPz_distr(LR_path, output_dir)
  
if __name__ == "__main__":
  main()