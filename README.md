# OpenHarmonics

OpenHarmonics is a tool for performing non-linear harmonic analysis studies of distribution systems. A key feature of the tool is that it captures the non-linearity of grid connected power electronics. The tool was developed as a part of my Masters research in the Energy Systems Group at the University of Toronto in 2013. OpenHarmonics has the ability to analyze the voltage and current harmonic spectra in locations of interest in the simulated power system. OpenHarmonics is built on-top of the OpenDSS open source distribution system simulator. Its primary purpose is to extend the harmonic analysis functionality of OpenDSS to include non-linear power electronic converter models. As such, the advantages of OpenHarmonics are realized when the distribution system to be analyzed has loads and generators interfaced to the grid through power electronic converters. The derived converter models are more computationally intensive than the converter models currently being used in the commercial harmonic analysis software programs. However, the derived converter models are orders of magnitude more accurate, allowing accurate prediction of non-characteristic harmonics (e.g. positive sequence 3rd harmonic) that are entirely neglected in linear commercial harmonic analysis software. OpenHarmonics was developed in MATLAB and the user interacts with OpenHarmonics through the MATLAB workspace.

For installation instructions and an overview of the tool please refer to the included User Manual. 

For a more detailed description of the tool, please refer to [1] and [2].

If OpenHarmonics is used for research purposes, please cite [1].

[1] P. Gray, “Harmonic Models of Common Converter Topologies for Accurate Harmonic Analysis of Distribution Systems,” Thesis, 2013.
[2] P.A. Gray, P.W. Lehn, and G. Kish, "Harmonic Analysis Software Tool for Steady-State Analysis of Distorted Grids with Multiple Distributed Generators," in Proc. Canada Conf.
CIGRE, Sep. 2014, pp. 1–8.  

In case of any questions or issues related to the tool please post on the "Issues" tab and I will try to answer.
