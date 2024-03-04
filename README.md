# VIV_sFOM_POD_inference
Data and codes for non-intrusive reduced-order modelling for VIV problems, via sparse FOM inference.
The code is written in MATLAB 2017a.

%---------------------------------------------%

The corresponding theory and analysis is given in "Adjacency-based, non-intrusive model reduction for Vortex-Induced Vibrations", 
by L. Gkimisis, T. Richter, P. Benner.
https://arxiv.org/abs/2302.09981

%---------------------------------------------%

The corresponding flowfield datasets for Re=90, Re=180 are uploaded in Zenodo, in zipped format and can be downloaded from https://zenodo.org/records/10779375.
Oscillation data are collected in "functionals.txt", flowfield data are given in .csv format for each timestep in "Flowfield data".

%---------------------------------------------%

VIV_sFOM_POD_inference.m : Main code to be executed. The user is prompted to provide with the directory paths to the dataset and the MESH2D tool.
MESH2D can be downloaded from https://matlab.mathworks.com/open/fileexchange/v1?id=25555.
The dataset can be downloaded and unzipped from https://zenodo.org/records/10779375.

%---------------------------------------------%

Functions called in Full_FSI_3A_proj.m :

FE_solver_v3b.m : Construct a mesh by Delaunay triangulation using MESH2D and solve a linear Laplace problem for fluid mesh displacement.
natsortfiles : Read snapshot data files sorted over time.
custom_cmap : Custom colormap used for figures.
MESH2D Package : Unstructured mesh generation. See https://www.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-delaunay-based-unstructured-mesh-generation

%---------------------------------------------%
