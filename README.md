
# <center>  mirror-state-processing</center>

This repository is served as the platform for sharing the codes of miror state processing in Polarization-sensitive optical coherence tomography.

**Mirror state constraint reveals that the polarization evovlement in an reciprocal and unitary system always tends to pass through a specific polarizatio state and this constraint simpliies the local retardation reconstruction.** This project demonstrates the process of local retardation reconstruction with single polarization state. Also, we compared the results with that of benchmark dual input processing.

#V1.0
Drafted by Qiaozhou Xiong, Nanyang Technological University. 

To run this version, you need corresponding cardiovascular dataset. Sample dataset will be uploaded in latter version. Tool box developed by Brett E. Bouma group is required to run the codes.

## Folder and file description

*Intravascular 2019-07-09dist*: results with reliability metric associated with distance and SNR

*Intravascular 2019-07-09sina*: results with reliability metric associated with sina and SNR

*Reliability metric generation*: Genearte the reliability metrics by Montel Carlo experiment

*single input processing toolbox*: toolbox needed for dual input processing, single input processing with mirror state constraint, comparison visualization

*FullSpectral Mirror state.mat,Binned Mirror state.mat*: mirror state corresponds to full spectrum and binned specturm

*FilesNameForComparison.mat, FilesR1.mat*: files name and setting for comparsion from the whole data set.

*Step1_mirrorStateEstimation.m*: The main scripts for the mirror state estimation, the output will be two mat file contains the full spectral mirror state and binned mirror state. The main fucntion would call the mirrorStateExtract function in the single input processing toolbox.

*Step2_mirrorpointConstructionAnalysis.m* The main scripts calls the function of  single input processing and dual input processing and results comparison visualization from the single input processing tool box.
