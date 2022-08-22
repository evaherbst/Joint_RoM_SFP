# Joint_RoM_SFP

:pencil:  Please read our papers:

[1](https://doi.org/10.1111/joa.13717): Herbst, E.C.\*, Eberhard, E.A.\*, Hutchinson, J.R. & Richards, C.T. (2022) Spherical frame projections for visualising joint range of motion, and a complementary method to capture mobility data. Journal of Anatomy. \* shared first authorship

[2](https://doi.org/10.1111/joa.13738): Herbst, E.C., Eberhard, E.A., Richards, C.T. & Hutchinson, J.R. (2022) In vivo and ex vivo range of motion in the fire salamander Salamandra salamandra. Journal of Anatomy.

If you use this method, please cite our papers (Herbst and Eberhard et al. 2002) [![DOI:https://doi.org/10.1111/joa.13717](http://img.shields.io/badge/DOI-10.1111/joa.13717-GREEN.svg)](https://doi.org/10.1111/joa.1371) and (Herbst et al. 2022) [![DOI:https://doi.org/10.1111/joa.13738](http://img.shields.io/badge/DOI-10.1111/joa.13738-GREEN.svg)](https://doi.org/10.1111/joa.13738) and the doi of the most recent Github release:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6914547.svg)](https://doi.org/10.5281/zenodo.6914547)


This repository contains MATLAB code to transform Qualisys motion capture data to joint centered movement, using a custom joint range of motion rig and the SFP method of visualization (Herbst et al. in review).

It contains example code to assign anatomical coordinate system to salamander bones (using the knee as a case study). 

Required inputs are motion capture data stored as .mat files, and a file of points placed on the perimeters of the articular surfaces of the bone, as well as on the plate that the bone was attached to in the rig (saved as a .txt file).

If joint-centric data is already available, we direct users directly to the SFP visualization repository (https://github.com/eeberhard/SFP).
