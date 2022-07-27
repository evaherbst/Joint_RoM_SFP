# Joint_RoM_SFP

This repository contains MATLAB code to transform Qualisys motion capture data to joint centered movement, using a custom joint range of motion rig and the SFP method of visualization (Herbst et al. in review).

It contains example code to assign anatomical coordinate system to salamander bones (using the knee as a case study). 

Required inputs are motion capture data stored as .mat files, and a file of points placed on the perimeters of the articular surfaces of the bone, as well as on the plate that the bone was attached to in the rig (saved as a .txt file).

If joint-centric data is already available, we direct users directly to the SFP visualization repository (https://github.com/eeberhard/SFP).
