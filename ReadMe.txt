Multi-Resolution Multi-Modal Sensor Fusion For Remote Sensing Data With Label Uncertainty       ----- ReadMe File

This folder Includes the paper and demo MATLAB code for The MIMRF fusion algorithm. 

***************************************************************
***NOTE: If the MIMRF algorithm is used in any publication or presentation, the following reference must be cited:

 X. Du and A. Zare, "Multi-Resolution Multi-Modal Sensor Fusion For Remote Sensing Data With Label Uncertainty," Under Review. Available: https://arxiv.org/abs/1805.00930

***************************************************************
***
The MIMRF Fusion Algorithm runs using the following function:

[Parameters] = learnCIMeasureParams();
[measure,initialMeasure, Analysis] = learnCIMeasure_minmax_multires(Bags, Labels, Parameters,trueInitMeasure)

The TrainBags input is a 1xNumTrainBags cell. Inside each cell, NumPntsInBag x nSources cell. Inside each cell, the "collection" of all possible combinations generated from the multi-resolution data set. Details please see Section 3 of the MIMRF paper.
The TrainLabels input is a 1xNumTrainBags double vector that only take values of "1" and "0" (two-class classfication problems). Training labels for each bag.

The parameters input is a struct with the following fields: 
  Parameters - struct - The struct contains the following fields:
                 These parameters are user defined (can be modified)
                    1. nPop: Size of population
                    2. sigma: Sigma of Gaussians in fitness function
                    3. nIterations: Number of iterations
                    4. eta:Percentage of time to make small-scale mutation
                    5. sampleVar: Variance around sample mean
                    6. mean: mean of ci in fitness function. Is always set to 1 if the positive label is "1".
                    7. analysis: if ="1", record all intermediate results
                 *Parameters can be modified in [Parameters] = learnCIMeasureParams() function.	
********************************************************************* 
****The root directory contains the following files:

demo_main.m                     - Run this. This is the main demo file for the MIMRF fusion algorithm. 
learnCIMeasureParams.m              - User-set Parameters
learnCIMeasure_minmax_multires   - MIMRF Fusion algorithm function with min-max objective function (Two-class classification).
evalFitness_minmax_multires.m            - Fitness function for the MIMRF fusion algorithm with min-max objective function.
ComputeTestMap.m      - Compute test map (MIMRF fusion results) after learning the fuzzy measures from training 

generateSimData_MU.m   - This function generates a simulated multi-resolution data set for demo purposes, based on "Img_MU.mat" as true labels.
demo_MultiRes_data_MU.mat 	          - A pre-defined multi-resolution image data set.
*The motivation behind this demo data set is that, suppose we have five sensors that produces different, multi-resolution detection results (Detector 1 and 2 can detect the letter shape "M", Detector 3 and 4 can detector the letter "U", and Detector 5 can detect the background). We use the MIMRF algorithm to fuse those sensor outputs and ideally obtain high confidence on both "M" and "U". All five detectors produces multi-resolution imagery. The "Seg1-Seg5" variables contains segmentation information used to construct bags.
** This demo will produce two images: figure 100 shows the five fusion sources with the bag construction, figure 101 shows the true labels and MIMRF fusion results.
*** Note that, the measure learned by the MIMRF algorithm will ideally have high values on g_{12} and g_{34} (the 6th and 13th elements in this example), indicating the combination of fusion source 1&2 and 3&4 are useful in the fusion task.

****The subfolders contains the following files:
MIMRF_paper.pdf   - The paper
./util            - All functions, as follows (Note: Many of the functions are used in the MICI algorithms also proposed by the authors. See https://github.com/GatorSense/MICI):
ChoquetIntegral_g_MultiSources.m   - Compute the Choquet integral output given a single data point "hx" and measure "g"
computeci.c                        - Compute the Choquet integral output given a data point "hx" and measure "g" in c code. *Need to run "mex computeci.c"
evalInterval      	-evaluate valid interval of a fuzzy measure
ismember_findrow_mex.c	- Find row index if vector A is part of a row in Vector B in c code. **Need to run in MATLAB "mex ismember_findrow_mex.c"
ismember_findrow_mex_my.m     - Output row index if vector A is part of a row in Vector B (uses above c code).
share.h    - global variable header to be used in computeci.c
invcdf_TruncatedGaussian.m       - compute inverse cdf for Truncated Gaussian.
sampleMeasure.m                  - either flip a coin and randomly sample a brand new measure from uniform between sampleMeasure_Above and sampleMeasure_Bottom; OR only sample and update one element in the measure.
sampleMeasure_Above.m            - sampling a new measure from "top-down".
sampleMeasure_Bottom.m           - sampling a new measure from "bottom-up".
sampleMultinomial_mat.m      - sample from a multinomial distribution.
rowcol.m 		-returns the row and column index given total pixel index number in an image


********************************************************************* 
***********************************************************************
Authors: Xiaoxiao Du, Alina Zare
Department of Electrical and Computer Engineering, University of Missouri
Department of Electrical and Computer Engineering, University of Florida
 Email Address: xdy74@mail.missouri.edu; azare@ece.ufl.edu
 Latest Revision: May 2018

This code uses MATLAB Statistics and Machine Learning Toolbox, 
MATLAB Optimization Toolbox and MATLAB Parallel Computing Toolbox. 
*Note: if you do not have MATLAB Parallel Toolbox, simply change "parfor" in the code to "for". 

% This product is Copyright (c) 2018 X. Du and A. Zare
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
% 1. Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 3. Neither the name of the University nor the names of its contributors
% may be used to endorse or promote products derived from this software
% without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF MISSOURI AND
% CONTRIBUTORS ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OR CONTRIBUTORS
% BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES,
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE
