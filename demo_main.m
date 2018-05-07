%% Main demo script for MIMRF multi-resolution multi-modal fusion algorithm
%
% Note: This demo shows results on a set of simulated multi-resolution 
% images. However, this algorithm is a general framework that can also 
% work with multi-modal data sets, such as rasterized hyperspectral 
% imagery versus LiDAR point clouds (as shown in the paper).
%
% Note: This code uses uses MATLAB Statistics and Machine Learning Toolbox, 
% MATLAB Optimization Toolbox and MATLAB Parallel Computing Toolbox. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by: X. Du and A. Zare
% Latest revision: May 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This product is Copyright (c) 2018 X. Du and A. Zare
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%   3. Neither the name of the University nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF MISSOURI AND
% CONTRIBUTORS ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED.  IN NO EVENT SHALL THE UNIVERSITY OR CONTRIBUTORS
% BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES,
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%

%% Demo: load a demo 3-source data

disp('This demo takes ~5 min.')
prompt = 'Running this code will clear your workspace.Would you like to continue?[Y/n]';
str = input(prompt,'s');
if (str=='N') || (str=='n')
    ;
elseif (str=='Y') || (str=='y')
    
clear;close all;clc
addpath('util')

%%%%%%% If no mex file existed in the ./util folder, run the following two lines to compile mex file %%%%%%%
% mex computeci.c
% mex ismember_findrow_mex.c
%% The MU Fusion demo

%%%%%%% Generate a simulated multiresolution "MU" data set.
%%%%%%% The goal is to detect/highlight the letters "M" and "U" in the scene.
[Bags, Labels, Seg] = generateSimData_MU();
%%%%%%%  Training Stage: Learn measures given training Bags and Labels
[Parameters] = learnCIMeasureParams(); %user-set MIMRF parameters
[measure_MIMRF, initialMeasure_MIMRF, Analysis_MIMRF] = learnCIMeasure_minmax_multires(Bags, Labels, Parameters);%noisy-or model
%%%%%%%  Testing Stage: Given the learned measures above, compute and plot fusion results
[TestConfMap] = computeTestMap(Bags, Labels, measure_MIMRF, Seg);

load('demo_MultiRes_data_MU.mat')
figure(101);
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,2,1);imagesc(Img);title('True Labels')
subplot(1,2,2);imagesc(TestConfMap);colorbar;title('MIMRF Fusion result')


end
