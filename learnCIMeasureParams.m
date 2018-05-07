function [Parameters] = learnCIMeasureParams()
% Parameters for learnCIMeasure_noisyor()
% OUTPUT 
%  Parameters - struct - The struct contains the following fields:
%                % These parameters are user defined (can be modified)
%                    1. nPop: Size of population
%                    2. sigma: Sigma of Gaussians in fitness function
%                    3. nIterations: Number of iterations
%                    4. eta:Percentage of time to make small-scale mutation
%                    5. sampleVar: Variance around sample mean
%                    6. mean: mean of ci in fitness function. Is always set to 1 if the positive label is "1".
%                    7. analysis: if ="1", record all intermediate results
%
               
Parameters.nPop = 30; %Size of population
Parameters.sigma = 0.1; %Sigma of Gaussians in fitness function
Parameters.maxIterations = 5000; %max number of iteration
Parameters.fitnessThresh = 0.0001;%stopping criteria: when fitness change is less than the threshold
Parameters.eta = 0.8; %Percentage of time to make small-scale mutation
Parameters.sampleVar = 0.1; %variance around sample mean
Parameters.mean = 1; % mean of ci in fitness function. Is always set to 1 if the positive label is "1".
Parameters.analysis = 1; % if ="1", record all intermediate results
end

 