function [measure] = sampleMeasure_Above(nSources,upperindex)
%sampleMeasure_Above - sampling a new measure from "top-down"
% The order of sampling a brand new measure is to first sample (nSource-1)-tuples (e.g. g_1234 for 5-source),
% then (nSource-2)-tuples (e.g. g_123), and so on until singletons (g_1)
% Notice it will satisfy monotonicity!
%
% INPUT
%   nSources - number of sources
%   upperindex - the cell that stores all the corresponding supersets (upper index) of measure elements
%
% OUTPUT
%   - measure - new measure after update
%
% Written by: X. Du 03/2018
%

%sample new measure, prior is uniform/no prior
Nmeasure = 2^nSources-1 ; %total length of measure
measure = zeros(1,Nmeasure);
measure(Nmeasure) = 1; %g_all
measure(Nmeasure-1:-1:Nmeasure-nSources) = rand(1,nSources); %g_1234,1235,1245,1345,2345 for 5 sources, for example (second to last tier)

for j = (Nmeasure-nSources-1):-1:1
    upperBound = min(measure(upperindex{j})); %upper bound
    measure(j)= 0 + (upperBound - 0)*rand(1);
end
        
end

