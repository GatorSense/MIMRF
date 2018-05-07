function [measure] = sampleMeasure_Bottom(nSources,lowerindex)
%sampleMeasure_Bottom - sampling a new measure from "bottom-up"
% The order of sampling a brand new measure is to first sample singleton (e.g. g_1..),
% then duplets (e.g. g_12..), then triples (e.g. g_123..), etc..
% Notice it will satisfy monotonicity!
%
% INPUT
%   nSources - number of sources
%   lowerindex - the cell that stores all the corresponding subsets (lower index) of measure elements
%
% OUTPUT
%   - measure - new measure after update
%
% Written by: X. Du 03/2018
%

%sample new measure, prior is uniform/no prior
measure = zeros(1,(2^nSources-1));
measure(1:nSources) = rand(1,nSources); %sample densities
measure(end) = 1; %g_all

for j = (nSources+1) : (size(measure,2)-1)
    lowerBound = max(measure(lowerindex{j})); %lower bound     
    measure(j)= lowerBound + (1-lowerBound)*rand(1);
end

end