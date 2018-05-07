
function [subsetInterval] = evalInterval(measure,nSources,lowerindex,upperindex)
% Evaluate the valid interval width of a measure, then sort in descending order.
%
% INPUT
%   measure -  measure to be evaluated after update
%   nSources - number of sources
%   lowerindex - the cell that stores all the corresponding subsets (lower index) of measure elements
%   upperindex - the cell that stores all the corresponding supersets (upper index) of measure elements
% OUTPUT
%   subsetInterval - the valid interval widths for the measure elements,before sorting by value
%
% Written by: X. Du 03/2018
%

Nmeasure = 2^nSources-1 ; %total length of measure
lb = zeros(1,Nmeasure-1);
ub = zeros(1,Nmeasure-1);
for j = 1:nSources  %singleton
    lb(j) = 0; %lower bound
    ub(j) = min(measure(upperindex{j})); %upper bound
end

for j = (nSources+1) : (Nmeasure-nSources-1)
    lb(j) = max(measure(lowerindex{j})); %lower bound
    ub(j) = min(measure(upperindex{j})); %upper bound
end

for j = (Nmeasure-nSources):(Nmeasure-1) %(nSources-1)- tuple
    lb(j) = max(measure(lowerindex{j})); %lower bound
    ub(j) = 1;%upper bound
end
    subsetInterval = ub - lb;   
end