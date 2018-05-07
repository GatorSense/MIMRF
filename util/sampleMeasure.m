function [measure] = sampleMeasure(nSources, lowerindex, upperindex,  prevSample, IndexsubsetToChange, sampleVar)
%sampleMeasure - sampling a new measure
% If the number of inputs are two, sampling a brand new measure (only used during initialization)
% If the number of inputs are all five, sampling a new value for only one element of the measure
% This code samples a new measure value from truncated Gaussian distribution.
%
% INPUT
%   nSources - number of sources
%   lowerindex - the cell that stores all the corresponding subsets (lower index) of measure elements
%   upperindex - the cell that stores all the corresponding supersets (upper index) of measure elements
%   prevSample (optional)- previous measure, before update
%   IndexsubsetToChange (optional)- the measure element index to be updated
%   sampleVar (optional) - sampleVar set in Parameters
% OUTPUT
%   - measure - new measure after update
%
% Written by: X. Du 03/2018
%

if(nargin < 4)
    % sample a brand new measure
    % Flip a coin to decide whether to sample from above or sample from bottom
    %sample new measure, prior is uniform/no prior
    coin = rand(1);
    if coin >= 0.5 % sample from bottom
        [measure] = sampleMeasure_Bottom(nSources,lowerindex);
    else   % sample from above
        [measure] = sampleMeasure_Above(nSources,upperindex);
    end
else
    %sample just one new element of measure
    
    Nmeasure = 2^nSources-1;
    measure = prevSample;
    if (IndexsubsetToChange<=nSources) && (IndexsubsetToChange>=1) %singleton
        lowerBound = 0; 
        upperBound = min(measure(upperindex{IndexsubsetToChange})); 
    elseif (IndexsubsetToChange>=(Nmeasure-nSources)) && (IndexsubsetToChange<=(Nmeasure-1)) %(nSources-1)-tuple
        lowerBound = max(measure(lowerindex{IndexsubsetToChange})); 
        upperBound = 1;
    else  %remaining elements
        lowerBound = max(measure(lowerindex{IndexsubsetToChange})); 
        upperBound = min(measure(upperindex{IndexsubsetToChange})); 
    end
    
    denom = upperBound - lowerBound;
    v_bar = sampleVar/(denom^2+eps);%-changed
    x_bar = measure(IndexsubsetToChange); %-changed
    
    %% sample from a Truncated Gaussian
    sigma_bar = sqrt(v_bar);
    c2 = rand(1);       % randomly generate a value between [0,1]
    [val] = invcdf_TruncatedGaussian (c2,x_bar,sigma_bar,lowerBound,upperBound); % new measure element value
    measure(IndexsubsetToChange) = val;  % new measure element value
end

end