function [y] = ChoquetIntegral_g_MultiSources(hx,g,nSources)
% Compute the Choquet integral output given a data point "hx" and measure "g"
% Can be used for multiple sources (>=3)
% This is correct! But only computes one data point at a time. If multiple
% data points, need to do a for-loop.
% 
% A faster implementation can be seen in "computeci.c" (mex file) for
% multiple data points.
%
% INPUT
%    hx       - 1 x nSources double  - data value  
%    g        - a fuzzy measure
%    nSources - number of sources
%
% OUTPUT
%   y         - 1x1 double  - the Choquet integral output for the data point
%
% Written by: X. Du 08/10/2015

%%
% sort the source values for data point hx
[vals, sort_inds] = sort(hx, 'descend');
y = 0;

%for first value
y = y + (vals(1) - vals(2))*g(sort_inds(1));

%for middle values
nElem_prev = 0;
sort_ind_temp = sort_inds(1);
 for i = 2:nSources-1 %sample rest
        nElem = nchoosek(nSources,i);%the number of combinations, e.g.,3
        elem = nchoosek([1:nSources], i);%the number of combinations, e.g., (1,2),(1,3),(2,3)
             val_all = [];
             tmp_all = [];
       nElem_prev = nElem_prev+nchoosek(nSources,i-1); 
        for j = 1:nElem
            elemSub = nchoosek(elem(j,:), i-1);       
            tmp = num2cell(elem(j,:));
            valg = g(nElem_prev+j);
              
            tmp_all = vertcat(tmp_all,cell2mat(tmp));
            val_all = horzcat(val_all,valg);
        end
        
        sort_ind_temp = horzcat(sort_ind_temp,sort_inds(i)); % this is the index of g, say g_124
        
        [vals_sort_ind_temp, sorted_sort_ind_temp] = sort(sort_ind_temp, 'ascend'); %in order to satisfy g_12 == g_21 for example
        
        [~, Locb] = ismember(vals_sort_ind_temp, tmp_all,'rows');
        g_temp = val_all(Locb);
        y = y + (vals(i) - vals(i+1))*g_temp;
 end

%for end value
y = y + (vals(end) - 0)*g(end);

end

