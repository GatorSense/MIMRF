function [y] = computeci(Bags,measure)
% Compute the Choquet integral output given data cells "Bags" and measure "measure"
% Can be used for multiple sources (>=3)
% This is correct! Takes in bag values. Fast in matrix forms.
% 
% A faster implementation can be seen in "computeci.c" (mex file) for
% multiple data points.
%
% INPUT
%    Bags       - 1 x nbags cell  - in each cell nPntPerBag x nSources double - data value  
%    measure  - 1 x (2^nSources-1) double - a fuzzy measure
%
% OUTPUT
%   y         - 1xN double  - the Choquet integral output for the data point
%
% Written by: X. Du 02/29/2016

%%

[nPntsBags, ~] = cellfun(@size, Bags);
nSources = size(Bags{1}, 2); %number of sources
n_bags = numel(Bags);


% set up MeasureSection (nchoosek)
for jjj = 1:nSources
    MeasureSection.NumEach(jjj) = nchoosek(nSources,jjj);%compute the cumulative number of measures of each tier. E.g. singletons,2s, 3s,..
    MeasureSection.Each{jjj} =  nchoosek([1:nSources],jjj);
end
MeasureSection.NumCumSum = cumsum(MeasureSection.NumEach);

% %Precompute differences and indices for each bag
diffM = cell(n_bags,1);
for i = 1:n_bags
    [v, indx] = sort(Bags{i}, 2, 'descend');
    vz = horzcat(zeros(size(v,1), 1), v) - horzcat(v, zeros(size(v,1), 1));
    diffM{i} = vz(:, 2:end);
    for j = 1:(size(diffM{i},2)-1) %# of sources in the combination (e.g. j=1 for g_1, j=2 for g_12)
        for n = 1:size(Bags{i},1)
            tmp{i,j}(n,:) = sort(indx(n,1:j), 2);
        end
    end
end

n_sources = size(Bags{1},2);
sec_start_inds = zeros(n_sources,1);
nElem_prev = 0;
for j = 1:(n_sources-1)
    if j == 1 %singleton        
        sec_start_inds(j) = 0;
    else  %non-singleton
        nElem_prev = nElem_prev+MeasureSection.NumEach(j-1);
        sec_start_inds(j) = nElem_prev;
    end
end

bag_row_ids = cell(n_bags,1); %the row indices of  measure used for each bag
for i = 1:length(Bags)
    [nPnts1, nSources] = size(diffM{i});
    bag_row_ids{i} = zeros(nPnts1,nSources-1);
    
    for n = 1:size(Bags{i},1)                
        for j = 1:(nSources-1) 
            if j == 1
                bag_row_ids{i}(n,j) = tmp{i,j}(n,1);
            else  %non-singleton
                elem = MeasureSection.Each{j};%the number of combinations, e.g., (1,2),(1,3),(2,3)
                [~,~,row_id] = ismember_findrow_mex_my(tmp{i,j}(n,:),elem);                
                bag_row_ids{i}(n,j) = sec_start_inds(j) + row_id;
            end
        end
    end
end
%Create oneV cell matrix
oneV = cell(n_bags,1);
for i  = 1:n_bags
    oneV{i} = ones(nPntsBags(i), 1);
end

y=[];
%%
singletonLoc = (nPntsBags == 1);

if sum(singletonLoc) %if there are singleton bags
diffM_s = vertcat(diffM{singletonLoc});
bag_row_ids_s = vertcat(bag_row_ids{singletonLoc});
% labels_s = Labels(singletonLoc)';

%Compute CI for singleton bags
ci = sum(diffM_s.*horzcat(measure(bag_row_ids_s),ones(sum(singletonLoc),1)),2);
y = vertcat(y,ci);
end

diffM_ns = diffM(~singletonLoc);
bag_row_ids_ns = bag_row_ids(~singletonLoc);
oneV = oneV(~singletonLoc);

%Compute CI for non-singleton bags
for i = 1:length(diffM_ns)
    ci = sum(diffM_ns{i}.*horzcat(measure(bag_row_ids_ns{i}), oneV{i}),2);
    y = vertcat(y,ci);
end

end

