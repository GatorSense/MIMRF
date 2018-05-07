
function [fitness] = evalFitness_minmax_multires(Labels, measure, n_bags, nPntsBags, oneV, bag_row_ids, diffM)
% Evaluate the fitness a measure, using min( sum(max((ci-0)^2)) + sum(min(ci-1)^2) ) for classification.
%
% INPUT
%    Labels         - 1xNumTrainBags double  - Training labels for each bag
%    measure        -  measure to be evaluated after update
%    n_bags         - number of bags
%    nPntsBags      - 1xNumTrainBags double    - number of points in each bag
%    oneV           - cell matrix n_bags x max(nPntsBags), convenient for computing CI output 
%    bag_row_ids    - the row indices of  measure used for each bag
%    diffM          - Precompute differences for each bag
% OUTPUT
%   fitness         - the fitness value using min(sum(min((ci-d)^2))) for regression.
% 
% Written by: X. Du 05/2018
%

fitness = 0;

%Compute CI 
for i = 1:n_bags
    for j = 1:nPntsBags(i) %for all points in each bag
    y = sum(diffM{i,j}.*horzcat(measure(bag_row_ids{i,j}), oneV{i,j}),2);
    
     %%%%%%%%%%%min-max pick points
        if Labels(i)==1
            [y_ci(j),Locb(j,i)] = max(y); 
        else
            [y_ci(j),Locb(j,i)] = min(y);
        end
       
    end
    
    if Labels(i)~=1 %negative bag label=0
        fitness = fitness - max(y_ci.^2);
    else
        fitness = fitness - min((y_ci-1).^2);
    end   
         
end

   % sanity check
    if isinf(fitness) || ~isreal(fitness)
        fitness = real(fitness);
        if fitness == Inf
            fitness = -10000000;
        elseif fitness == -Inf
            fitness = -10000000;
        end
    end


    

end
