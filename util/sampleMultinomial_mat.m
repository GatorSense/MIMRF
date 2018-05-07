function [outputIndex] = sampleMultinomial_mat(PDFdistr, NumSamples, MethodStr )
% Generate random samples from a multinomial. Return index of the sample.
%
% INPUT
% PDFdistr - distribution, can be fitness or interval or any distribution
%            input. Must be 1x#of distrbution values
% NumSamples - Number of samples
% MethodStr - ='descend', sort by descending order, the larger the distribution
% the more likely it is to be sampled. ='ascend', sort by ascending order,
% the smaller the distribution value the more likely it is to be sampled
% (to be used in picking samples). Case insensitive!
% OUTPUT
% Indxccc - index of the samples
%
% Written by: X. Du 03/2018
%

% %% Method 1: correct! Vector form. Uses a for-loop. Slow.
% [PDFdistr_sortV,PDFdistr_sortIdx] = sort(PDFdistr,'descend');
% PDFdistr_resort = PDFdistr(PDFdistr_sortIdx);
% PDFdistr_resort_norm = PDFdistr_resort./sum(PDFdistr_resort);
% PDFdistr_resort_normCDF = cumsum(PDFdistr_resort_norm);
%             ccc = rand(1,NumSamples);
%     for pciter = 1:NumSamples
%         Indxccc(pciter) = find(ccc(pciter)<=PDFdistr_resort_normCDF,1); %draw from Multinomial CDF, the first index/element that is larger than ccc =rand(1)
%     end
%
%    outputIndex =  PDFdistr_sortIdx(Indxccc);
%% Method 2: correct! Matrix form. Faster
if strcmpi(MethodStr,'descend')
    if  sum(PDFdistr<0)==0 %all positive distribution values
        [~,PDFdistr_sortIdx] = sort(PDFdistr,'descend');
        PDFdistr_resort = PDFdistr(PDFdistr_sortIdx);
    elseif  sum(PDFdistr<0)==numel(PDFdistr) %all negative
        PDFdistr = -1./PDFdistr;
        [~,PDFdistr_sortIdx] = sort(PDFdistr,'descend');
        PDFdistr_resort = PDFdistr(PDFdistr_sortIdx);
    else
        PDFdistr = PDFdistr-max(PDFdistr)-eps;
        PDFdistr = -1./PDFdistr;
        [~,PDFdistr_sortIdx] = sort(PDFdistr,'descend');
        PDFdistr_resort = PDFdistr(PDFdistr_sortIdx);
    end
elseif strcmpi(MethodStr,'ascend')
    if  sum(PDFdistr<0)==0 %all positive distribution values
        PDFdistr = -1./PDFdistr;
        [~,PDFdistr_sortIdx] = sort(PDFdistr,'ascend');
        PDFdistr_resort = PDFdistr(PDFdistr_sortIdx);
    elseif  sum(PDFdistr<0)==numel(PDFdistr) %all negative     
        [~,PDFdistr_sortIdx] = sort(PDFdistr,'ascend');
        PDFdistr_resort = PDFdistr(PDFdistr_sortIdx);
    else
        PDFdistr = PDFdistr-max(PDFdistr)-eps;
        [~,PDFdistr_sortIdx] = sort(PDFdistr,'ascend');
        PDFdistr_resort = PDFdistr(PDFdistr_sortIdx);
    end   
else
    error('Please input a sorting order: ''descend'' or ''ascend''!') 
end

PDFdistr_resort_norm = PDFdistr_resort./sum(PDFdistr_resort);
PDFdistr_resort_normCDF = cumsum(PDFdistr_resort_norm);
ccc = rand(1,NumSamples);
cccrepmat = repmat(ccc',1,numel(PDFdistr));
PDFdistr_resort_normCDFrepmat = repmat(PDFdistr_resort_normCDF,size(ccc,2),1);
diffrepmat = cccrepmat - PDFdistr_resort_normCDFrepmat;
diffb = (diffrepmat<0);
[~,Indxccc]= max( diffb~=0, [], 2 );
outputIndex =  PDFdistr_sortIdx(Indxccc);
end