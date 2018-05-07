function [TestConfMap] = computeTestMap(Bags, Labels, measure, Seg)
% Test stage of the MIMRF algorithm
% Input: 
%     - Bags (test): in cell form, multi-resolution. 1xNumBags cell. 
%               Within each Bag, NumInstances x NumSources cell.
%               Within each instance cell of the bag, contains an array of
%               corresponding points from all sources.
%     - Labels (test): Bag-level labels for test
%     - measure: the fuzzy measure learned from training. 1x(2^C-1) double
%     where C is the number of sources to be fused.
%     - Seg: this is specific to this demo data set. We used SLIC
%     segmentation to generate bags from the image and the Seg info is just
%     used to determine the index so as to map the test image back.
% Output:
%     - TestConfMap: ImgSize1 x ImgSize2 double. Can be plotted as an
%     image. This is the fusion result for each pixel of the image.
%
% Written by: X. Du 05/2018
%
nSources = size(Bags{1},2);
for ib = 1:numel(Bags)
    y_ci = zeros(1,size(Bags{ib},1));
    idx_all = find(Seg==ib-1); %all the pixel indices in each segment
        for j = 1:size(Bags{ib},1)
            hx = combvec(Bags{ib}{j,1},Bags{ib}{j,2});
            for jj = 3:nSources
            hx = combvec(hx,Bags{ib}{j,jj});
            end
            hx_unique = unique(hx','rows');
            Bags_hx{1} = hx_unique;
            y =  computeci(Bags_hx,measure);
            %          [aaa,bbb]=sort(y,'descend');
            if numel(y)==1 %only single point combination
                Locb(j,ib)=1;
                y_ci(j) = y;
            else
            %%%%%%%%%%%Use min-max operators to select the ``correct'' points
                if Labels(ib)==1
                    [y_ci(j),Locb(j,ib)] = max(y);
                else
                    [y_ci(j),Locb(j,ib)] = min(y);
                end
            end
        end
    TestConfMapr(idx_all) = y_ci';
end
    TestConfMap = reshape(TestConfMapr,[size(Seg,1),size(Seg,2)]);

end