function [Bags, Labels, Seg_trainImg] = generateSimData_MU()
% This function generates a simulated multi-resolution, five-source,
% "MU" data set for the MIMRF algorithm.
% Based on a set of five pre-generated multi-resolution images, construct
% five-source multi-resolution training bags.
%
% Written by: X. Du 05/2018
%
%% Set up parameters to construct multi-resolution bags
Parameters.windowsize  = 1.5; %determine the pixel correspondence between multi-resolution images
Parameters.ifPlot = 1; %if=1, plot the fusion sources
Parameters.dsratio_Img1 = 1; % the downsampling ratio used to create each fusion source
Parameters.dsratio_Img2 = 2;
Parameters.dsratio_Img3 = 1;
Parameters.dsratio_Img4 = 2;
Parameters.dsratio_Img5 = 1;
%% Generate data set

% Read in an "MU" logo template
load('demo_MultiRes_data_MU.mat') % This is a pre-generated "MU" lettering patterns

ImgM = Img; %letter "M"
ImgM(39:end,:)=0;

ImgU = Img;  %letter "U"
ImgU(1:39,:)=0;



%%% Generate five images (as fusion sources) with different resolutions 
%%% (by downsamplng in different ratios)
%%% The images were pre-segmented by the SLIC algorithm (Achanta et al., 2012) is used to generate
%%% bags. Each superpixel is a bag (as marked with red boundaries in figure 100). 
%%% The labels for the bags are also generated. If a superpixel has any
%%% part of the "M" or "U" letter, it is considered positive and otherwise
%%% negative.
figure(100);set(gcf, 'Position', get(0, 'Screensize'));
%%%%% Image 1: High confidence on letter "M"
idx_Img1=find(Img1==1);
J1=Img1;
imwrite(J1,'murgb.PNG');
I = imread('murgb.PNG');
B = [];
Labels1 = [];
for i = 1:numel(unique(Seg1)) %number of Seg
    BW=double(Seg1==(i-1));
    BWSegize(i) = sum(sum(BW));
    [B{i},L{i}] = bwboundaries(BW,'noholes');
end
if Parameters.ifPlot==1
    subplot(1,5,1);imagesc(I);hold on
    for i = 1:length(B)
    visboundaries(B{i});hold on
    end
    axis equal tight
    title('Source 1')
end

for i = 1:length(B)
	if ~isempty(intersect(idx_Img1,find(Seg1==i-1)));
        Labels1(i) = 1;
    else
        Labels1(i) = 0;
    end
end
%sum(Labels1)

%%%%% Image 2: High confidence on letter "M"
idx_Img2=find(Img2==1);
J2=Img2;
imwrite(J2,'murgb.PNG');
I = imread('murgb.PNG');
B = [];
Labels2 = [];
for i = 1:numel(unique(Seg2)) %number of Seg
    BW=double(Seg2==(i-1));
    BWSegize(i) = sum(sum(BW));
    [B{i},L{i}] = bwboundaries(BW,'noholes');
end
if Parameters.ifPlot==1
    subplot(1,5,2);imagesc(I);hold on
    for i = 1:length(B)
    visboundaries(B{i});hold on
    end
    title('Source 2')
    axis equal tight
end


for i = 1:length(B)
	if ~isempty(intersect(idx_Img2,find(Seg2==i-1)));
        Labels2(i) = 1;
    else
        Labels2(i) = 0;
    end
end
%sum(Labels2)


%%%%% Image 3: High confidence on letter "U"
idx_Img3=find(Img3==1);
J3=Img3;
imwrite(J3,'murgb.PNG');
I = imread('murgb.PNG');
B = [];
Labels3 = [];
for i = 1:numel(unique(Seg3)) %number of Seg
    BW=double(Seg3==(i-1));
    BWSegize(i) = sum(sum(BW));
    [B{i},L{i}] = bwboundaries(BW,'noholes');
end
if Parameters.ifPlot==1
    subplot(1,5,3);imagesc(I);hold on
    for i = 1:length(B)
    visboundaries(B{i});hold on
    end
    title('Source 3')
    axis equal tight
end


for i = 1:length(B)
	if ~isempty(intersect(idx_Img3,find(Seg3==i-1)));
        Labels3(i) = 1;
    else
        Labels3(i) = 0;
    end
end
%sum(Labels3)

%%%%% Image 4: High confidence on letter "U"
idx_Img4=find(Img4==1);
J4=Img4;
imwrite(J4,'murgb.PNG');
I = imread('murgb.PNG');
B = [];
Labels4 = [];
for i = 1:numel(unique(Seg4)) %number of Seg
    BW=double(Seg4==(i-1));
    BWSegize(i) = sum(sum(BW));
    [B{i},L{i}] = bwboundaries(BW,'noholes');
end
if Parameters.ifPlot==1
    subplot(1,5,4);imagesc(I);hold on
    for i = 1:length(B)
    visboundaries(B{i});hold on
    end
    title('Source 4')
    axis equal tight
end


for i = 1:length(B)
	if ~isempty(intersect(idx_Img4,find(Seg4==i-1)));
        Labels4(i) = 1;
    else
        Labels4(i) = 0;
    end
end
%sum(Labels4)


%%%%% Image 5. High confidence on "background"
idx_Img5=find(Img5==0);
J5=Img5;
imwrite(J5,'murgb.PNG');
I = imread('murgb.PNG');
B = [];
Labels5 = [];
for i = 1:numel(unique(Seg5)) %number of Seg
    BW=double(Seg5==(i-1));
    BWSegize(i) = sum(sum(BW));
    [B{i},L{i}] = bwboundaries(BW,'noholes');
end
if Parameters.ifPlot==1
    subplot(1,5,5);imagesc(I);hold on
    for i = 1:length(B)
    visboundaries(B{i});hold on
    end
    title('Source 5')
    axis equal tight
end


for i = 1:length(B)
	if ~isempty(intersect(idx_Img5,find(Seg5==i-1)));
        Labels5(i) = 1;
    else
        Labels5(i) = 0;
    end
end
%sum(Labels5)

%%%%%%%%%%% Generate multi-resolution bags

Easting0 = 558264; %randomly generated numbers to simulate a geo-tagged image with easting northing info for each pixel
Northing0 = 4310377;

Image1.Easting  = Easting0:Parameters.dsratio_Img1:(Easting0+size(J1,2)-1);
Image1.Northing  = Northing0:Parameters.dsratio_Img1:(Northing0+size(J1,1)-1);
Image1.Data = J1;
Image1.Datar = reshape(J1,[size(J1,1)*size(J1,2),1]);
Image1.Seg = Seg1;
count = 1;
for i = 1:numel(Image1.Easting)
    for j = 1:numel(Image1.Northing)
        Image1.EN(count,1) =  Image1.Easting(i);
        Image1.EN(count,2) =  Image1.Northing(j);
        count = count+1;
    end
end

Image2.Easting  = Easting0:Parameters.dsratio_Img2:(Easting0+size(J1,2)-1);
Image2.Northing  = Northing0:Parameters.dsratio_Img2:(Northing0+size(J1,1)-1);
Image2.Data = J2;
Image2.Datar = reshape(J2,[size(J2,1)*size(J2,2),1]);
Image2.Seg = Seg2;
count = 1;
for i = 1:numel(Image2.Easting)
    for j = 1:numel(Image2.Northing)
        Image2.EN(count,1) =  Image2.Easting(i);
        Image2.EN(count,2) =  Image2.Northing(j);
        count = count+1;
    end
end

Image3.Easting  = Easting0:Parameters.dsratio_Img3:(Easting0+size(J1,2)-1);
Image3.Northing  = Northing0:Parameters.dsratio_Img3:(Northing0+size(J1,1)-1);
Image3.Data = J3;
Image3.Datar = reshape(J3,[size(J3,1)*size(J3,2),1]);
Image3.Seg = Seg3;
count = 1;
for i = 1:numel(Image3.Easting)
    for j = 1:numel(Image3.Northing)
        Image3.EN(count,1) =  Image3.Easting(i);
        Image3.EN(count,2) =  Image3.Northing(j);
        count = count+1;
    end
end

Image4.Easting  = Easting0:Parameters.dsratio_Img4:(Easting0+size(J1,2)-1);
Image4.Northing  = Northing0:Parameters.dsratio_Img4:(Northing0+size(J1,1)-1);
Image4.Data = J4;
Image4.Datar = reshape(J4,[size(J4,1)*size(J4,2),1]);
Image4.Seg = Seg4;
count = 1;
for i = 1:numel(Image4.Easting)
    for j = 1:numel(Image4.Northing)
        Image4.EN(count,1) =  Image4.Easting(i);
        Image4.EN(count,2) =  Image4.Northing(j);
        count = count+1;
    end
end

Image5.Easting  = Easting0:Parameters.dsratio_Img5:(Easting0+size(J1,2)-1);
Image5.Northing  = Northing0:Parameters.dsratio_Img5:(Northing0+size(J1,1)-1);
Image5.Data = J5;
Image5.Datar = reshape(J5,[size(J5,1)*size(J5,2),1]);
Image5.Seg = Seg5;
count = 1;
for i = 1:numel(Image5.Easting)
    for j = 1:numel(Image5.Northing)
        Image5.EN(count,1) =  Image5.Easting(i);
        Image5.EN(count,2) =  Image5.Northing(j);
        count = count+1;
    end
end

Bags = [];
%%%%%%%%%Method 1: corresponds to all pixels
for i = 1:(max(max(Image5.Seg))+1)
    idx_all = [];
    idx_all = find(Image5.Seg==i-1); %all the pixel indices in each segment
    [indexRowCol] = rowcol(idx_all,Image5.Seg);
    for j = 1:numel(idx_all)  
      easting_min =  Image5.Easting(indexRowCol(j,2));
      northing_min =  Image5.Northing(indexRowCol(j,1));
     idxtemp_Img1 =   find(Image1.EN(:,1)>=easting_min & Image1.EN(:,1)<=easting_min+Parameters.windowsize & Image1.EN(:,2)>=northing_min & Image1.EN(:,2)<=northing_min+Parameters.windowsize);
     idxtemp_Img2 =    find(Image2.EN(:,1)>=easting_min & Image2.EN(:,1)<=easting_min+Parameters.windowsize & Image2.EN(:,2)>=northing_min & Image2.EN(:,2)<=northing_min+Parameters.windowsize);
     idxtemp_Img3 =    find(Image3.EN(:,1)>=easting_min & Image3.EN(:,1)<=easting_min+Parameters.windowsize & Image3.EN(:,2)>=northing_min & Image3.EN(:,2)<=northing_min+Parameters.windowsize);
     idxtemp_Img4 =    find(Image4.EN(:,1)>=easting_min & Image4.EN(:,1)<=easting_min+Parameters.windowsize & Image4.EN(:,2)>=northing_min & Image4.EN(:,2)<=northing_min+Parameters.windowsize);

        if isempty(idxtemp_Img1)
            dist = pdist2(Image1.EN,[easting_min northing_min]);
            [~,idxtemp_Img1] = min(dist);
        end
        if isempty(idxtemp_Img2)
            dist = pdist2(Image2.EN,[easting_min northing_min]);
            [~,idxtemp_Img2] = min(dist);
        end
        if isempty(idxtemp_Img3)
            dist = pdist2(Image4.EN,[easting_min northing_min]);
            [~,idxtemp_Img4] = min(dist);
        end
        if isempty(idxtemp_Img4)
            dist = pdist2(Image4.EN,[easting_min northing_min]);
            [~,idxtemp_Img4] = min(dist);
        end

          Bags{i}{j,1} = Image1.Datar(idxtemp_Img1)' ;
          Bags{i}{j,2} = Image2.Datar(idxtemp_Img2)' ;
          Bags{i}{j,3} = Image3.Datar(idxtemp_Img3)' ;
          Bags{i}{j,4} = Image4.Datar(idxtemp_Img4)' ;
          Bags{i}{j,5} = Image5.Datar(idx_all(j)) ; 
    end 
end
Seg_trainImg = Image5.Seg;
Labels = Labels5; %we are generating labels based on the bags from Img5 

delete 'murgb.PNG'
end