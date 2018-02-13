clear all; 
close all; 
%  = double(imread ('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\test\octa_digital_3x3_phantom_v00_super.png'));
%I  = double(imread ('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\test\PCZMI3968362_Angiography 3x3 mm_11-29-2016_13-0-18_OD_sn5572_FlowCube_z.img_original_super.png'));
I  = (imread ('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\test\PCZMI3968362_Angiography 3x3 mm_11-29-2016_13-0-18_OD_sn5572_FlowCube_z.img_original_super - Copy.png'));
I = (histeq(I));
I = double(imguidedfilter(I));
%.FrangiScaleRange : The range of sigmas used, default [1 8]
%       .FrangiScaleRatio : Step size between sigmas, default 2
%       .FrangiBetaOne : Frangi correction constant, default 0.5
%       .FrangiBetaTwo : Frangi correction constant, default 15
%       .BlackWhite : Detect black ridges (default) set to true, for
%                       white ridges set to false.

Options.FrangiScaleRange = [1 4];
Options.FrangiScaleRatio = 1;
Options.BlackWhite = false;
sigmas=Options.FrangiScaleRange(1):Options.FrangiScaleRatio:Options.FrangiScaleRange(2);
[outIm2, whatScale2,Direction2] = Hessian_Vesselness(I,Options,sigmas);
angles = sin(deg2rad(abs(Direction2)));
%figure; imshow(outIm2,[]), title('out'); colormap(jet);
%figure; imshow(Direction2,[]);title('Dirc');colormap(jet);
%figure; imshow(whatScale2,[]);title('scale');colormap(jet);

angles(find(outIm2<0.1))=-1;
whatScale2(find(outIm2<0.1))=-1;
%figure; imshow(angles,[]);title('th angles');
%figure; imshow(whatScale2,[]);title('th sclae');
%edge detection 
BW = edge(outIm2,'Canny');
figure; imshow(BW);title('canny');
total = BW + outIm2; 
totalGray = mat2gray(total); 
figure; imshow(totalGray);

imBiniary = imbinarize(totalGray,0.1);
imBiniary = bwareaopen(imBiniary, 150);
figure; imshow(~imBiniary); title ( 'Binary3');
D = bwdist(imBiniary); figure; imshow(D,[]);
DBinary = imbinarize(D,5);
figure; imshow(DBinary,[]);title('Binary');

figure; imshow(imfuse(DBinary, BW),[]);
% cluster based on properties. 
cc = bwconncomp(DBinary); 
L = labelmatrix(cc);
pixelIndexList = label2idx(L);
stats = regionprops(cc, 'Area','Orientation','Eccentricity','MajorAxisLength','MinorAxisLength'); 

[idx,C] = kmeans(transpose(cell2mat(struct2cell(stats))),5);
%%%%
% = DBinary;
[m n] = size(DBinary);
cluster1  = zeros(m);
cluster2  = zeros(m);
cluster3  = zeros(m);
cluster4  = zeros(m);
cluster5  = zeros(m);

for i=1:length(idx)
    id = idx(i);
    if (id==1)
    cluster1(cell2mat(pixelIndexList(i)))=id+1;
    end
    if (id==2)
    cluster2(cell2mat(pixelIndexList(i)))=id+1;
    end
    if (id==3)
    cluster3(cell2mat(pixelIndexList(i)))=id+1;
    end
    if (id==4)
    cluster4(cell2mat(pixelIndexList(i)))=id+1;
    end
    if (id==5)
    cluster5(cell2mat(pixelIndexList(i)))=id+1;
    end
end


figure; imshow([cluster1, cluster2, cluster3,cluster4,cluster5],[]); title ('clustered');
figure; imshow(imfuse(~imBiniary, cluster1),[]);
figure; imshow(imfuse(~imBiniary, cluster2),[]);
figure; imshow(imfuse(~imBiniary, cluster3),[]);
