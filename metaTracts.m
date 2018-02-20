clear all; 
close all; 
%  = double(imread ('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\test\octa_digital_3x3_phantom_v00_super.png'));
%I  = double(imread ('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\test\PCZMI3968362_Angiography 3x3 mm_11-29-2016_13-0-18_OD_sn5572_FlowCube_z.img_original_super.png'));
I  = (imread ('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\test\PCZMI3968362_Angiography 3x3 mm_11-29-2016_13-0-18_OD_sn5572_FlowCube_z.img_original_super - Copy.png'));
I = (histeq(I));
I = double(imguidedfilter(I));


%Generate paths 
%path1  = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\totaldata\';
%path1  = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\diseased\'
path1 = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\diseased2Repeats\'

outpath = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\metaTractsTests1\';
outpath = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\diseased2Repeats\garb\'

fileID = fopen('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\metaTractsTests1\train_data_id.txt','w');

%paths1 = allImagePaths(path1,'3mm');
%paths2 = allImagePaths(path1,'3x3');
%paths  = [paths1 paths2];

paths = allImagePaths(path1,'CASE1');


%make file paths 



%for each image 
for id=1:length(paths)
    dirpath = paths(id).folder; 
    fname   = paths(id).name; 
    pathname = strcat(dirpath,"\");
    pathname = strcat(pathname,fname);
    disp(pathname)
    %readimage
    I=imread(convertStringsToChars(pathname));
    if size(I,3)==3
        I = rgb2gray(I);
    end
    %ComputeNonVascularProps(double(I));
    outname = strcat(outpath,fname);
    %imwrite(I,outname);
    ComputeAvascular(I, outpath, fname, id );
    
    
    fprintf(fileID,'%s %s\n', num2str(id), fname);
end
function ComputeAvascular(I, outpath, fname, idit )

I = (histeq(I));
I = imguidedfilter(I);
J=I;
I = double(I);
Options.FrangiScaleRange = [1 4];
Options.FrangiScaleRatio = 1;
Options.BlackWhite = false;
sigmas=Options.FrangiScaleRange(1):Options.FrangiScaleRatio:Options.FrangiScaleRange(2);
[outIm2, whatScale2,Direction2] = Hessian_Vesselness(I,Options,sigmas);
angles = sin(deg2rad(abs(Direction2)));
%figure; imshow(outIm2,[]), title('out'); colormap(jet);
%figure; imshow(Direction2,[]);title('Dirc');colormap(jet);
%figure; imshow(whatScale2,[]);title('scale');colormap(jet);

%angles(find(outIm2<0.1))=-1;
whatScale2(find(outIm2<0.1))=0;
%figure; imshow(angles,[]);title('th angles');
%figure; imshow(whatScale2,[]);title('th sclae');
%edge detection 
BW = edge(outIm2,'Canny');
%figure; imshow(BW);title('canny');
total = BW + outIm2; 
totalGray = mat2gray(total); 
%figure; imshow(totalGray);

imBiniary = imbinarize(totalGray,0.1);
%imBiniary = bwareaopen(imBiniary, 150);
%figure; imshow(~imBiniary); title ( 'Binary3');

%write out the images

    img3(:,:,1)  = J;
    img3(:,:,2)  = J;
    img3(:,:,3)  = J;
    img1(:,:,1)  = uint8(~imBiniary).*255;
    img1(:,:,2)  = uint8(~imBiniary).*255;
    img1(:,:,3)  = uint8(~imBiniary).*255;
%wideImage = [img3, img1];
wideImage =[uint8(I), uint8(outIm2.*255), uint8(mat2gray(angles.*255)),...
    uint8(mat2gray(whatScale2.*255))];
outnamePix2pix = strcat(outpath, strcat(num2str(idit),'.png'));
%outnamePix2pix = strcat('C:\Users\u6abhatt\Documents\Work\Code\pix2pix-tensorflow\retina\val\', strcat(num2str(id),'.png'));
imwrite(wideImage, convertStringsToChars(outnamePix2pix));


D = bwdist(imBiniary); %figure; imshow(D,[]);
DBinary = imbinarize(D,5);
%figure; imshow(DBinary,[]);title('Binary');

%figure; imshow(imfuse(DBinary, BW),[]);
% cluster based on properties. 
cc = bwconncomp(DBinary); 
L = labelmatrix(cc);
pixelIndexList = label2idx(L);
stats = regionprops(cc, 'Area','Orientation','Eccentricity','MajorAxisLength','MinorAxisLength'); 

[idx,C] = kmeans(transpose(cell2mat(struct2cell(stats))),5);
%%%%
% = DBinary;
[m n] = size(DBinary);
allCluster = zeros(m); 
cluster1  = zeros(m);
cluster2  = zeros(m);
cluster3  = zeros(m);
cluster4  = zeros(m);
cluster5  = zeros(m);

for i=1:length(idx)
    id = idx(i);
    if (id==1)
    cluster1(cell2mat(pixelIndexList(i)))=id+1;
    allCluster(cell2mat(pixelIndexList(i)))=id+1;
    end
    if (id==2)
    cluster2(cell2mat(pixelIndexList(i)))=id+1;
    allCluster(cell2mat(pixelIndexList(i)))=id+1;
    end
    if (id==3)
    cluster3(cell2mat(pixelIndexList(i)))=id+1;
    allCluster(cell2mat(pixelIndexList(i)))=id+1;
    end
    if (id==4)
    cluster4(cell2mat(pixelIndexList(i)))=id+1;
    allCluster(cell2mat(pixelIndexList(i)))=id+1;
    end
    if (id==5)
    cluster5(cell2mat(pixelIndexList(i)))=id+1;
    allCluster(cell2mat(pixelIndexList(i)))=id+1;
    end
end
%% print the clustered image 
    img3(:,:,1)  = J;
    img3(:,:,2)  = J;
    img3(:,:,3)  = J;
    img1(:,:,1)  = uint8(mat2gray(allCluster)).*255;
    img1(:,:,2)  = uint8(mat2gray(allCluster)).*255;
    img1(:,:,3)  = uint8(mat2gray(allCluster)).*255;
wideImage = [img3, img1];
outnamePix2pix = strcat(strcat(outpath,"\\clustered\\"), strcat(num2str(idit),'.png'));
%outnamePix2pix = strcat('C:\Users\u6abhatt\Documents\Work\Code\pix2pix-tensorflow\retina\val\', strcat(num2str(id),'.png'));
imwrite(wideImage, convertStringsToChars(outnamePix2pix));


%{
figure; imshow([cluster1, cluster2, cluster3,cluster4,cluster5],[]); title ('clustered');
figure; imshow(imfuse(~imBiniary, cluster1),[]);
figure; imshow(imfuse(~imBiniary, cluster2),[]);
figure; imshow(imfuse(~imBiniary, cluster3),[]);
figure; imshow(imfuse(~imBiniary, cluster4),[]);
figure; imshow(imfuse(~imBiniary, cluster5),[]);
%}
end

function [allpaths] = allImagePaths (path, pattern)
    disp(['Find all images in the path- ',path])
    myFiles = dir(fullfile(path)); %gets all wav files in struct
    tmpPaths=[];
    for d=1:length(myFiles)
        %disp(['name d',num2str(d), myFiles(d).name])
        if ~strcmp(myFiles(d).name,'.') && ~strcmp( myFiles(d).name,'..') 
            if (myFiles(d).isdir) 
                %disp( myFiles(d));
                newpath =strcat(myFiles(d).folder,'\\');
                newpath =strcat(newpath, myFiles(d).name);
                tmpPaths = [tmpPaths allImagePaths(newpath, pattern)];
            else
                %disp(myFiles(d).name)
                [filepath,name,ext] = fileparts(myFiles(d).name);
                if contains(ext,'.png')
                    if pattern == ""
                        tmpPaths= [tmpPaths myFiles(d)];
                    else
                        if contains(myFiles(d).name, pattern)
                        % Many files have super in them
                        % looking for the ones t
                            if regexp(myFiles(d).name, '\w*super.png', 'match', 'once')
                                tmpPaths= [tmpPaths myFiles(d)];
                            end
                        end
                    end
                end
            end
            
        end
    end
    allpaths = tmpPaths;
end