%Compute and return the non vasular properties. 
clear all; 
close all; 
addpath('H:\Code\Hessian');
addpath('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\frangi_filter_version2a\');

%path1 = '\\zusdqbsfs006.usdqb.zeiss.org\advshare\Data\IschemiaTesting\DR_DensityResults';
%path1 = '\\zusdqbsfs006.usdqb.zeiss.org\advshare\Data\IschemiaTesting\Normals';
%path1 = 'C:\\Users\\u6abhatt\\Documents\\Work\\Code\\AriMatlab\\data';
%path1  = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\normal';
path1 = 'N:\Portugal_Normal_DR\OCTA_SRL';
path1 = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\totaldata';
path1 = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\test';
path1 = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\super';

xx=11;
alpha = 1 ;
paths = allImagePaths(path1,'3x3');

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
    outpath = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\weightTests\';
    outname = strcat(outpath,fname);
    %imwrite(I,outname);
    ComputeNonVascularProps(double(I), outpath, fname, id );
end
    



function [allpaths] = allImagePaths (path, pattern)
    global xx; 
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

function ComputeNonVascularProps(Img, outpath, fname,id)
% Return the properties. 
% List of descriptors without the FAZ 
% Print out the associated images.
Img   =imresize(Img,[256 256]);
NoiseLevel=min(Img(:));
Img    =(Img-NoiseLevel)/(max(Img(:))-NoiseLevel);
sigma  =1.5;    % scale parameter in Gaussian kernel
G      =fspecial('gaussian',9,sigma); % Gaussian kernel
ImgSmooth=conv2(1-Img,G,'same');  % smooth image by Gaussiin convolution
Y_tile=10;
X_tile=10;
Is=ImgSmooth;
Ih=(adapthisteq(1-Img,'NumTiles',[Y_tile X_tile]));
IH=((Img*2)+(1-Ih))/3;
ComputeHessian(IH, outpath, fname, id );
end

function ComputeHessian ( img , outpath, fname, id)
    global xx
    Options.FrangiScaleRange = [0.1 16];
    Options.FrangiScaleRatio = 0.1;
    sigmas=Options.FrangiScaleRange(1):Options.FrangiScaleRatio:Options.FrangiScaleRange(2);
    %sigmas1=sigmas(1:1:20); %(0.1 to 2)
    %sigmas2=sigmas(51:1:100);%(5 to 10)
    %sigmas4=sigmas(131:1:end);%(13 to 16)
    sigmas1=sigmas(1:0.5:5);
    sigmas2=sigmas(5:1:10);
    sigmas4=sigmas(10:1:20);
	[outIm, whatScale,Direction] = Hessian_Vesselness(1-img,Options,sigmas2); %figure; imshow(outIm);title('outIm(5-10)');%figure; imshow(Direction);
	[outImL,whatScale,Direction] = Hessian_Vesselness(1-img,Options,sigmas4); %figure; imshow(outImL);title('outImL(13-16)');%figure; imshow(Direction);
	[outImh,whatScale,Direction] = Hessian_Vesselness(1-img,Options,sigmas1); %figure; imshow(outImh);title('outImh');%figure;imshow(Direction);
	[outImR,whatScale,Direction] = Hessian_Vesselness(img,Options,1); %figure; imshow(outImR);title('outImR');%figure; imshow(Direction);
     
     IR3=outImR;
     HessianThres=0;
     Lochessian = find(IR3<=HessianThres);
     Hochessian = find(IR3>HessianThres);
     IR3(Lochessian) = 1; %figure; imshow(IR3);title('IR3 1');
     IR3(Hochessian) = 0; %figure; imshow(IR3);title('IR3 2');
	outIm=max(outIm,outImL); %figure; imshow(outIm);title('max(outIm,outImL)');
	Im=max(outIm.*(outIm>0.4*(max(outIm(:)))),outImh); %figure; imshow(Im);title('Im');
	I3=Im;
	II=Im;
	IIM=Im(Im>0);
	TH=median(IIM(:));
    I3=Im>0.6*TH;
    %figure; imshow(I3,[]);
    I33=bwareaopen(I3,1);
    IR3=IR3.*I33;
	I3=I3.*IR3;
	I3=bwareaopen(I3,9);
	BW=I3;
	BW=1-BW;
    %figure, imshow(BW);title('BW');
	BW=bwareaopen(BW,9);
	BW=1-BW;
    %figure; imshow(BW);title('1-BW');
	BW=bwmorph(BW,'bridge',1);
	BW=1-BW;
	BW=bwareaopen(BW,9);
	BW=1-BW;
	BW=bwmorph(BW,'diag',10);
	BW=1-BW;
	BW=bwareaopen(BW,5);
	BW=1-BW;
% 	BW=bwmorph(BW,'bridge',1);
	BW=bwmorph(BW,'diag',10);
	BW=bwareaopen(BW,15);
%	figure(12), imshow(BW,[0,1],'Border','Tight'), title('Adaptive Threshold plus Hessian filter')
	BW3 = bwmorph(BW,'skel',Inf);
    
    BWopenSkel3(:,:,1)  = BW3;
    BWopenSkel3(:,:,2)  = BW3;
    BWopenSkel3(:,:,3)  = BW3;
    
	BWDi=bwmorph(BW3,'dilate',1);
	BW3=BWDi;
	BW3=bwmorph(BW3,'bridge',1);
	BW3=bwmorph(BW3,'diag',10);
	BW3=bwareaopen(BW3,9);
	BW3 = bwmorph(BW3,'skel',Inf);
	BW3=bwmorph(BW3,'dilate',1);
	BW3=bwmorph(BW3,'bridge',10);
	BW3=bwmorph(BW3,'diag',10);
	BW3=bwareaopen(BW3,9);

	BWopen=bwareaopen(BW3,2048);
    
    BWopenSkel = bwmorph(BWopen,'skel',Inf);
    
    
    img3(:,:,1)  = img;
    img3(:,:,2)  = img;
    img3(:,:,3)  = img;
    wideImage = [img3, BWopenSkel3];
    %%
    %global xx;
    global alpha
    %outnamePix2pix = strcat('C:\Users\u6abhatt\Documents\Work\Code\pix2pix-tensorflow\retina\train\', strcat(num2str(id),'.png'));
    outnamePix2pix = strcat('C:\Users\u6abhatt\Documents\Work\Code\pix2pix-tensorflow\retina\train2\', strcat(num2str(id),'.png'));
    xx=xx+1;
    imwrite(wideImage, convertStringsToChars(outnamePix2pix));
    %%
	Isolate=BW3-BWopen;
    %This is a distance map, which can be used for weighting not for contouring.
	D = bwdist(BWopen);%figure; imshow(D,[]); title('D');f
    
    
    Dopen=D>20;
	if sum(Dopen(:))>625
			Dopen=bwareaopen(Dopen,625);
	else
% 			Dopen=D>15;

			Dopen=bwareaopen(Dopen,36);

	end
	
	Dopen=bwmorph(Dopen,'dilate',2);
	Dopen=1-bwareaopen(1-Dopen,8100);

	Isolate=and(Dopen,Isolate);
	
	BW(Isolate==1)=0;
	BW=bwareaopen(BW,16);
 
	BBW=bwmorph(BW,'dilate',3);
	BBW=1-BBW;
	BBW=bwmorph(BBW,'dilate',1);
	BBW=bwareaopen(BBW,100);
	
	BBI(:,:,1)=img+BBW/2;
	BBI(:,:,2)=img;
	BBI(:,:,3)=img;
    %figure; imshow(BBI,[]);title ('BBI');
    
    %Compute characteristics from the image 
    DMat = mat2gray(D); 
    %figure; imshow(DMat, []);
    
    %find  large regions near FAZ 
    se = strel('disk',2);
    FAZbase = imopen(~BW,se); %figure; imshow(FAZbase,[]);
    
    labeledImage = bwlabel(FAZbase, 4);  
    blobMeasurements = regionprops(labeledImage, 'Area');
    allBlobAreas = [blobMeasurements.Area];
    largestBlobIndex = (allBlobAreas >= max(allBlobAreas));
    largestIndex = find(largestBlobIndex);
    keeperBlobsImage = ismember(labeledImage, largestIndex);

    se = strel('disk',2);
    keeperBlobsImageT = imerode(keeperBlobsImage,se);
    %remove the FAZ 
    CC = bwconncomp(keeperBlobsImageT, 4);
    S  = regionprops(CC, 'Area');
    L  = labelmatrix(CC);
    keeperBlobsImage = ismember(L, find([S.Area] <max([S.Area])));
    keeperBlobsImage = imdilate(keeperBlobsImage,se);
    %figure; imshow(BW2);



    level = adaptthresh(DMat);
    imTh2  = imbinarize (DMat, level);
    %figure; imshow(imTh2,[]);
    
    se = strel('disk',2);
    imTh2afterOpening = imopen(imTh2,se);
    %figure; imshow(imTh2afterOpening,[]);
    %imTh2afterOpening(keeperBlobsImage)=1;
    %nothing should over lap the vessels.
    imTh2afterOpening(BW)=0;
    %Determine the connected components:
    CC = bwconncomp(imTh2afterOpening, 4);
    %Compute the area of each component:
    S = regionprops(CC, 'all');
    numberOfBlobs = size(S, 1);
    
    %Remove small objects:
    L = labelmatrix(CC);
    imTh2afterOpening = ismember(L, find([S.Area] <max([S.Area])));
    %REcompute everything after removing center 
    CC = bwconncomp(imTh2afterOpening, 4);
    %Compute the area of each component:
    S = regionprops(CC, 'all');
    numberOfBlobs = size(S, 1);
    %Remove small objects:
    L = labelmatrix(CC);
    
    %imTh2afterOpeningLargest20 = ismember(L, find([S.Area] >= P));
    coloredLabels = label2rgb (L, 'hsv', 'k', 'shuffle'); % pseudo random color labels
    %%
    %Generate images
    %figure;
    %imshow(BW);
    %title('Outlines, from bwboundaries()'); 
    %axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
    %hold on;
    %boundaries = bwboundaries(imTh2afterOpening);
    %numberOfBoundaries = size(boundaries, 1);
    %for k = 1 : numberOfBoundaries
    %    thisBoundary = boundaries{k};
    %    plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
    %end
    %hold off;
    %%%
    %figure; imshow(coloredLabels,[]);
    
    
    %color by area
    %init
    outim = zeros(size(L));
    sArea = [S.Area];
    uniqueArea = unique(sArea);
    
    
    %for  k = 1:(length(uniqueArea)-1) % removing the largest value (FAZ)
    for  k = ((length(uniqueArea)-1)-100):((length(uniqueArea))) %top 100
        if k > 1
        ar                   = uniqueArea(k);
        allBlobAreas         = sArea;
        allowableAreaIndexes = (allBlobAreas == ar);
        keeperIndexes        = find(allowableAreaIndexes);
        keeperBlobsImage = ismember(L, keeperIndexes);
        outim(keeperBlobsImage) = ar;
        
        end
    end
   
    normalizedWeights = (mat2gray(outim));
   
    
    outname = strcat(outpath,"\");
    outname = strcat(outname, fname);
    
    indimage = uint8(normalizedWeights.*128);
    
    rgbImage = ind2rgb(indimage, jet);
    
    
    BBI(:,:,2)=img/2 + double(indimage);
	BBI(:,:,3)=img/2;
	BBI(:,:,1)=img/2;
    imwrite(BBI, convertStringsToChars(outname));;
    [filepath,name,ext] = fileparts(fname);
    outname =  strcat(outpath,"\");
    outname =  strcat(outname, name);
    outname1 =  strcat(outname, '_NormalizedArea.png');
    outname2 =  strcat(outname, '_AllArea.png');
    imwrite(rgbImage, convertStringsToChars(outname1));
    imwrite(coloredLabels, convertStringsToChars(outname2));

    
    %%% weights %% by distance. 
    weights1 = zeros(size(L));
    [m, n]  = size (weights1); 
    weights1(m/2,n/2)=1; 
    Dist = mat2gray(bwdist(weights1,'euclidean'));
    for ii=1:m
        for jj=1:n
            %indimage(ii,jj) = indimage(ii,jj)*(1-Dist(ii,jj))+img(ii,jj).*255;
            if indimage(ii,jj) > 1
            BBI2(ii,jj,1) = (1-Dist(ii,jj))+img(ii,jj);
            else 
            BBI2(ii,jj,1) = img(ii,jj);
            end;
            BBI2(ii,jj,2) = img(ii,jj);
            
            if indimage(ii,jj) > 1
            BBI2(ii,jj,3) = (Dist(ii,jj))+img(ii,jj);
            else 
            BBI2(ii,jj,3) = img(ii,jj);
            end;
        end
    end
    outname4 =  strcat(outname, '_distanceFromCenter.png');
    imwrite( BBI2,convertStringsToChars(outname4));
    
    
    options = struct(...
    'FrangiScaleRange',[4 10],...
    'FrangiScaleRatio',1,...
    'FrangiBetaTwo', 15,...
    'BlackWhite',      false);
[J,Scale,Direction]=FrangiFilter2D(Im, options);
J=imadjust(J);
%figure; imshow(J);
largevesselsBase=imbinarize(J);
%figure;imshow(largevesselsBase);
%remove small parts
se = strel('disk',2);
largevessels = imerode(largevesselsBase,se);
%figure; imshow(largevessels);
%%%%  distance from large vessels %%%%
Dist = mat2gray(bwdist(largevessels,'euclidean'));
    for ii=1:m
        for jj=1:n
            %indimage(ii,jj) = indimage(ii,jj)*(1-Dist(ii,jj))+img(ii,jj).*255;
            if indimage(ii,jj) > 1
            BBI3(ii,jj,1) = (Dist(ii,jj))+img(ii,jj);
            else 
            BBI3(ii,jj,1) = img(ii,jj);
            end
            BBI3(ii,jj,2) = img(ii,jj);
            
            if indimage(ii,jj) > 1
            BBI3(ii,jj,3) = (1-Dist(ii,jj))+img(ii,jj);
            else 
            BBI3(ii,jj,3) = img(ii,jj);
            end
        end
    end
    outname5 =  strcat(outname, '_distanceFromLargeVessels.png');
    imwrite( BBI3,convertStringsToChars(outname5));
%%%%  distance from large vessels %%%%
    %figure; imshow(imfuse(outim,img));
    %%% Orientation of the top 100%%%%
    CC = bwconncomp(normalizedWeights, 4);
    %Compute the area of each component:
    S = regionprops(CC, 'Orientation');
    %Remove small objects:
    L = labelmatrix(CC);
    outAngle = zeros(size(L));
    sAngle = [S.Orientation];
    for i=1:length(sAngle)
        angle = sAngle(i);
        allBlobAngles         = sAngle;
        allowableAngleIndexes = (allBlobAngles == angle);
        keeperIndexes         = find(allowableAngleIndexes);
        keeperBlobsImage = ismember(L, keeperIndexes);
        outAngle(keeperBlobsImage) = angle;
    end
    normalizedAngles = (mat2gray(outAngle));
    indimage = uint8(normalizedWeights.*128);
    rgbImage = ind2rgb(indimage, jet);
    outname3 =  strcat(outname, '_AllAngle.png');
    imwrite(rgbImage, convertStringsToChars(outname3));

    
end 