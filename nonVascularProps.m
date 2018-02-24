%Compute and return the non vasular properties. 
clear all; 
close all; 
addpath('H:\Code\Hessian');
addpath('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\frangi_filter_version2a\');

%path1 = '\\zusdqbsfs006.usdqb.zeiss.org\advshare\Data\IschemiaTesting\DR_DensityResults';
%path1 = '\\zusdqbsfs006.usdqb.zeiss.org\advshare\Data\IschemiaTesting\Normals';
%path1 = 'C:\\Users\\u6abhatt\\Documents\\Work\\Code\\AriMatlab\\data';
%path1  = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\normal';
%path1 = 'N:\Portugal_Normal_DR\OCTA_SRL';
path1 = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\test';
%path1 = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\diseased2Repeats\';
%path1 = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\test';
%path1 = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\super';
path1 = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\portugal\OCTA_SRL\';
path1 = 'C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\portugal\portugalSmall\'

xx=11;
alpha = 1 ;
paths1 = allImagePaths(path1,'3x3');
paths2 = allImagePaths(path1,'3mm');
paths = [paths1 paths2];
fileID = fopen(strcat(path1,'\train_data_id.txt'),'w');
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
    %ComputeNonVascularProps(I, outpath, fname, id );
    stats=ComputeSkeleton2 ( I, outpath, fname, id);
    if id ==1
        totalstats = stats ;
    else
        totalstats = vertcat(totalstats, stats);
    end
    
    
    fprintf(fileID,'%s %s\n', num2str(id), fname);
end
% write total stats info 
writetable(struct2table(totalstats), 'AreaInfo.txt')

[B,ind] = sort( transpose([totalstats.Area]));
newdata = transpose(struct2cell(totalstats));

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
%{
Img   =imresize(Img,[512 512],'bicubic');
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
%}
%I = (histeq(Img));
%I = imguidedfilter(I);
%I = double(I);
ComputeHessian(Img, outpath, fname, id );
end

function ComputeHessian ( img , outpath, fname, id)
    
    imgCopy = img; 
    img = (histeq(img));
    img = imguidedfilter(img);
    img = double(img);
    Options.FrangiScaleRange = [1 6];
    Options.FrangiScaleRatio = 1;
    sigmas=Options.FrangiScaleRange(1):Options.FrangiScaleRatio:Options.FrangiScaleRange(2);
    %sigmas1=sigmas(1:1:20); %(0.1 to 2)
    %sigmas2=sigmas(51:1:100);%(5 to 10)
    %sigmas4=sigmas(131:1:end);%(13 to 16)
    % trying the existing code; 
    
    %H = 512;
    %W = 512;
    %I=imresize(img, [H W]);
   
    IM=sort(nonzeros(img(:)));
    NoiseLevel=mean(IM(end-4:end))*.11;
    %dims.mmX = 3;
    %dims.mmY = 3;
%calculate binary images
[BinaryImage, preHessianImage] = VasculatureBinarization(img, NoiseLevel);

    %existing code try end; 
    
    %sigmas1=sigmas(1:0.5:20);
    %sigmas2=sigmas(21:1:50);
    %sigmas4=sigmas(21:1:40);
    sigmas1 = [1:1:4];
    sigmas2 = [21:1:50];
    sigmas4 = [21:1:40];
    
    img=(img-NoiseLevel)/(max(img(:))-NoiseLevel);
    %smooth 
    %sigma=0.5;    % scale parameter in Gaussian kernel
    %G=fspecial('gaussian',12,sigma); % Caussian kernel
    %img=conv2(img,G,'same');  % smooth image by Gaussiin convolution
    
	[outIm2, whatScale2,Direction2] = Hessian_Vesselness(1-img,Options,sigmas2); %figure; imshow(outIm);title('outIm(5-10)');%figure; imshow(Direction);
	[outIm4, whatScale4,Direction4]  = Hessian_Vesselness(1-img,Options,sigmas4); %figure; imshow(outImL);title('outImL(13-16)');%figure; imshow(Direction);
	[outIm1, whatScale1,Direction1]  = Hessian_Vesselness(1-img,Options,sigmas1); %figure; imshow(outImh);title('outImh');%figure;imshow(Direction);
	[outImR, whatScale,Direction]    = Hessian_Vesselness(img,Options,1); %figure; imshow(outImR);title('outImR');%figure; imshow(Direction);
     
    %imLast (:,:) = max(outIm2(:,:),outIm4(:,:));
    %imLast (:,:) = max(imLast(:,:),outIm1(:,:));
    %imLast2 (:,:) = max(imLast(:,:),outImR(:,:)); 
    %figure; imshow([outImR,outIm2,outIm4]);
    level1 = adaptthresh(mat2gray(outIm1));
    level2 = adaptthresh(mat2gray(outIm2));
    level4 = adaptthresh(mat2gray(outIm4));
    
    %figure; imshow([whatScale1, whatScale2, whatScale4]);
     %figure; imshow([Direction1, Direction2, Direction4]);
    
    imBinary1 = imbinarize(mat2gray(outIm1),level1);
    %imBinary2 = imbinarize(mat2gray(outIm2),level2+0.3);
    imBinary2 = imbinarize(mat2gray(outIm1),level1+0.3);
    
    
    %imBinary1=bwareaopen(imBinary1,10);
    %imBinary2=bwareaopen(imBinary2,10);
    %figure; imshow([imBinary1, imBinary2,imBinary4]);
    
    %imSkel1 = bwmorph(imBinary1,'open');
    %imSkel2 = bwmorph(imBinary2,'open');
    %imSkel4 = bwmorph(imBinary4,'open');
    %figure; imshow([imSkel1, imSkel2,imSkel4]);
    
    imSkel1 =  bwmorph(bwmorph(imBinary1,'thin', inf),'clean');
    imSkel2 =  bwmorph(bwmorph(imBinary2,'thin', inf),'clean');
    %imSkel4 =  bwmorph(bwmorph(imBinary4,'thin', inf),'clean');
    %figure; imshow([imSkel1, imSkel2,imSkel4]);
    
    se = strel('disk',2);
    %imSkel11 = imdilate(imSkel1, se);
    imSkel2 = imdilate(imSkel2, se);
    imSkel2 = bwareaopen(imSkel2,50);
    %imSkel44 = imdilate(imSkel4, se);
    %figure; imshow([imSkel11, imSkel22, imSkel44]);
    %extend lines 
    %imSkel11 =  bwmorph(bwmorph(imSkel1,'diag',2),'thin', inf);
    %imSkel22 =  bwmorph(bwmorph(imSkel2,'diag',2),'thin', inf);
    %imSkel44 =  bwmorph(bwmorph(imSkel4,'diag',2),'thin', inf);
    %figure; imshow([imSkel11, imSkel22,imSkel44]);title('diag');
    
    %'bridge'
    %extend lines 
    %imSkel111 =  bwmorph(bwmorph(imSkel11,'bridge',2),'thin', inf);
    %imSkel222 =  bwmorph(bwmorph(imSkel22,'bridge',2),'thin', inf);
    %imSkel444 =  bwmorph(bwmorph(imSkel44,'bridge',2),'thin', inf);
    %figure; imshow([imSkel111, imSkel222,imSkel444]);title('bridge');
     %imSkel111 =  bwmorph(imSkel111,'close');
     %imSkel222 =  bwmorph(imSkel222,'close');
     %imSkel444 =  bwmorph(imSkel444,'close');
     %figure; imshow([imSkel111, imSkel222,imSkel444]);title('thick');
    %im = max(imSkel111,imSkel222);
    %im = max(im,imSkel444);
    %im = bwmorph(im, 'fill');
    %figure; imshow(im);
    %im = bwmorph(im, 'thin',inf);
    %figure; imshow(im);
    %figure; imshow([img,im]);
    %{
    imjoined = max(imSkel1,imSkel2);
    imjoined = bwmorph(imjoined, 'thin', inf);
    
    imBinaryD = imbinarize(mat2gray(outIm1));
    imBinaryD = bwmorph(imBinaryD, 'thin', inf);
    D = bwdist(imBinaryD);%figure; imshow(D,[]); title('D');f
    imBinaryD = imbinarize(mat2gray(D),0.5);
    ind = find (imBinaryD>0);
    imjoined (ind)=0; 
    imjoined = imdilate(imjoined, se); 
    
    BWopenSkel3(:,:,1)  = imjoined;
    BWopenSkel3(:,:,2)  = imjoined;
    BWopenSkel3(:,:,3)  = imjoined;
    
    %}
    %BWopenSkel3(:,:,1)  = max (uint8(BinaryImage.Skeleton).*255, uint8(imSkel2).*255);
    %BWopenSkel3(:,:,2)  = max (uint8(BinaryImage.Skeleton).*255, uint8(imSkel2).*255);
    %BWopenSkel3(:,:,3)  =  max (uint8(BinaryImage.Skeleton).*255, uint8(imSkel2).*255);
    
    %skeleton  = uint8(BinaryImage.Skeleton).*128;
    %skeleton(find(uint8(imSkel2).*128 == 128))=0;
    %BWopenSkel3(:,:,1) = uint8(imSkel2).*128;
    %BWopenSkel3(:,:,2) = skeleton;
    %BWopenSkel3(:,:,3) = 0;
    
    
    %skeleton  = BinaryImage.Skeleton;
    skeleton  = imSkel1;
    skeleton(find(imSkel2 == 1))=0;
    BWopenSkel3(:,:,1) = imSkel2;
    BWopenSkel3(:,:,2) = skeleton;
    BWopenSkel3(:,:,3) = 0;
    
    % add scale information 
    scaleSpace = uint8(mat2gray(whatScale1.*imSkel2).*255);
    finalSkeleton(:,:,1) = scaleSpace;
    finalSkeleton(:,:,2) = uint8(skeleton.*255);
    finalSkeleton(:,:,3) = 0;
    
    %for i=1:20000
    %    x = int16(rand*255)+1;
    %    y = int16(rand*255)+1;
    %    BWopenSkel3(x,y,1) =1; 
    %    BWopenSkel3(x,y,2) =0;
    %    BWopenSkel3(x,y,3) =0; 
    %end    
    
    img3(:,:,1)  = imgCopy;
    img3(:,:,2)  = imgCopy;
    img3(:,:,3)  = imgCopy;
    %wideImage = [uint8(img3), uint8(BWopenSkel3).*255];
    wideImage = [img3, finalSkeleton];
    %figure; imshow([img3, finalSkeleton]);
    %outnamePix2pix = strcat('C:\Users\u6abhatt\Documents\Work\Code\pix2pix-tensorflow\retina\train3b\', strcat(num2str(id),'.png'));
    outnamePix2pix = strcat('C:\Users\u6abhatt\Documents\Work\Code\pix2pix-tensorflow\retina\val3\', strcat(num2str(id),'.png'));
    imwrite((wideImage), convertStringsToChars(outnamePix2pix));
    
    %UINT TEST
    %{
    outnameTrainA = strcat('C:\Users\u6abhatt\Documents\Work\Code\UINT\skel2super\trainA\', strcat(num2str(id),'.png'));
    outnameTrainB = strcat('C:\Users\u6abhatt\Documents\Work\Code\UINT\skel2super\trainB\', strcat(num2str(id),'.png'));
    imwrite(double(img3), convertStringsToChars(outnameTrainA));
    imwrite(double(BWopenSkel3), convertStringsToChars(outnameTrainB));
    %}
    
    %{
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

    %}
end 

function stats=ComputeSkeleton2(Img, outpath, fname,id)
%I  = (imread ('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\test\PCZMI3968362_Angiography 3x3 mm_11-29-2016_13-0-18_OD_sn5572_FlowCube_z.img_original_super - Copy.png'));
%I  = (imread ('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\diseased2Repeats\CASE1-P610625SG_Angiography 3x3 mm_7-11-2016_10-36-55_OS_sn0370_super.png'));
%I = (imread ('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\normal\P5000-00006-DensityRR-313_Angiography 3x3 mm_3-4-2016_11-20-18_OD_sn0044_super.png'));
%I = imresize(I, [512, 512],'lanczos3');
%figure; imshow( I, []);

%figure; imshow(I,[]);

%IM=sort(nonzeros(I(:)));
%NoiseLevel=mean(IM(end-4:end))*.11;
%[BinaryImage, preHessianImage] = VasculatureBinarization(double(I), 0.1);
I = Img; 
I = double(histeq(I));
%I = padarray(I,[50 50],'both');
%I = double(imguidedfilter(I));


Options.FrangiScaleRange = [1 4];
Options.FrangiScaleRatio = 1;
Options.BlackWhite = false;
sigmas=Options.FrangiScaleRange(1):Options.FrangiScaleRatio:Options.FrangiScaleRange(2);
[outIm1, whatScale1,Direction1] = Hessian_Vesselness(I,Options,sigmas);
angles = sin(deg2rad(abs(Direction1)));

%figure; imshow(1-outIm1);
%figure; imshow(exp((1-outIm1).^3),[]);
%figure; imshow(imbinarize(mat2gray(exp((1-outIm1).^3))),[]);title('3');colormap(jet);
%figure; imshow(imbinarize(mat2gray(exp((1-outIm1).^2))),[]);title('2');colormap(jet);

%remove small things
invI = imbinarize(mat2gray(exp((1-outIm1).^2)));
%figure; imshow(~invI);title('~invI');
withoutsmall1  = bwareaopen(~invI, 50);
%figure; imshow(withoutsmall1);title('removed small');

% remove small blacks;
% this shows a nice set of surrounding edges
%{
invwithoutsmall1 = imcomplement(withoutsmall1); figure; imshow(invwithoutsmall1);
SE = strel('disk',2);
invwithoutsmall1SE = imerode(invwithoutsmall1,SE);figure; imshow(invwithoutsmall1SE);
invwithoutsmall1SE = imdilate(invwithoutsmall1,SE);figure; imshow(invwithoutsmall1SE);
diffI = imabsdiff(invwithoutsmall1SE, invwithoutsmall1);
figure; imshow(diffI);
%}


imSkel1 =  bwmorph(bwmorph(withoutsmall1,'thin', inf),'clean');
%figure; imshow(imSkel1);title('imskel');



SE = strel('disk',2);
imSkel1Extend = imdilate(imSkel1,SE);
%figure; imshow(imSkel1Extend,[]);title('imskelExtend');
%figure; imshow(imSkel1Extend,[]);

%remove intersection
imSkel1WOintersection = bwareaopen(1-imSkel1Extend, 100);% figure; imshow(imSkel1WOintersection);

intersectionMask = ~imabsdiff(imSkel1WOintersection,imSkel1Extend);
%figure; imshow(intersectionMask);title('intersectionMask');
%figure; imshow(imfuse(intersectionMask, imSkel1Extend),[]);
%figure; imshow(I.*intersectionMask,[]);
addMask = (I.*intersectionMask) > 150; % IMP I range goes from 0-255 check this./ 
% add the mask skeleton4 has no intersection/
skeleton4 = imSkel1Extend + addMask;
%figure; imshow(skeleton4);
%redo with the new one 
skeleton4Thin =  bwmorph(bwmorph(skeleton4,'thin', inf),'clean');
%figure; imshow(skeleton4Thin);
 outnamePix2pix = strcat('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\portugal\smalltestsGAN\', strcat(num2str(id),'.png'));
 imwrite((uint8(skeleton4Thin.*255)), convertStringsToChars(outnamePix2pix));


output ( :,:,1) = uint8(I.*0.5) + uint8(skeleton4Thin.*255); 
output (:,:,2)  = uint8(I.*0.5);
output (:,:,3)  = uint8(I.*0.5);
 outnamePix2pix = strcat('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\portugal\smalltests\', strcat(num2str(id),'.png'));
 imwrite((uint8(output)), convertStringsToChars(outnamePix2pix));
%figure; imshow( output,[]);
%{
%%Comparison%%
diskSize = 51;
[binaryImage_SRL, smoothMap_SRL, additionalImages] = VascularDensityCalculateSRL(I, diskSize, diskSize, -1);
%figure; imshow ( binaryImage_SRL.Skeleton);
comp ( :,:,1) = uint8(I.*0.5) + uint8( binaryImage_SRL.Skeleton.*255); 
comp (:,:,2)  = uint8(I.*0.5);
comp (:,:,3)  = uint8(I.*0.5);
%figure; imshow( comp,[]);
 outnamePix2pix = strcat('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\portugal\comparison\', strcat(num2str(id),'.png'));
 imwrite( [(uint8(output)),(uint8(comp))], convertStringsToChars(outnamePix2pix));
%}
SE = strel('disk',1);
skeleton4Extend = imdilate(skeleton4Thin,SE);
[L,NUM] = bwlabeln(1-skeleton4Extend);
[m,n] = size( skeleton4Extend );
test2 = zeros( m);
image = mat2gray(exp((1-outIm1).^2)); 
%figure; imshow ( L,[]);


for ii=1:NUM
individualMask = L==ii;
individualMask = imfill(individualMask,'holes');
numberOfTruePixels = sum(individualMask(:));
if(numberOfTruePixels > 400 )
grayvals = uint8(I);
grayvalInvert = uint8(ones (m).*255);
grayvalInvert = grayvalInvert - grayvals;
grayvalInvert = mat2gray(grayvalInvert);
grayvalInvert = mat2gray( exp(grayvalInvert));
%figure; imshow (grayvalInvert,[]);

    
Iindi = grayvalInvert.*individualMask;
%figure; imshow ( Iindi,[]);
Ibin = imbinarize(Iindi);
%figure; imshow (Ibin,[]);
test2 = test2 | Ibin;
end
end
% cluster 
cc = bwconncomp(test2);  
L = labelmatrix(cc);



indexNums = [];
fnames    = [];
for ii=1:cc.NumObjects
indexNums = vertcat(indexNums , ii);
fnames = vertcat(fnames , {fname});
end
stats = regionprops(cc, 'Area','Centroid');

C = num2cell(indexNums);
[stats(:).index] = deal(C{:});

C = num2cell(fnames);
[stats(:).fn] = fnames{:};



RGB2 = label2rgb(L, 'gray'); 
close all;
figure, imshow(RGB2);title  ('Area');
%figure; imshow ( imfuse(RGB2,imfuse(skeleton4Thin,I),'blend','Scaling','joint'));
hold on
for k = 1:numel(stats)
    x = stats(k).Centroid(1);
    y = stats(k).Centroid(2);
    text(x, y, sprintf('%d, %d', stats(k).Area,k), 'Color', 'g', ...
        'FontWeight', 'bold');
end
saveas(gcf,strcat('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\portugal\smalllabelImages\',fname));

end