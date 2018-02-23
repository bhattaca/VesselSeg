function [BinaryImage, smoothMap] = VascularDensityCalculateNN( image, noiseLevel, diskH, diskW, mmX, mmY )
    %% ARINDAM ALGO 
    %% G.A. tracker
    gaTracker(); %Portal tracker function, please dont delete. Contact Alex 510-542-7614 if this line gives you any problems!

    %%
    imageNET = NET.convertArray(image(:), 'System.Byte', length(image(:)));
    densityNETObj = Czm.Oct.CommonAlgorithms.Density.VasculatureDensity(imageNET, size(image, 1), size(image, 2), diskH, diskW, mmX, mmY);

    densityNETObj.CalculateEdge(noiseLevel);
    densityNETObj.Calculate(noiseLevel);
    
    outImage = densityNETObj.OutImage;
    outWidth = densityNETObj.OutWidth;
    outHeight = densityNETObj.OutHeight;
    
    BinaryImage.Perf = reshape(densityNETObj.Area.double, [outWidth, outHeight]);
    BinaryImage.Vessel = reshape(densityNETObj.Skeleton.double, [outWidth, outHeight]);
    BinaryImage.Edge = reshape(densityNETObj.Edge.double, [outWidth, outHeight]);
    BinaryImage.PerfForFAZ = reshape(densityNETObj.AreaForFAZ.double, [outWidth, outHeight]);
    
    smoothMap.RescaledGray = reshape(outImage.double, [outWidth, outHeight]);
    smoothMap.Perf = reshape(densityNETObj.AreaMap.double, [outWidth, outHeight]);
    smoothMap.Vessel = reshape(densityNETObj.LengthMap.double, [outWidth, outHeight]);
    smoothMap.Edge = reshape(densityNETObj.EdgeMap.double, [outWidth, outHeight]);
    smoothMap.ImageNormalized = reshape(densityNETObj.ImageNormalized.double, [outWidth, outHeight]);
    
    %%%
    I = image; 
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



imSkel1 =  bwmorph(bwmorph(withoutsmall1,'thin', inf),'clean');
%figure; imshow(imSkel1);title('imskel');



SE = strel('disk',2);
imSkel1Extend = imdilate(imSkel1,SE);
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
%imwrite((uint8(skeleton4Thin.*255)), convertStringsToChars(outnamePix2pix));


output ( :,:,1) = uint8(I.*0.5) + uint8(skeleton4Thin.*255); 
output (:,:,2)  = uint8(I.*0.5);
output (:,:,3)  = uint8(I.*0.5);
outnamePix2pix = strcat('C:\Users\u6abhatt\Documents\Work\Code\AriMatlab\data\portugal\smalltests\', strcat(num2str(id),'.png'));
%imwrite((uint8(output)), convertStringsToChars(outnamePix2pix));

SE = strel('disk',1);
skeleton4Extend = imdilate(skeleton4Thin,SE);
BinaryImage.Vessel = skeleton4Thin;
end

