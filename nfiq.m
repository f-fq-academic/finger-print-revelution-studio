function [ imgQualityScore ] = nfiq( inputAddress )
%NFIQ Summary of this function goes here
%   Detailed explanation goes here

[~,~,ext] = fileparts(inputAddress);

%Check if image format is not wsq so it's raw image
if(strcmp(ext,'.wsq')==0)
    %Get Image Info
    info = imfinfo(inputAddress);
    %Create Command to run nfiq algorithm
    fname = strcat('G:\Thesis\NBISBuild\bin\nfiq.exe -raw',{' '},int2str(info.Width),',',int2str(info.Height),',',int2str(info.BitDepth),{' "'},inputAddress,'"');
else
    fname = strcat('G:\Thesis\NBISBuild\bin\nfiq.exe ',{' "'},inputAddress,'"');
   
end

[status,imgQualityScore] = system(fname{1:1});

%if nfiq.exe reutrn error then we couldnt identify image quality scroe fot this image
if(status~=0)
    imgQualityScore = -1;
end

end

