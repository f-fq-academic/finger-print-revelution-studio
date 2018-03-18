clear all;
close all;
clc;

%inputAddress = 'G:\Thesis\RawDatabase\sample_102_1.tif';
inputAddress = '/home/fingerproject/Desktop/FingerprintFolder/Databases/FVC2000/DB1_B/101_8.tif';

%Quality Estimation
%imgQualityScore = nfiq(inputAddress);
%imgQualityScore




%read image from input
Iinitial = imread(inputAddress);

% make image gray-scale
% Get the number of rows and columns,
% and, most importantly, the number of color channels.
[rows, columns, numberOfColorChannels] = size(Iinitial);
if numberOfColorChannels > 1
    % It's a true color RGB image.  We need to convert to gray scale.
    Iinitial = rgb2gray(Iinitial);
end



show(Iinitial,1,'Original Image');

%Image Segmentation
%segmentedImage = fingerprint_segmentation(Iinitial);
%clear Iinitial;   %remove the original image from the memory
%show(segmentedImage,2);

% Identify ridge-like regions and normalise image
blksze = 16; thresh = 0.1; margin=20;
[normim, mask] = ridgesegment(Iinitial, blksze, thresh, margin);
clear segmentedImage;   %remove the original image from the memory


% Determine ridge orientations
[orientim, reliability] = ridgeorient(normim, 1, 5, 5);
%plotridgeorient(orientim, 20, im, 2)
show(reliability,3,'Reliablity')

% Determine ridge frequency values across the image
blksze = 36;
[freq, medfreq] = ridgefreq(normim, mask, orientim, blksze, 5, 5, 15);
show(freq,4,'Frequency Image');

% Actually I find the median frequency value used across the whole
% fingerprint gives a more satisfactory result...
freq = medfreq.*mask;
%midfreq can be appox to 0.1

newim = ridgefilter(normim, orientim, freq, 0.5, 0.5, 0);
show(newim,5,'Enhanced Image');

% Binarise, ridge/valley threshold is 0
binim = newim > 0;
clear newim;
%binim = newim > 150;
%binim = adaptiveThres(newim,32,0);

show(binim,6,'Binary Image');

% Display the skeletonized fingerprint from the binary image
thinim = ridgethin(1-binim, margin);
show(thinim.*mask.*(reliability>0.5),7,'Thin image');
clear binim;

%display all minutae
%1st stage - find all the minutae
[ridgeEnd, ridgeBifurcation, ridgerOrderMap] = findminutae(thinim, orientim, mask.*(reliability>0.5));
%show_minutia(normim, ridgeEnd, ridgeBifurcation);

%2nd stage - remove spurious minutae
[ridgeEnd, ridgeBifurcation] = remove_spurious_minutae(ridgeEnd, ridgeBifurcation, ridgerOrderMap);


% print out x,y,theta(indegrees) to console
minutae = [ridgeEnd;ridgeBifurcation];
% [minutae(:,[1,2])';(minutae(:,3)*180/pi)']'  %displays


show_minutia(thinim.*mask.*(reliability>0.5), ridgeEnd, ridgeBifurcation, 'Fingerprint Minutaes');
%show_minutia(thinim, ridgeEnd, ridgeBifurcation, 'Fingerprint Minutaes');


%Calculate distance of each point from each other
counter=1;
size = (length(minutae)*(length(minutae)-1)/2) +1;
dist =size:5;

for k = 1 : length(minutae)
    for m=1:length(minutae)
        if(k~=m)
        dist(counter,1) = minutae(k,1);
        dist(counter,2) = minutae(k,2);
        dist(counter,3) = minutae(k,3);
        dist(counter,4) = minutae(m,1);
        dist(counter,5) = minutae(m,2);
        dist(counter,6) = minutae(m,3);
        dist(counter,7)=  sqrt( (minutae(k,2) - minutae(m,2))^2 + (minutae(k,1) - minutae(m,1))^2);
        counter = counter+1;
        end;
    end
end


%Create figure to show all distance calculated lines
figure(1025);
plot(minutae(:,2),minutae(:,1),'om',...
    'MarkerSize',5,...
    'LineWidth',1);
hold on;

for x=1:length(dist)
    plot([dist(x,2) dist(x,5)],[dist(x,1) dist(x,4)]);
    hold on;
end

max(dist(:,7))
mean(dist(:,7))
min(dist(:,7))

figure(1026);
plot(minutae(:,2),minutae(:,1),'om',...
    'MarkerSize',5,...
    'LineWidth',1);
hold on;


%Calculate Degree between All Two Lines
%arccos((AB^2+CB^2-AC^2)/(2*AB*CB))
% counterOfLines=1;
%sizeOfLines = (length(minutae)*(length(minutae)-1)/2) +1;
% %degs;
% 
% for k = 1 : length(dist)
%     for m=k+1:length(dist)
%         CBIndex = find(ismember(dist,[dist(k,4) dist(k,5) dist(k,6) dist(m,4) dist(m,5) dist(m,6)],'R2012a'),1);
%         
% %         AB = dist(k,7)
% %         CB = dist(CBIndex,7)
% %         AC = dist(m,7)
%         
%        degs(counterOfLines,1) = acosd(dist(k,7)^2+dist(CBIndex,7)^2-dist(m,7)^2)/(2*dist(k,7)*dist(CBIndex,7));
%         counterOfLines = counterOfLines+1;
%     end
%     
% end
% 
% counterOfLines


