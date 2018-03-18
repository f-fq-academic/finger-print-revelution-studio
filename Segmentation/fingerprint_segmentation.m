function F3=fingerprint_segmentation(H)

% Matlab Code for paper (version 1.5): 
% [1] M. F. Fahmy, and M. A. Thabet, "A fingerprint segmentation technique
% based on Morphological processing," ISSPIT 2013.
% 
% Code written by M. A. Thabet.
% 
% Code released for **research purposes only**
% 
% Feel free to modify/distribute but please cite [1].
% 
% Contact m.a.t.bishay@qmul.ac.uk


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to segment fingerprint image using morphological processing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=im2double(H);

%%% Step 1: Feature Extraction
H1=rangefilt(H); 

%%% Convert to binary image
H2=adaptivethreshold(H1,16,0.05,0);

%%% Step 2: Morphological Processing
SE = strel('disk', 6, 4);
H3=imclose(~H2,SE);
H4=imopen(H3,SE);
H4 = imfill(H4,'holes');
SE1=strel('square',12);
H4=imerode(H4,SE1);
[m1,n1] = size(H);
F=zeros(m1,n1);
[L, num] = bwlabel(H4,8);
max1=0;
for j=1:num
[r, c] = find(bwlabel(L)==j);
m=size(r,1);
if (max1<m)
    id=j;
    max1=m;
end
end
[r, c] = find(bwlabel(L)==id);
s3=size(r,1);
for i=1:s3
    F(r(i),c(i))=1;
end

%%% Step 3: Contour Smoothing
erf=boundaries(F);
erf=erf{1};
z=frdescp(erf);
z1=ifrdescp(z,50);
x = round(z1(:, 1)); 
y = round(z1(:, 2));
if(min(x)<=0)
      x = x - min(x) + 1;
end
if (min(y)<=0)
      y = y - min(y) + 1;
end

ZZ=zeros(m1,n1);
for ii=1:size(z1,1)
    ZZ(x(ii),y(ii))=1;
end

SE1=strel('square',3);
ZZ=imdilate(ZZ,SE1);

[r1, c1] = find(ZZ==1);
ZZ=imfill(ZZ,'holes');

for i=1:m1
  for j=1:n1
      if(ZZ(i,j)==1)
          F3(i,j)=H(i,j);
      else
          F3(i,j)=0.5;
      end
  end
end



function bw=adaptivethreshold(IM,ws,C,tm)

if (nargin<3)
    error('You must provide the image IM, the window size ws, and C.');
elseif (nargin==3)
    tm=0;
elseif (tm~=0 && tm~=1)
    error('tm must be 0 or 1.');
end

IM=mat2gray(IM);

if tm==0
    mIM=imfilter(IM,fspecial('average',ws),'replicate');
else
    mIM=medfilt2(IM,[ws ws]);
end
sIM=mIM-IM-C;
bw=im2bw(sIM,0);
bw=imcomplement(bw);

function B = boundaries(BW, conn, dir)	

if nargin < 3	
   dir = 'cw';
end	
if nargin < 2	
   conn = 8;
end
L = bwlabel(BW, conn);
numObjects = max(L(:));
if numObjects > 0
   B = {zeros(0, 2)};
   B = repmat(B, numObjects, 1);
else
   B = {};	
end

Lp = padarray(L, [1 1], 0, 'both');

M = size(Lp, 1);
if conn == 8	
   
   offsets = [-1, M - 1, M, M + 1, 1, -M + 1, -M, -M-1];
else	
   offsets = [-1, M, 1, -M];	
end	
	
if conn == 8	
   next_search_direction_lut = [8 8 2 2 4 4 6 6];	
else	
   next_search_direction_lut = [4 1 2 3];	
end	
	
if conn == 8	
   next_direction_lut = [2 3 4 5 6 7 8 1];	
else
   next_direction_lut = [2 3 4 1];	
end	
START    = -1;
BOUNDARY = -2;

scratch = zeros(100, 1);
	
[rr, cc] = find((Lp(2:end-1, :) > 0) & (Lp(1:end-2, :) == 0));
rr = rr + 1;	
for k = 1:length(rr)	
   r = rr(k);	
   c = cc(k);	
   if (Lp(r,c) > 0) & (Lp(r - 1, c) == 0) & isempty(B{Lp(r, c)})
      
      idx = (c-1)*size(Lp, 1) + r;	
      which = Lp(idx);	
      scratch(1) = idx;
      Lp(idx) = START;	
      numPixels = 1;	
      currentPixel = idx;
      initial_departure_direction = [];
      done = 0;	
      next_search_direction = 2;	
      while ~done	
         direction = next_search_direction;	
         found_next_pixel = 0;
         for k = 1:length(offsets)
            neighbor = currentPixel + offsets(direction);
            if Lp(neighbor) ~= 0	
               if (Lp(currentPixel) == START) & ...	
                      isempty(initial_departure_direction)
                  
                  initial_departure_direction = direction;            
               elseif (Lp(currentPixel) == START) & ...	
                      (initial_departure_direction == direction)
                 
                  done = 1;
                  found_next_pixel = 1;
                  break;
               end
               next_search_direction = ...	
                   next_search_direction_lut(direction);
               found_next_pixel = 1;
               numPixels = numPixels + 1;	
               if numPixels > size(scratch, 1)
                  scratch(2*size(scratch, 1)) = 0;
               end
               scratch(numPixels) = neighbor;
               if Lp(neighbor) ~= START
                  Lp(neighbor) = BOUNDARY;	
               end	
               currentPixel = neighbor;
               break;	
            end
            direction = next_direction_lut(direction);
         end	
         if ~found_next_pixel
           
            numPixels = 2;
            scratch(2) = scratch(1);	
            done = 1;
         end	
      end
      
      [row, col] = ind2sub(size(Lp), scratch(1:numPixels));
      B{which} = [row - 1, col - 1];	
   end	
end
if strcmp(dir, 'ccw')
   for k = 1:length(B)	
      B{k} = B{k}(end:-1:1, :);	
   end	
end

function z = frdescp(s)	

[np, nc] = size(s);	
if nc ~= 2 	
   error('S must be of size np-by-2.'); 
end	
if np/2 ~= round(np/2);
   s(end + 1, :) = s(end, :);	
   np = np + 1;	
end

x = 0:(np - 1);	
m = ((-1) .^ x)';		

s(:, 1) = m .* s(:, 1);	
s(:, 2) = m .* s(:, 2);	
s = s(:, 1) + i*s(:, 2);	
z = fft(s);

function s = ifrdescp(z, nd)	

np = length(z);	

if nargin == 1 | nd > np 
   nd = np; 	
end

x = 0:(np - 1);	
m = ((-1) .^ x)';	

d = round((np - nd)/2); 
z(1:d) = 0;	
z(np - d + 1:np) = 0;	

zz = ifft(z);	
s(:, 1) = real(zz);	
s(:, 2) = imag(zz);	

s(:, 1) = m.*s(:, 1);	
s(:, 2) = m.*s(:, 2);