% This code is used to generate the multiscale threshold corresponding to
% LSP

% It is from the thre_N,dp in the paper CMDSP: A completed multi-domain shrinkage pattern for texture image classification by Bin Li, Xiaochun Xu,
% which is submitted to PeerJ Computer Science



clear all;
clc;

% Add the path for the dataset such as Outex_TC_00020,UMD,UIUC
rootpic = 'E:\TextureClassification\outex\Outex_TC_00010\';
picNum = 4320;

radius = 2; 
neighbors = 8;

nh = 128-2*radius;
nl = 128-2*radius;

Outex10_MG28th3 =  zeros(nh,nl,picNum);

% Outex10_M28th3 =  zeros(nh,nl,picNum);
% Outex10_M2Ga28th3 =  zeros(nh,nl,picNum);
% Outex10_M2Ga316th3 =  zeros(nh,nl,picNum);
% Outex10_M18th3 =  zeros(nh,nl,picNum);
for i=1:picNum;
    filename = sprintf('%s\\images\\%06d.ras', rootpic, i-1);
    fprintf('No.%d\n',i);
    
    Gray = imread(filename); 

   
    Grayy = im2double(Gray);
   Gray = (Grayy-mean(Grayy(:)))/std(Grayy(:))*20+128; 


    % The multiscale threshold used to generate the center pixel of the gradient image
     Gradient_im2 = sobel2_grad0(Gray);
     Gray = (Gradient_im2-mean(Gradient_im2(:)))/std(Gradient_im2(:))*20+128;


   % % The multiscale threshold used to generate the center pixel of the entropy image
   %  h = [1,1,1,1,1;
   %      1,1,1,1,1;
   %      1,1,1,1,1;
   %      1,1,1,1,1;
   %      1,1,1,1,1];   
   %  Gray = entropyfilt(Gray,h);
% 
    % The multiscale threshold used to generate the center pixel of the Gaussian image
    image_scale =  Gray;
   imgExt = padarray(image_scale,[3 3],'symmetric','both');
   sigma = 2^0.25;
   scale = 4;
   xx = 2*ceil(2*sigma)+1;
   Image_S(:,:,1) = Gray;
   for gaussianconv = 1:(scale-1)
        h = fspecial('gaussian', 2*ceil(2*sigma)+1, sigma);
        Image_filter = imfilter(imgExt,h);
        Image_S(:,:,gaussianconv+1) = Image_filter(4:end-3,4:end-3);
        imgExt = Image_filter;
%         figure(gaussianconv+1);
%         imshow(imgExt,[]);     
   end
     Gray = Image_S(:,:,3);

%     
    [a,b] = size(Gray); 
    [ysize xsize] = size(Gray);
    
   
     spoints = zeros(neighbors,2); 
     ang = 2*pi/neighbors;
     for n = 1:neighbors
        spoints(n,1) = -radius*sin((n-1)*ang);
        spoints(n,2) = radius*cos((n-1)*ang);
     end
     miny=min(spoints(:,1));
     maxy=max(spoints(:,1));
     minx=min(spoints(:,2));
     maxx=max(spoints(:,2));

    bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;
    bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;

    origy=1-floor(min(miny,0));
    origx=1-floor(min(minx,0));

    dx = xsize - bsizex;
    dy = ysize - bsizey;

   bins = 2^neighbors;

   d_image = Gray;
   d_C = d_image(origy:origy+dy,origx:origx+dx);


th22_Gray= zeros((a-2*radius),(b-2*radius),neighbors);
th44_Gray_o= zeros(a,b,neighbors);
th44_Gray = zeros((a-2*radius),(b-2*radius),neighbors);
over_Gray = zeros((a-2*radius),(b-2*radius),neighbors);

for j = 1:neighbors
  y = spoints(j,1)+origy;
  x = spoints(j,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    N = d_image(ry:ry+dy,rx:rx+dx);
    D{j} = N >= d_C;   
%     Diff{j} = abs(N-d_C); 
    Diff{j} = 1./(1+exp(N-d_C)); 
    [h l]= size(Diff{j});
    
%   %2*2 local mean
    Gray_re22 = im2col(Diff{j},[(h/2) (l/2)],'distinct');
    Gray_re22_M  = mean(Gray_re22);
    Gray22_1 =Gray_re22_M(1).*ones((h/2),(l/2));
    Gray22_2 =Gray_re22_M(2).*ones((h/2),(l/2));
    Gray22_3 =Gray_re22_M(3).*ones((h/2),(l/2));
    Gray22_4 =Gray_re22_M(4).*ones((h/2),(l/2));
    th22_Gray(:,:,j) = [Gray22_1,Gray22_3;
                 Gray22_2,Gray22_4];  
    %4*4 local mean          
    Gray_Ext = padarray(Diff{j},[radius radius],'symmetric','both');        
    Gray_re44 = im2col(Gray_Ext,[(a/4) (b/4)],'distinct');
    Gray_re44_M  = mean(Gray_re44);
    Gray44_1 =Gray_re44_M(1).*ones((a/4),(b/4));
    Gray44_2 =Gray_re44_M(2).*ones((a/4),(b/4));
    Gray44_3 =Gray_re44_M(3).*ones((a/4),(b/4));
    Gray44_4 =Gray_re44_M(4).*ones((a/4),(b/4));
    Gray44_5 =Gray_re44_M(5).*ones((a/4),(b/4));
    Gray44_6 =Gray_re44_M(6).*ones((a/4),(b/4));
    Gray44_7 =Gray_re44_M(7).*ones((a/4),(b/4));
    Gray44_8 =Gray_re44_M(8).*ones((a/4),(b/4));
    Gray44_9 =Gray_re44_M(9).*ones((a/4),(b/4));
    Gray44_10 =Gray_re44_M(10).*ones((a/4),(b/4));
    Gray44_11 =Gray_re44_M(11).*ones((a/4),(b/4));
    Gray44_12 =Gray_re44_M(12).*ones((a/4),(b/4));
    Gray44_13 =Gray_re44_M(13).*ones((a/4),(b/4));
    Gray44_14 =Gray_re44_M(14).*ones((a/4),(b/4));
    Gray44_15 =Gray_re44_M(15).*ones((a/4),(b/4));
    Gray44_16 =Gray_re44_M(16).*ones((a/4),(b/4));
    th44_Gray_o(:,:,j) = [Gray44_1,Gray44_5,Gray44_9,Gray44_13;
                        Gray44_2,Gray44_6,Gray44_10,Gray44_14;
                        Gray44_3,Gray44_7,Gray44_11,Gray44_15;
                        Gray44_4,Gray44_8,Gray44_12,Gray44_16]; 
    th44_Gray(:,:,j) = th44_Gray_o((1+radius):(a-radius),(1+radius):(b-radius),j); 
  
    MeanDiff(j) = mean(mean(Diff{j}));
    over_Gray(:,:,j) = MeanDiff(j).*ones((a-2*radius),(b-2*radius));
%     
  else
    % Interpolation needed, use double type images 
    ty = y - fy;
    tx = x - fx;

    % Calculate the interpolation weights.
    w1 = (1 - tx) * (1 - ty);
    w2 =      tx  * (1 - ty);
    w3 = (1 - tx) *      ty ;
    w4 =      tx  *      ty ;
    % Compute interpolated pixel values
    N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
        w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    D{j} = N >= d_C;     
%     Diff{j} = abs(N-d_C);
    Diff{j} = 1./(1+exp(N-d_C)); 
    [h l]= size(Diff{j});
%    %2*2 local mean 
    Gray_re22 = im2col(Diff{j},[(h/2) (l/2)],'distinct');%按照2*2的网格划分图像
    Gray_re22_M  = mean(Gray_re22);
    Gray22_1 =Gray_re22_M(1).*ones((h/2),(l/2));
    Gray22_2 =Gray_re22_M(2).*ones((h/2),(l/2));
    Gray22_3 =Gray_re22_M(3).*ones((h/2),(l/2));
    Gray22_4 =Gray_re22_M(4).*ones((h/2),(l/2));
    th22_Gray(:,:,j) = [Gray22_1,Gray22_3;
                 Gray22_2,Gray22_4];  
      %4*4 local mean       
    Gray_Ext = padarray(Diff{j},[radius radius],'symmetric','both');        
    Gray_re44 = im2col(Gray_Ext,[(a/4) (b/4)],'distinct');
    Gray_re44_M  = mean(Gray_re44);
    Gray44_1 =Gray_re44_M(1).*ones((a/4),(b/4));
    Gray44_2 =Gray_re44_M(2).*ones((a/4),(b/4));
    Gray44_3 =Gray_re44_M(3).*ones((a/4),(b/4));
    Gray44_4 =Gray_re44_M(4).*ones((a/4),(b/4));
    Gray44_5 =Gray_re44_M(5).*ones((a/4),(b/4));
    Gray44_6 =Gray_re44_M(6).*ones((a/4),(b/4));
    Gray44_7 =Gray_re44_M(7).*ones((a/4),(b/4));
    Gray44_8 =Gray_re44_M(8).*ones((a/4),(b/4));
    Gray44_9 =Gray_re44_M(9).*ones((a/4),(b/4));
    Gray44_10 =Gray_re44_M(10).*ones((a/4),(b/4));
    Gray44_11 =Gray_re44_M(11).*ones((a/4),(b/4));
    Gray44_12 =Gray_re44_M(12).*ones((a/4),(b/4));
    Gray44_13 =Gray_re44_M(13).*ones((a/4),(b/4));
    Gray44_14 =Gray_re44_M(14).*ones((a/4),(b/4));
    Gray44_15 =Gray_re44_M(15).*ones((a/4),(b/4));
    Gray44_16 =Gray_re44_M(16).*ones((a/4),(b/4));
    th44_Gray_o(:,:,j) = [Gray44_1,Gray44_5,Gray44_9,Gray44_13;
                        Gray44_2,Gray44_6,Gray44_10,Gray44_14;
                        Gray44_3,Gray44_7,Gray44_11,Gray44_15;
                        Gray44_4,Gray44_8,Gray44_12,Gray44_16];  
    th44_Gray(:,:,j) = th44_Gray_o((1+radius):(a-radius),(1+radius):(b-radius),j); 
  
%     
     MeanDiff(j) = mean(mean(Diff{j}));
    over_Gray(:,:,j) = MeanDiff(j).*ones((a-2*radius),(b-2*radius));
   
  end  
end
% mean_M_neighbors = th44_Gray;
 mean_M_neighbors = (th22_Gray+th44_Gray+over_Gray)./3;
sum = zeros((a-2*radius),(b-2*radius));
   for k= 1:neighbors
    sum = sum + mean_M_neighbors(:,:,k);
   end 
   sum_M = sum./neighbors;
 
 Outex10_MG28th3(:,:,i) = sum_M;
 
 % Outex10_ME18th3(:,:,i) = sum_M;

 % Outex10_M2Ga28th3(:,:,i) = sum_M;

  % Outex10_M28th3(:,:,i) = sum_M;
end

 save('Outex10_MG28th3.mat','Outex10_MG28th3');


 % save('Outex10_ME18th3.mat','Outex10_ME18th3');


 % save('Outex10_M2Ga28th3.mat','Outex10_M2Ga28th3');

% save('Outex10_M28th3.mat','Outex10_M28th3');
