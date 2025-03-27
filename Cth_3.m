% This code is used to generate the multiscale threshold corresponding to
% LCP

% It is from the thre_N,gc in the paper Enhanced Texture Classification through a Completed Multi-Domain Shrinkage Pattern by Bin Li*, Xiaochun Xu,* , and Q.M.Jonathan Wu,
% which is submitted to The Visual Computer

clear all;
clc;

% Add the path for the dataset such as Outex_TC_00020,UMD,UIUC
rootpic = 'E:\TextureClassification\outex\Outex_TC_00010\';
picNum = 4320;%the number of images

for i=1:picNum;
    filename = sprintf('%s\\images\\%06d.ras', rootpic, i-1);
    fprintf('No.%d\n',i);
    Gray = imread(filename); %read the image


 
    Grayy = im2double(Gray);
     Gray = (Grayy-mean(Grayy(:)))/std(Grayy(:))*20+128; 
    % 
    % The multiscale threshold used to generate the center pixel of the gradient image
    % Gradient_im2 = sobel2_grad0(Gray);
    % Gray = (Gradient_im2-mean(Gradient_im2(:)))/std(Gradient_im2(:))*20+128;


    % The multiscale threshold used to generate the center pixel of the entropy image
    % h = [1,1,1,1,1;
    %     1,1,1,1,1;
    %     1,1,1,1,1;
    %     1,1,1,1,1;
    %     1,1,1,1,1];   
    % Gray = entropyfilt(Gray,h);


% %   The multiscale threshold used to generate the center pixel of the Gaussian image
%     image_scale =  Gray;
%    imgExt = padarray(image_scale,[3 3],'symmetric','both');
%    sigma = 2^0.25;
%    scale = 4;
%    xx = 2*ceil(2*sigma)+1;
%    Image_S(:,:,1) = Gray;
%    for gaussianconv = 1:(scale-1)
%         h = fspecial('gaussian', 2*ceil(2*sigma)+1, sigma);%h=fspecial('gaussian',5,3);
%         Image_filter = imfilter(imgExt,h);
%         Image_S(:,:,gaussianconv+1) = Image_filter(4:end-3,4:end-3);
%         imgExt = Image_filter;
% %         figure(gaussianconv+1);
% %         imshow(imgExt,[]);     
%    end
%      Gray = Image_S(:,:,3);


    [a,b] = size(Gray); 
    
    %compute the mean
    Gray_re11 = im2col(Gray,[a ,b],'distinct');
    Gray_re11_M = mean(Gray_re11);
    th11_Gray = Gray_re11_M.*ones(a,b); 
    
    %2*2 local mean
    Gray_re22 = im2col(Gray,[(a/2) (b/2)],'distinct');
    Gray_re22_M  = mean(Gray_re22);
    Gray22_1 =Gray_re22_M(1).*ones((a/2),(b/2));
    Gray22_2 =Gray_re22_M(2).*ones((a/2),(b/2));
    Gray22_3 =Gray_re22_M(3).*ones((a/2),(b/2));
    Gray22_4 =Gray_re22_M(4).*ones((a/2),(b/2));
    th22_Gray = [Gray22_1,Gray22_3;
                 Gray22_2,Gray22_4];  
    
      %4*4 local mean
    Gray_re44 = im2col(Gray,[(a/4) (b/4)],'distinct');
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
    th44_Gray = [Gray44_1,Gray44_5,Gray44_9,Gray44_13;
                 Gray44_2,Gray44_6,Gray44_10,Gray44_14;
                 Gray44_3,Gray44_7,Gray44_11,Gray44_15;
                 Gray44_4,Gray44_8,Gray44_12,Gray44_16];
   
       % multi-scale hierarchical threshold
    Outex10_Gth3(:,:,i) = (th11_Gray + th22_Gray + th44_Gray )./3;     

    
    % Outex10_Cth3(:,:,i) = (th11_Gray + th22_Gray + th44_Gray )./3;

 
    % Outex10_Eth3(:,:,i) = (th11_Gray + th22_Gray + th44_Gray )./3;

   
    % Outex10_2Gath3(:,:,i) = (th11_Gray + th22_Gray + th44_Gray )./3;
end


% save('Outex10_Gth3.mat','Outex10_Gth3');

save('Outex10_Cth3.mat','Outex10_Cth3');

% save('Outex10_Eth3.mat','Outex10_Eth3');
% save('Outex10_2Gath3.mat','Outex10_2Gath3');
  
