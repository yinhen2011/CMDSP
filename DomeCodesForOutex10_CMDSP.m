clear all;
clc;


% extract Outex_TC_00010 to the "rootpic" folder
% Load the dataset path, such as Outex_TC_00020, UMD, and UIUC
rootpic = 'E:\TextureClassification\outex\Outex_TC_00010\';
picNum = 4320;


% % Radius and Neighborhood
% R=2; P=8;
% load patternMappingriu2_1_8;
R=3;
P=16;
load patternMappingriu2_3_16;
% R=3;
% P=24;
% load patternMappingriu2_3_24;

% Load the multi-scale threshold
load Outex10_Cth3;
load Outex10_Gth3;
load Outex10_2Gath3;

% 
% load Outex10_M18th3;
% load Outex10_MG18th3;
% load Outex10_M2Ga18th3;

% 
% load Outex10_M28th3;
% load Outex10_MG28th3;
% load Outex10_M2Ga28th3;

% % 
load Outex10_M316th3;
load Outex10_MG316th3;
load Outex10_M2Ga316th3;

% 
% load Outex10_M216th3;
% load Outex10_MG216th3;
% load Outex10_M2Ga216th3;


% load Outex10_M324th3;
% load Outex10_MG324th3;
% load Outex10_M2Ga324th3;

% load Outex10_M524th3;
% load Outex10_MG524th3;
% load Outex10_M2Ga524th3;



for i=1:picNum;
    filename = sprintf('%s\\images\\%06d.ras', rootpic, i-1);
    fprintf('No.%d\n',i);
    Grayy = imread(filename);
    Grayx = im2double(Grayy);
    Gray = (Grayx-mean(Grayx(:)))/std(Grayx(:))*20+128; % image normalization, to remove global intensity

%  the corresponding multi-scale threshold of original image
       threC = Outex10_Cth3(:,:,i);
       % threM =Outex10_M28th3(:,:,i);
       % threM =Outex10_M18th3(:,:,i);
       threM =Outex10_M316th3(:,:,i);
       % threM =Outex10_M216th3(:,:,i);
       % threM =Outex10_M324th3(:,:,i);
        % threM =Outex10_M524th3(:,:,i);
       [LSP_I,LCP_I] = CMDSP(Gray,R,P,patternMappingriu2,'x',threC,threM);  
       [UM0] = FindUniform(LSP_I,P);


%  %  the corresponding multi-scale threshold of gradient image
       threG = Outex10_Gth3(:,:,i);
       % threMG =Outex10_MG28th3(:,:,i);
       % threMG =Outex10_MG18th3(:,:,i);
        % threMG =Outex10_MG216th3(:,:,i);
          threMG =Outex10_MG316th3(:,:,i);
       % % threMG =Outex10_MG524th3(:,:,i);
        % threMG =Outex10_MG324th3(:,:,i);
       Gradient_im2 = sobel2_grad0(Gray);
       Gradient_Gray = (Gradient_im2-mean(Gradient_im2(:)))/std(Gradient_im2(:))*20+128;
       [LSP_IG,LCP_G] = CMDSP(Gradient_Gray,R,P,patternMappingriu2,'x',threG,threMG);
       [UMG] = FindUniform(LSP_IG,P);


%the corresponding multi-scale threshold of Gaussian image
       threGa = Outex10_2Gath3(:,:,i);
       % threMGa =Outex10_M2Ga28th3(:,:,i);
       % threMGa =Outex10_M2Ga18th3(:,:,i);
       % threMGa =Outex10_M2Ga216th3(:,:,i);
       threMGa =Outex10_M2Ga316th3(:,:,i);
       % threMGa =Outex10_M2Ga324th3(:,:,i);
         % threMGa =Outex10_M2Ga524th3(:,:,i);

      image_scale = Gray;
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
   
      end
        GrayGa = Image_S(:,:,3);
        [LSP_Ga,LCP_Ga] = CMDSP(GrayGa,R,P,patternMappingriu2,'x',threGa,threMGa);

       % Determines if the pattern is uniform 
         [UMGa] = FindUniform(LSP_Ga,P);


         % multi-domain central pattern
     MDCP = 2^0*LCP_I + 2^1*LCP_Ga + 2^2*LCP_G;
     %  multi-domain shrinkage assessment pattern
     MDSAP = 2^0*UM0 + 2^1*UMG + 2^2*UMGa;
    % *********************************************************************      
    % Generate histogram of SAP_IOS
    LSP_IH(i,:) = hist(LSP_I(:),0:patternMappingriu2.num-1);
     % Generate histogram of SAP_IGrS
    LSP_GH(i,:) = hist(LSP_IG(:),0:patternMappingriu2.num-1);
     % Generate histogram of SAP_IGaS
    LSP_GaH(i,:) = hist(LSP_Ga(:),0:patternMappingriu2.num-1);
    
  % Generate histogram of multi-domain central pattern
     MDCPH= hist(MDCP(:),0:8-1);
     MDCPHH(i,:) = MDCPH;
 
 % Generate histogram of multi-domain shrinkage assessment pattern
     MDSAPH= hist(MDSAP(:),0:8-1);
     MDSAPHH(i,:) = MDSAPH;

     %The shrinkage pattern (SelectP) with the most unified patterns is selected
    % OSP in the paper
       CM0 = find(LSP_I==P+1);
       [xM0,yM0]=size(CM0);
       CMG = find(LSP_IG==P+1);
       [xMG,yMG]=size(CMG);
       CMGa = find(LSP_Ga==P+1);
       [xMGa,yMGa]=size(CMGa);          
       MinF = min(min(xM0,xMG),xMGa);
       if   CM0 == MinF
           SelectP = LSP_I;
       else
           if xMG == MinF
               SelectP = LSP_IG;
           else
               if xMGa == MinF
                   SelectP = LSP_Ga;
               end
           end
       end


% Generate histogram of OSP
 OSPH(i,:) = hist(SelectP(:),0:patternMappingriu2.num-1);


% The finally fused features Multi_M which is defined as CMDSP_OS/A/C in
% the paper; Generate histogram (CMDSPH) of Multi_M
    Multi_M = [SelectP(:),MDCP(:),MDSAP(:)];
    num1 = patternMappingriu2.num;
    num2 = 8;
    num3 = 8;
    Hist3D = zeros(num1,num2,num3);
    num_pix =size(LSP_I(:));
       for idx = 1:num_pix
          Hist3D(Multi_M(idx,1)+1,Multi_M(idx,2)+1,Multi_M(idx,3)+1) = Hist3D(Multi_M(idx,1)+1,Multi_M(idx,2)+1,Multi_M(idx,3)+1)+1;
       end
     Multi_MHH = reshape(Hist3D,1,numel(Hist3D));
     CMDSPH(i,:) = Multi_MHH;



end

% read picture ID of training and test samples, and read class ID of
% training and test samples
trainTxt = sprintf('%s000\\train.txt', rootpic);
testTxt = sprintf('%s000\\test.txt', rootpic);
[trainIDs, trainClassIDs] = ReadOutexTxt(trainTxt);
[testIDs, testClassIDs] = ReadOutexTxt(testTxt);
% 


% the classification result of MDCP
% we reported the result in table 2 in the paper 
trains = MDCPHH(trainIDs,:);
tests = MDCPHH(testIDs,:);
trainNum = size(trains,1);
testNum = size(tests,1);
DistMat = zeros(P,trainNum);
DM_MDCP = zeros(testNum,trainNum);
for i=1:testNum;
    test = tests(i,:);        
    DM_MDCP(i,:) = distMATChiSquare(trains,test)';
end
CP_MDCP=ClassifyOnNN(DM_MDCP,trainClassIDs,testClassIDs)


% the classification result of MDSAP
% we reported the result in table 2 in the paper 
trains = MDSAPHH(trainIDs,:);
tests = MDSAPHH(testIDs,:);
trainNum = size(trains,1);
testNum = size(tests,1);
DistMat = zeros(P,trainNum);
DM_MDSAP = zeros(testNum,trainNum);
for i=1:testNum;
    test = tests(i,:);        
    DM_MDSAP(i,:) = distMATChiSquare(trains,test)';
end
CP_MDSAP=ClassifyOnNN(DM_MDSAP,trainClassIDs,testClassIDs)



% the classification result of OSP
% we reported the result in table 2 in the paper 
trains = OSPH(trainIDs,:);
tests = OSPH(testIDs,:);
trainNum = size(trains,1);
testNum = size(tests,1);
DistMat = zeros(P,trainNum);
DM_OSP = zeros(testNum,trainNum);
for i=1:testNum;
    test = tests(i,:);        
    DM_OSP(i,:) = distMATChiSquare(trains,test)';
end
CP_OSP=ClassifyOnNN(DM_OSP,trainClassIDs,testClassIDs)



% the classification result of CMDSP
% we reported the result in table 2 and table 3 in the paper
trains = CMDSPH(trainIDs,:);
tests = CMDSPH(testIDs,:);
trainNum = size(trains,1);
testNum = size(tests,1);
DistMat = zeros(P,trainNum);
DM_CMDSP = zeros(testNum,trainNum);
for i=1:testNum;
    test = tests(i,:);        
    DM_CMDSP(i,:) = distMATChiSquare(trains,test)';
end
CP_CMDSP=ClassifyOnNN(DM_CMDSP,trainClassIDs,testClassIDs)

