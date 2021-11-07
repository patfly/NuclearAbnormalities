%% Code written by Patrick Flynn, Mitchison Lab, Harvard Medical School
% 20211107

% This code uses the output labeled images from the free NucleAIzer segemntation
% software and quantifies the number of micronuclei.

% Prior to submitting images to the NucleAIzer software, images were 
% pre-processed in ImageJ as follows: unsharpened mask (Radius: 8 pixels, Max Weight: 0.6) 
% followed by a FFT Bandpass filter (Large structures 60 pixels, Small structures 3 pixels, 
% autoscale after filtering). This fitlering was empirically determined to
% improve micronuclei detection.



%% Import file images prodcued by nucleAIzer
% Input both the filtered original image and the labeled image from
% NucleAIzer

%% Prepare Images
clear all
close all



[filename, pathname] = uigetfile('*.*', 'Select a file','MultiSelect','on'); 

%%
%  Ensure even number of images imported,  Original and Labeled
file_length = length(filename);
xx = mod(file_length,2);
if xx ~= 0
    disp('Wrong File Number');
    return;
end

%% Loop through imported images
counter= 0;
for i=2:2:length(filename)
    counter = counter+1;
    L = imread([pathname filename{i}]);
    
   

    origIm = imread([pathname filename{i-1}]);

    % Clear border

    stats = regionprops(L,'Centroid',...
                    'MajorAxisLength','MinorAxisLength','ConvexImage','Orientation','Eccentricity','Perimeter','Extent','ConvexArea');

    % Remove NaN
    idx = find([stats.MajorAxisLength] > 2 & [stats.MajorAxisLength] > 0);
    L_binary = uint16(ismember(L, idx));  

    L_clear = imclearborder(L_binary);

    L_final = L_clear.*L;

% Empirically determined size criteria for different cell lines and
% magnifications

% **************   20x:  <40   and >50,    10x:  <15 and >25  ****  

%  ********  For BJ Drug Tests:  ([stats.ConvexArea] < 700 & [stats.MajorAxisLength] < 40);

%  ********  For CANCERS Drug Tests:  ([stats.ConvexArea] < 500 & [stats.MajorAxisLength] < 40);

 stats = regionprops(L_final,'Centroid',...
                    'MajorAxisLength','MinorAxisLength','ConvexImage','Orientation','Eccentricity','Perimeter','Extent','ConvexArea');


    % FIRST PASS TO Identify micornculei
    idx = find([stats.ConvexArea] < 500 & [stats.MajorAxisLength] < 40 & [stats.ConvexArea] > 0);
    L_mn = uint16(ismember(L_final, idx));  
    L_final_mn = L_mn.*L_final;
     
%reset the L-final_mn to 1-xxx for combination later
    L_new = zeros(size(L_final_mn));
     [val_0,idx] = unique(L_final_mn) ;
     val = val_0(2:end);
     nodes_1 = L_final_mn;
     for i=1:length(val)
         nodes_1(L_final_mn==val(i))=i;
         L_new=nodes_1;
     end
%****************************************    
      
     
      if length(idx)>0
    
           II_mn=labeloverlay(origIm,L_final_mn);
%            figure
%            imshow(II_mn,[]);
%            title('First Catch MN');
%           imwrite(II_mn,'firstPASS_MN.tif');
      end
     
    cc = bwconncomp(L_mn,8);
    num_MN(counter)  = cc.NumObjects;
    
   
    
%  IDENTIFY LARGER OBJECTS AS PROXY FOR CELL #


     idx = find([stats.ConvexArea] > 600 & [stats.MajorAxisLength] > 15); 
    L_cell = uint16(ismember(L_final, idx));  
    L_final_cell = L_cell.*L_final;
    
     if length(idx)>0
      
         II=labeloverlay(origIm,L_final_cell);
%          figure
%           imshow(II,[]);
%           title('Cells Raw');
%           
%           imwrite(II,'cellsRAW.tif');
      end

    
    
    
% ******************************************
%      % With closure and Fill  % Only use for muiltilobed nuclei produced
%      by specific anti-mitotics
    L_bin = (L_final_cell>0);
     se = strel('disk',4);
     IB  = imclose(L_bin,se);
     IB2 = imfill(IB,'holes');
      LL_cell = uint16(bwlabel(IB2,4));


% %*************************************************

 
% *****************************************
     % Without
   %  LL_cell=uint16(L_final_cell);
%*****************************************************


 %%%   REMOVE POSSIBLE MN AT THIS POINT AFTER CLOSURE
        
        
        stats = regionprops(LL_cell,'Centroid',...
                    'MajorAxisLength','MinorAxisLength','ConvexImage','Orientation','Eccentricity','Perimeter','Extent','ConvexArea');
% BJ:  >850 and >15
% Cancers: >650 and >15

         idx = find([stats.ConvexArea] > 650 & [stats.MajorAxisLength] > 15);
        LL_cell_2 = uint16(ismember(LL_cell, idx));  
        cellLabel = bwlabel(LL_cell_2,4);
         if length(idx)>0    
             II=labeloverlay(origIm,cellLabel);
             figure
              imshow(II,[]);
               title('Cells Closed');
         end
        
    

% With Closure
   cc = bwconncomp(LL_cell_2,4);
   num_Cell_Closed(counter)  = cc.NumObjects;
    

 % WITH NO CLOSURE
    num_Cell_Orig(counter) = length(unique(L_final_cell));
    
    

Cell_save = label2rgb(LL_cell_2,'jet','w','shuffle');
 title('Measured Cells');
 imwrite(II,'cellsMeasured.tif');
 


% Find remaining MN and add to first MN count

% BJ:find([stats.ConvexArea] < 850 & [stats.MajorAxisLength] < 40& [stats.MajorAxisLength] > 2 & [stats.ConvexArea] > 10)

idx = find([stats.ConvexArea] < 850 & [stats.MajorAxisLength] < 10& [stats.MajorAxisLength] > 2 & [stats.ConvexArea] > 10);
LL_mn_2 = uint16(ismember(LL_cell, idx));  
L_final_mn_2 = LL_mn_2.*L_final;

 if length(idx)>0

     II=labeloverlay(origIm,L_final_mn_2);
%      figure
%       imshow(II,[]);
%       title('Second Catch MN');
%       imwrite(II,'secondPassMN.tif');
  end

    cc = bwconncomp(L_final_mn_2,8);
     num_MN_2(counter)  = cc.NumObjects;


       L_new_big = zeros(size(L_final_mn_2));
     [val_0,idx] = unique(L_final_mn_2) ;
     val = val_0(2:end);
     nodes_1 = L_final_mn_2;
     for i=1:length(val)
         nodes_1(L_final_mn_2==val(i))=i;
         L_new_big=nodes_1;
     end


%% Seprate detetcion for very small objects (20x)

    
    %  REMOVE ALL DETECTED OBJECTS, THEN FIND REMAINING MNi BASED ON SIZE
  
     maskObj = imbinarize(L,0);
     
     %Expand objects to get rid of edges
%      se=strel('disk',2);
%      maskObjDil = imdilate(maskObj,se);
    % maskHoles = imfill(maskObj,'holes');

     
     maskRemove = ~maskObj;
    
%     figure
%     imshow(origIm,[min(origIm(:)),0.5*max(origIm(:))]);
    
    threshVal = adaptthresh(origIm);
    threshIm = imbinarize(origIm,threshVal);
    
%     se = strel('disk',4);
    imRemoved = threshIm.*maskRemove;
%     imDil = imdilate(imRemoved,se);
%     imRemovedHole = imfill(imDil,'holes');
    
    
    imMN_final = imclearborder(imRemoved);
    
%      figure
%      imshow(imMN_final,[]);
     labelMN = bwlabel(imMN_final);
     
     
     % Watershed
%      D = bwdist(~imMN_final);
%      D = -D;
%      L = watershed(D,4);
%     L(~imMN_final) = 0;
%     rgb = label2rgb(L,'jet',[.5 .5 .5]);
% imshow(rgb)
% title('Watershed Transform')
% 
%    labelMN = L;
    
    
    
    
    % run regionprops on remaining objects
    
     stats = regionprops(labelMN,origIm,'Centroid',...
                    'MajorAxisLength','MinorAxisLength','ConvexImage','Orientation','Eccentricity','Perimeter','Extent','ConvexArea','MeanIntensity');
              % 1.5 for A549 because of shmutz
              % 1.35 for rest
      intensityThresh = 1.35*median(origIm(:));

   % Throw out objects that are not MN]
     idx = find([stats.MinorAxisLength] > 6 & [stats.MajorAxisLength] < 25 & [stats.Eccentricity] < 0.67 &  [stats.MeanIntensity] > intensityThresh);% & [stats.MeanIntensity]>median(origIm(:)) ); 
    imMN = uint16(ismember(labelMN, idx));  
  
  imMNL = bwlabel(imMN);
  
  if length(idx)>0
      
     II=labeloverlay(origIm,imMNL);
%      figure
%   imshow(II,[]);
%   title('Last Catch MN');
%     imwrite(II,'LASTPassMN.tif');
  end

  
  num_MN_3(counter) = length(idx);
  
  
  %combining mn label matrices
  %reset the new mn to 1xxx-yyy for combination later
    L_new_2 = zeros(size(imMNL));
     [val_0,idx] = unique(imMNL) ;
     val = val_0(2:end);
     nodes_1 = uint16(imMNL);
     for j=1:length(val)
         offset=j+length(unique(L_new));
         nodes_1(imMNL==val(j))=offset;
         L_new_2=nodes_1;
     end
     
     
     L_all_mn_1 = L_new;
     
     if num_MN_3>0
      L_all_mn_1 = bitor(L_all_mn_1,L_new_2);
     end
     
     if num_MN_2>0
      L_all_mn = bitor(L_all_mn_1,L_new_big);
     else
         L_all_mn = L_all_mn_1;
     end
      II=labeloverlay(origIm,L_all_mn,'Colormap','hsv');
     figure
    imshow(II,[]);
    title('All mn');
     imwrite(II,'AllPassMN.tif');
  
  
 end
  
% Sum up micornculei detected at multiple iterations
    num_MN_final = num_MN_2+num_MN+num_MN_3;
    num_MN_final
    
% Average cell number based on nuclei detected before and after closure
% call ***especially important for the heavily multi-lobed nuclei***
    num_cell =   mean([num_Cell_Closed;num_Cell_Orig],1);
    num_cell