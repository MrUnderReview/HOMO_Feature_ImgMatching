%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CVPR 2025 Submission Paper ID #16689
% This code is only for the purpose of reviewing

% HOMO-feature cross-modal image matching algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc;
clear;
addpath('functions','func_Math','func_HOMO')
save_path = '.\save_image\';

%% Program parameters
rot_flag = 1;  % Is there any obvious rotation difference, yes:1, no:0
               %   If the images do not have significant rotation (roughly less than 10°),
               %   we suggest setting rot_flag=0 to ensure fairness in comparison
               %   with non-rotation and network-based methods, and to save computing time.
scl_flag = 1;  % Is there any obvious scale difference
               %   If the images do not have significant scaling (roughly less than 1.2x),
               %   we suggest setting scl_flag=0 to save computing time.
par_flag = 1;  % Do you want parallel computing in multi-scale strategy
               %   (This may cause problem depending on computer system,
               %    for example, the RAM burden. From testing, 64GB is enough for 14-core)
pool_num = 8;  % You can set the max number of used CPU-core here.

trans_form = 'affine';  % What spatial transformation model do you need at the end:
                        % 'similarity', 'affine', 'projective'
out_form = 'Union';  % What output form of the image pair do you need at the end: 
                     % 'Reference', 'Union', 'Inter'
chg_scale = 1; % Do you want the resolution of sensed image be changed to the reference
Is_flag   = 1; % Do you want the visualization of registration results
I3_flag   = 1; % Do you want the Overlap form of registration results
I4_flag   = 1; % Do you want the Mosaic form of registration results

%% Algorithm parameters
nOctaves = 4;    % Gaussian pyramid octave number, default: 3 or 4
nLayers  = 3;    % Gaussian pyramid  layer number, default: 3 or 4
G_resize = 1.2;  % Gaussian pyramid downsampling ratio, default: 1.2 or 2
G_sigma  = 1.6;  % Gaussian blurring standard deviation, default: 1.6
key_type = 'PC-ShiTomasi'; % What kind of feature point do you want as the keypoint:
                           % 'Harris', 'ShiTomasi', 'PC-Harris', 'PC-ShiTomasi', 'Detector-free'
                           % 'Detector-free' means directly locating evenly distributed keypoints
                           %  from the images without considering image information,
                           %  it is to benchmark against methods such as LoFTR, XoFTR, and ASpanFormer.
radius   = 1;    % Local non-maximum suppression radius, default: 1 or 2
thresh   = 0;    % Keypoints response threshold, default: 0 or 50
Npoint   = 5000; % Keypoints number threshold, default: 5000 or 2000
patchsize= 72;   % GPolar patchsize, default: 72 or 96
NBA      = 12;   % GPolar spatial division, default: 12
NBO      = 12;   % GPolar angle division, default: 12

%% Image input and preprocessing
if ~exist('file1','var')
    file1 = [];
end
if ~exist('file2','var')
    file2 = [];
end
[image_1,file1,DataInfo1] = Readimage(file1);
[image_2,file2,DataInfo2] = Readimage(file2);
[image_1,resample1] = Deal_Extreme(image_1,64,512,1);
[image_2,resample2] = Deal_Extreme(image_2,64,512,1);
[I1_s,I1] = Preproscessing(image_1,1,[]); figure,imshow(I1_s); drawnow
[I2_s,I2] = Preproscessing(image_2,1,[]); figure,imshow(I2_s); drawnow

%%
if par_flag
    pool_num = min(pool_num, maxNumCompThreads);
    pool = gcp('nocreate');
    if ~isempty(pool)
        if pool_num ~= pool.NumWorkers
            delete(pool); pool = [];
        end
    end
    if isempty(pool)
        parpool(pool_num);  % Start parallel computing, time needed.
    end
end
warning off
    fprintf('\n** Image matching starts, have fun\n\n');



%% Build HOMO feature pyramids ---------- MOM and DoM
tic,[MOM_pyr1,DoMOM_pyr1] = Build_Homo_Pyramid(I1,...
    nOctaves,nLayers,G_resize,G_sigma,patchsize,NBA,1);
    t(1)=toc; fprintf(['Done: HOMO-feature extraction of reference image, time cost: ',num2str(t(1)),'s\n']);
%     Display_Pyramid(MOM_pyr1,'MOM Pyramid of Reference Image',0);
%     Display_Pyramid(DoMOM_pyr1,'DoMOM Pyramid of Reference Image',0);

tic,[MOM_pyr2,DoMOM_pyr2] = Build_Homo_Pyramid(I2,...
    nOctaves,nLayers,G_resize,G_sigma,patchsize,NBA,1);
    t(2)=toc; fprintf(['Done: HOMO-feature extraction of sensed image, time cost: ',num2str(t(2)),'s\n']);
%     Display_Pyramid(MOM_pyr2,'MOM Pyramid of Sensed Image',0);
%     Display_Pyramid(DoMOM_pyr2,'DoMOM Pyramid of Sensed Image',0);

%% Keypoints detection ---------- DoM enhancing
ratio = sqrt((size(I1,1)*size(I1,2))/(size(I2,1)*size(I2,2)));
if ratio>=1
    r2 = radius; r1 = round(radius*ratio);
else
    r1 = radius; r2 = round(radius/ratio);
end

tic,keypoints_1 = Detect_Homo_Keypoint(I1,DoMOM_pyr1,6,thresh,r1,Npoint,G_resize,key_type);
    t(3)=toc; fprintf(['Done: Keypoints detection of sensed image, time cost: ',num2str(t(3)),'s\n']);
    figure,imshow(I1_s); hold on; plot(keypoints_1(:,1),keypoints_1(:,2),'r.'); 
    title(['Reference image —— ',num2str(size(keypoints_1,1)),' keypoints']); drawnow
    clear DoMOM_pyr1

tic,keypoints_2 = Detect_Homo_Keypoint(I2,DoMOM_pyr2,6,thresh,r2,Npoint,G_resize,key_type);
    t(4)=toc; fprintf(['Done: Keypoints detection of sensed image, time cost: ',num2str(t(4)),'s\n']);
    figure,imshow(I2_s); hold on; plot(keypoints_2(:,1),keypoints_2(:,2),'r.');
    title(['Sensed image —— ',num2str(size(keypoints_2,1)),' keypoints']); drawnow
    clear DoMOM_pyr2
    
%% Keypoints description and matching (Multiscale strategy) ---------- GPolar and MsS
tic,[cor1,cor2] = Multiscale_Strategy(keypoints_1,keypoints_2,MOM_pyr1,MOM_pyr2,...
    patchsize,NBA,NBO,G_resize,5,1,rot_flag,scl_flag,par_flag);
    t(5)=toc; fprintf(['Done: Keypoints description and matching, time cost: ',num2str(t(5)),'s\n']);
    clear MOM_pyr1 MOM_pyr2
    matchment = Show_Matches(I1_s,I2_s,cor1,cor2,0);



%% Image transformation (Geography enable)
tic,[I1_r,I2_r,I1_rs,I2_rs,I3,I4,t_form,pos] = Transformation(image_1,image_2,...
    cor1,cor2,trans_form,out_form,chg_scale,Is_flag,I3_flag,I4_flag);
    t(6)=toc; fprintf(['Done: Image tranformation, time cost: ',num2str(t(6)),'s\n\n']);
    figure,imshow(I3); title('Overlap Form'); drawnow
    figure,imshow(I4); title('Mosaic Form'); drawnow

%%
T=num2str(sum(t)); fprintf([' Done image registration, total time: ',T,'s\n\n']);
    
%% Save results
if (exist(save_path,'dir')==0)  % If file folder does not exist
    mkdir(save_path);
end
Date = datestr(now,'yyyy-mm-dd_HH-MM-SS__'); tic
if exist('matchment','var') && ~isempty(matchment) && isvalid(matchment)
    saveas(matchment, [save_path,Date,'0 Matching Result','.jpg']);
end
if strcmpi(out_form,'reference')
    Imwrite(image_1, [save_path,Date,'1 Reference Image','.tif']);
    if Is_flag
        [I1_s,~] = Preproscessing(image_1,1,[]); 
        Imwrite(I1_s, [save_path,Date,'3 Reference Image Show','.png']);
    end
else
    Imwrite(I1_r , [save_path,Date,'1 Reference Image','.tif']);
    Imwrite(I1_rs, [save_path,Date,'3 Reference Image Show','.png']);
end
Imwrite(I2_r , [save_path,Date,'2 Registered Image','.tif']);
Imwrite(I2_rs, [save_path,Date,'4 Registered Image Show','.png']);
Imwrite(I3   , [save_path,Date,'5 Overlap of results','.png']);
Imwrite(I4   , [save_path,Date,'6 Mosaic of results','.png']);
T=num2str(toc); fprintf([' Registration results are saved in the save_image folder, time: ',T,'s\n']);