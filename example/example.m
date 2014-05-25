%This is a example to use the I-KSVD!
%use "sparsity" mode
%============================================================
clear all;
clc;
current_path=pwd;
im_path=strcat(current_path,'\image');   %%data source
path_two=strcat(current_path,'\Mid-variables\');
path_three=strcat(current_path,'\Mid-variables\IMinblock.mat');
path_four=strcat(current_path,'\Mid-variables\IMin.mat');
dirname = fullfile(im_path,'*.TIF');
imglist = dir(dirname);
%% set parameters of training %%
params.blocksize       =       8;  
params.ReIncAtom       =      10;        %%£¨not training the initial dictionary£©the number of atoms of every image when training the restricted dictionary.
params.IncAtom         =      20;        %% This is the number of the atoms of every image when training the initial dictionary.

params.codemode        =     'sparsity';    %% 'error' or 'sparsity' two modes
params.sigma           =      10;        %% The parameter of 'error' mode. The error of every pixel.
params.control_sparse  =      4.8;       %% The parameter of 'error' mode. 
                                           %This is the threshold of controlling sparsity in 'error' mode. 
params.T               =      0.05024;    %%The parameter of 'sparsity' mode. e.g. 0.0157 equals one non-zero coefficient in every column of a sparsity matrix.  
params.control_rmse    =      5;         %% The parameter of 'sparsity'mode. This is the threshold of controlling RMSE in 'sparsity' mode. 
 
%% train dictionary  %%
Dicttrain(params,im_path);

%% set parameters of reconstructing%%
params.codemode    =    'sparsity';  %%%%%%%%%  'error' or 'sparsity'%%%%%%%%%
params.T           =     3.2;        %Notice! The value is different from before. 
                                     %e.g. 3.2 equals 3.2 non-zero coefficient in every column of a sparsity matrix.
%% reconstruct image  %%
reconstruction(params,im_path,path_two,path_three,path_four);

%% show results  %%
IKSVD_path=strcat(current_path,'\result\IKSVD_psnr_1.mat');
ODL_path=strcat(current_path,'\result\ODL_psnr_1.mat');
RLS_path=strcat(current_path,'\result\RLS_psnr_1.mat');

load(IKSVD_path);
x=(1:length(imglist));
y_ksvd=PSNR;
load(ODL_path);
y_odl=PSNR;
load(RLS_path);
y_rls=PSNR;

figure(1);
plot(x,y_ksvd,'g--+',x,y_odl,'b--+',x,y_rls,'r--+');
legend('I-KSVD','ODL','RLS');
axis on;
xlabel('image');
ylabel('PSNR(:dB)');
title('Comparing the PSNR')
grid on;

IKSVD_path=strcat(current_path,'\result\D.mat');
ODL_path=strcat(current_path,'\result\D_ODL_sparse.mat');
RLS_path=strcat(current_path,'\result\D_RLS_sparse.mat');

figure(2);
load(IKSVD_path);
subplot(1,3,1);
ImD=displayPatches(D);
imagesc(ImD); colormap('gray');
title('I-KSVD');
load(ODL_path);
subplot(1,3,2);
ImD=displayPatches(D);
imagesc(ImD); colormap('gray');
title('ODL');
subplot(1,3,3);
load(RLS_path);
ImD=displayPatches(D);
imagesc(ImD); colormap('gray');
title('RLS');