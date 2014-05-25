%This is a example to use the I-KSVD!
%use "error" mode
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
%% set parameters %%
params.blocksize       =      8;  
params.ReIncAtom       =      10;        %%£¨not initial dictionary£©the number of atoms of every image when training the resrt
params.IncAtom         =      20;        %% This is the number of the atoms of every image when training the initial dictionary.

params.codemode        =     'error';    %% 'error' or 'sparsity' two modes
params.sigma           =      20;        %% the parameter of 'error' mode. The error of every pixel.
params.control_sparse  =      4.8;       %% the parameter of 'error' mode. 
                                         %This is the threshold of controlling sparsity in 'error' mode. 
params.T               =      0.0157;    %%  %%The parameter of 'sparsity' mode. e.g. 0.0157 equals one non-zero coefficient in every column of a sparsity matrix.  
params.control_rmse    =      8;         %% The parameter of 'sparsity'mode. This is the threshold of controlling RMSE in 'sparsity' mode.

%% train dictionary  %%
Dicttrain(params,im_path);

%% set parameters of reconstructing%%
params.codemode    =    'error';  %%%%%%%%%  'error' or 'sparsity'%%%%%%%%%
params.sigma       =     20;
%% reconstruct image  %%
reconstruction(params,im_path,path_two,path_three,path_four);
%% show results  %%
IKSVD_path=strcat(current_path,'\result\IKSVD_sparse_1.mat');
ODL_path=strcat(current_path,'\result\ODL_20_sparse_1.mat');
RLS_path=strcat(current_path,'\result\RLS_20_sparse_1.mat');

load(IKSVD_path);
x=(1:length(imglist));
y_ksvd=sparse;
load(ODL_path);
y_odl=sparse;
load(RLS_path);
y_rls=sparse;

figure(1);
plot(x,y_ksvd,'g--+',x,y_odl,'b--+',x,y_rls,'r--+');
legend('I-KSVD','ODL','RLS');
axis on;
xlabel('image');
ylabel('number of nonzero coefficients');
title('Comparing the sparsity of decompostion')
grid on;

IKSVD_path=strcat(current_path,'\result\D.mat');
ODL_path=strcat(current_path,'\result\D_ODL_psnr.mat');
RLS_path=strcat(current_path,'\result\D_RLS_psnr.mat');

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
