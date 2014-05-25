function y=reconstruction(params,path,path_two,path_three,path_four)
%============================================================
%               demo2 - reconstruct an image
% This demo reads an image and reconstruct the image using 
% K-SVD dictionarie. The function displays the original and
% reconstructed image£¬and shows the resulting trained dictionary. 
%============================================================
tic

%*********************************
  
dirname = fullfile(path,'*.tif');
imglist = dir(dirname);

%% set parameters %%
load_path=path;
load_path(end-5:end)=[];
load_path=strcat(load_path,'\result\D');
load (load_path);

sparse=[];
 PSNR=[];
 RMSE=[];
BLOCK              =     2050;
params.blocksize   =     8;
params.dict        =     D;
%params.codemode    =    'error';  %%%%%%%%%  'error' or 'sparsity'%%%%%%%%%
%params.sigma       =     20;
% params.T           =     1;
params.maxval      =     255;
params.memusage    =    'low';
%% test  %%
 for k=1:length(imglist)

         disp( sprintf('process images %d',k));
         imgname=fullfile(path, imglist(k).name); 
         [num1,num2]=imreadImage(imgname,BLOCK,path_two);
         squaresum=0;
         spa=0;

      for i = 1 : num1
         for   j = 1 : num2
               m=(i-1)*num2+j;
                disp( sprintf('process subimages %d-%d-%d',i,j,m));
               [IMinblock] = getblock(m,path_three);
               params.x=IMinblock;
               [Ioutblock,nz]=ompreconstruct_mode(params,5);
                squaresum =squaresum+sum(sum((IMinblock-Ioutblock).^2));
                spa=spa+nz; 
           

         end
      end   
             [Psnr,Rmse]=psnr(squaresum,params.maxval,path_four); 
             RMSE=[RMSE,Rmse];       
             PSNR=[PSNR,Psnr];         
             sparse=[sparse,spa];
   path_copy=path;
   path_copy(end-5:end)=[];
   s_path_one=strcat(path_copy,'\result\IKSVD_psnr_1.mat');
   s_path_two=strcat(path_copy,'\result\IKSVD_sparse_1.mat');
   save (s_path_one,'PSNR');
   save (s_path_two,'sparse');
 end

  toc
end
  
  
  
