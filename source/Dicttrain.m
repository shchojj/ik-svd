function y=Dicttrain(params,path)
%============================================================
%              Dicttrain - Get a dictionaryset
% This demo reads an image and get the image using 
% IK-SVD dictionarie. The function displays the original and   
% reconstructed images and shows the resulting trained dictionary. 
%============================================================
clc
tic

%  
%********************************
%
dirname = fullfile(path,'*.TIF');
imglist = dir(dirname);

disp('  Available test images:');
disp(' ');
for k = 1:length(imglist)
  printf('  %d. %s', k, imglist(k).name);
end                                                                                                          
disp(' ');
imnum = 0;

 while (~isnumeric(imnum) || ~iswhole(imnum) || imnum<1 || imnum>length(imglist))
  imnum = input(sprintf(' How many base images selected   (%d-%d): ', 2, length(imglist)), 's');
  imnum = sscanf(imnum, '%d');
end
%% set parameters %%

params.blocksize       =      8;  
params.flag            =      0;  
params.ReIncAtom       =      10;   
params.IncAtom         =      20;   
params.maxval          =      255;   
params.trainnum        =      60000;  

%params.codemode        =     'error';   
%params.sigma           =      20;      
%params.control_sparse  =      4.8;   
                                       
                                     
                                                            
%params.T               =      0.0157;     
%params.control_rmse    =      8;     
                                      


params.ksvd_iternum    =       10;   %%  The iterations of KSVD
params.IPR_iternum     =       10;   %%  The iterations of IPR    
params.mu_coherence    =       0.9;  %% Controlling the degree of association of IPR.

params.exact           =       0;       %% exact K-SVD update instead of approximate (0)
params.memusage        =      'low';

D=[];    

    
%% select and read base image %%

 BaseImgNum = imnum;
 ImgNum  = length(imglist);
 SelectImgNum=randperm(ImgNum);
 SelectImgNum=SelectImgNum(1:BaseImgNum);
%MedVar = ceil(ImgNum/BaseImgNum);
%ImgNum  =1:ImgNum;
%SelectImgNum=find(mod(ImgNum,MedVar)==1);
%SelectImgNum  =  SelectImgNum(1:BaseImgNum);

disp('Performing K-SVD get a base dictionary...');

for k = 1:BaseImgNum
   i=SelectImgNum(k);
   
         disp( sprintf('train images %d', i));
         imgname=fullfile(path, imglist(i).name);   
         IMin = imread(imgname);
         IMin = double( IMin);
         if (length(size(IMin))>2)
            IMin = rgb2gray(IMin);
         end
        if (max(IMin(:))<2)
            IMin = IMin*255;
        end
        params.imagname=i;
        D = ksvdtrain(IMin,params,D);
       
end


%% Train restricted dictionary %%
params.flag=1;
disp('Performing train restricted dictionary....');

 for i=1:length(imglist)
      JudgeMat=find(SelectImgNum==i);
      if (isempty(JudgeMat)==0)
         continue
      else 
          
         disp( sprintf('train images %d', i));
          
         imgname=fullfile(path, imglist(i).name);   
         IMin = imread(imgname);
         IMin = double( IMin);
         if (length(size(IMin))>2)
        IMin = rgb2gray(IMin);
         end
        if (max(IMin(:))<2)
        IMin = IMin*255;
        end
         params.imagname=i;
                  D = ksvdtrain(IMin,params,D);   
     
      end
         
 end
path(end-5:end)=[];
s_path=strcat(path,'\result\D');
save (s_path,'D');
 toc
end
 