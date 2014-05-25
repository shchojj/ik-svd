function [PSNR,RMSE]= psnr (squaresum,maxval,path)
 load(path,'IMin');
%load(path);
RMSE=sqrt(squaresum/numel(IMin));
PSNR=20*log10(maxval * sqrt(numel(IMin)) / sqrt(squaresum));
end