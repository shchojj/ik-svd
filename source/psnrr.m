function [PSNR]= psnrr (squaresum,maxval)

% IMin=uint8(IMin);
%  figure; imshow(IMin); title('Original image');
 %PSNR=20*log10(maxval * sqrt(numel(IMin)) / norm(IMin(:)-Iout(:)));
PSNR=20*log10(maxval * sqrt(numel(IMin)) / sqrt(squaresum));

end