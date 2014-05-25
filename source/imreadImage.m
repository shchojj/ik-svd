
function [num1,num2,pixel1,pixel2]=imreadImage(imgname,BLOCK,path)

IMin = imread(imgname);
IMin= double(IMin);
[NN1,NN2] = size(IMin);
pixel1=NN1*NN2;
pixel2=nnz(IMin);
num1=ceil(NN1/BLOCK); num2 =ceil(NN2/BLOCK);
num=max(num1,num2);
t1 = (0:num-1)*BLOCK + 1; t2 = (1:num)*BLOCK;
t3 = (0:num-1)*BLOCK+ 1; t4 = (1:num)*BLOCK;
IMinblock=cell(1,num1*num2);

% figure;
for i = 1 : num1
    for j = 1 : num2
         m=(i-1)*num2+j;
        IMinblock{m}=IMin(t1(i):min(t2(i),NN1), t3(j):min(t4(j),NN2));

    end
end
path_one=strcat(path,'IMinblock.mat');
path_two=strcat(path,'IMin.mat');
save (path_one,'IMinblock');
save (path_two,'IMin');
end
