function [Iout,N,dictsize]=getblock(num,path)
load (path,'IMinblock')
Iout=IMinblock{num};
[NN1,NN2] = size(Iout);
N=NN1*NN2;

%if (nargin>1)
%nz =nnz(Iout);
%dictsize=round(nz*dictsize/(BLOCK*BLOCK));
%end
end

