function [D] = ksvdtrain(x,params,dict,msgdelta)

%%%%% parse input parameters %%%%%

params.x  = x;
blocksize = params.blocksize;
params.initdict = dict;   
trainnum = params.trainnum;  

% blocksize %
if (numel(blocksize)==1)
  blocksize = ones(1,2)*blocksize;
end

% maxval %               
if (isfield(params,'maxval'))
  maxval = params.maxval;
else
  maxval =255;
  params.maxval = maxval;
end

% msgdelta %            
if (nargin<4)
  msgdelta = 5;
end

verbose = 't';
if (msgdelta <= 0)
  verbose='';
  msgdelta = -1;
end


if (isfield(params,'sigma'))
  sigma = params.sigma;
end
if (isfield(params,'T'))
  T = params.T;
end

params.Edata    = sqrt(prod(blocksize)) * sigma ;   % target error for omp
% params.reconstructmode = 'error';

params.Tdata    = T*params.blocksize*params.blocksize;   % arget sparsity for omp  *params.blocksize*params.blocksize
% params.reconstructmode = 'sparsity';

%%%% create training data %%%
ids = cell(2,1);
[ids{:}] = reggrid(size(x)-blocksize+1, trainnum, 'eqdist');  % 'eqdist'：按等间隔取训练样本的位置信息；'eqnum'：按每维等数量取训练样本的位置信息

params.data = sampgrid(x,blocksize,ids{:});  %按上面获取的位置获得训练样本数据

%%%保存训练样本信息%%%
m=params.imagname;
data=params.data;
s=['save trainsam_'  num2str(m)  '.mat ','data'];
eval(s);

% remove dc in blocks to conserve memory % 
blocksize = 2000;
for n = 1:blocksize:size(params.data,2)
  blockids = n : min(i+blocksize-1,size(params.data,2));
  params.data(:,blockids) = remove_dc(params.data(:,blockids),'columns');
end


%%%%% KSVD training %%%%%
if (params.flag==0)
    
  disp('IncKSVD training...');
    D = IncKsvd(params,verbose,msgdelta);    %进去非选择样本训练
else 
    D = ResIncKsvd1(params,verbose,msgdelta);  %进入选择性样本训练
end