function [y,nz] = ompreconstruct_mode(params,msgdelta)

% parse input arguments %
global MEM_LOW MEM_NORMAL MEM_HIGH memusage
global ompfunc ompparams exactsvd

CODE_SPARSITY = 1;
CODE_ERROR = 2;

x = params.x;
D = params.dict;
blocksize = params.blocksize;
stepsize = 8;

% blocksize %
if (numel(blocksize)==1)
  blocksize = ones(1,2)*blocksize;
end

ompparams = {'checkdict','off'};

%% coding mode %
if (isfield(params,'codemode'))
  switch lower(params.codemode)
    case 'sparsity'
      codemode = CODE_SPARSITY;
      thresh = params.T;
    case 'error'
      codemode = CODE_ERROR;
      thresh = sqrt(prod(blocksize)) * params.sigma;
    otherwise
      error('Invalid coding mode specified');
  end
end

% max number of atoms %
if (isfield(params,'maxatoms'))
  maxatoms = params.maxatoms;
else
  maxatoms = floor(prod(blocksize)/2);
end

if (codemode==CODE_ERROR && isfield(params,'maxatoms'))
  ompparams{end+1} = 'maxatoms';
  ompparams{end+1} = params.maxatoms;
end



% maxval %
if (isfield(params,'maxval'))
  maxval = params.maxval;
else
  maxval = 1;
end


% stepsize %

  if (numel(stepsize)==1)
    stepsize = ones(1,2)*stepsize;
  end



% msgdelta %
if (nargin <2)
  msgdelta = 5;
end
if (msgdelta<=0)
  msgdelta = -1;
end

% memory usage %
MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;

if (isfield(params,'memusage'))
  switch lower(params.memusage)
    case 'low'
      memusage = MEM_LOW;
    case 'normal'
      memusage = MEM_NORMAL;
    case 'high'
      memusage = MEM_HIGH;
    otherwise
      error('Invalid memory usage mode');
  end
else
  memusage = MEM_NORMAL;
end

% omp function %

if (codemode == CODE_SPARSITY)
  ompfunc = @omp;
  %min  |X - D*GAMMA|_2     s.t.  |GAMMA|_0 <= T
else
  ompfunc = @omp2;
  %min  |GAMMA|_0           s.t.  |X - D*GAMMA|_2 <= EPSILON
end


XtX = []; XtXg = [];
% compute G %

G = [];
if (memusage >= MEM_NORMAL)
  G = D'*D;
end


% verify dictionary normalization %

if (isempty(G))
  atomnorms = sum(D.*D);
else
  atomnorms = diag(G);
end
if (any(abs(atomnorms-1) > 1e-2))
  error('Dictionary columns must be normalized to unit length');
end




%%  reconsruct images %%%%%%%
nz                 =    0;         % count non-zeros in block representations
blocknum           =    prod(floor((size(x)-blocksize)./stepsize) + 1);
processedblocks    =    0;
tid                =    timerinit('ompdenoise', blocknum);

blocks             =    im2col(x,blocksize,'distinct');
[blocks, dc]       =    remove_dc(blocks,'columns');   % remove DC

  if (memusage == MEM_LOW)
      gamma = ompfunc(D,blocks,G,thresh,ompparams{:});
%     gamma = omp2(D,blocks,[],epsilon,'maxatoms',maxatoms,'checkdict','off');
%   else
%     gamma = omp2(D'*blocks,sum(blocks.*blocks),G,epsilon,'maxatoms',maxatoms,'checkdict','off');
%   end
  else
      if (codemode == CODE_SPARSITY)
       gamma = ompfunc(D'*blocks,G,thresh,ompparams{:});
    
      else
       gamma = ompfunc(D'*blocks,XtX,G,thresh,ompparams{:});
      end
  end
  
  
  nz              =     nz + nnz(gamma);   %Ï¡Êè¶È
 cleanblocks      =     add_dc(D*gamma, dc, 'columns');    
  y               =     col2im(cleanblocks,blocksize,size(x),'distinct');

if (msgdelta>0)
    processedblocks = processedblocks + size(blocks,2);
    timereta(tid, processedblocks, msgdelta);
end

if (msgdelta>0)
  timereta(tid, blocknum);
end




