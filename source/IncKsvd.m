function [D,Gamma,err,gerr] = IncKsvd(params,varargin)

global MEM_LOW MEM_NORMAL MEM_HIGH memusage
global ompfunc ompparams exactsvd

CODE_SPARSITY = 1;
CODE_ERROR = 2;

MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;


%%%%% parse input parameters %%%%%


data = params.data;       %
IncAtom = params.IncAtom;
mu_coherence=params.mu_coherence;

ompparams = {'checkdict','off'};

% coding mode %

if (isfield(params,'codemode'))
  switch lower(params.codemode)
    case 'sparsity'
      codemode = CODE_SPARSITY;
      thresh = params.Tdata;
    case 'error'
      codemode = CODE_ERROR;
      thresh = params.Edata;
    otherwise
      error('Invalid coding mode specified');
  end
elseif (isfield(params,'Tdata'))
  codemode = CODE_SPARSITY;
  thresh = params.Tdata;
elseif (isfield(params,'Edata'))
  codemode = CODE_ERROR;
  thresh = params.Edata;

else
  error('Data sparse-coding target not specified');
end


% max number of atoms %

if (codemode==CODE_ERROR && isfield(params,'maxatoms'))
  ompparams{end+1} = 'maxatoms';
  ompparams{end+1} = params.maxatoms;
end


% memory usage %

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


% iteration count %

if (isfield(params,'ksvd_iternum'))
  ksvd_iternum = params.ksvd_iternum;
else
  ksvd_iternum = 10;
end

if (isfield(params,'IPR_iternum'))
 IPR_iternum=params.IPR_iternum;;
else
  IPR_iternum = 10;
end


% omp function %

if (codemode == CODE_SPARSITY)
  ompfunc = @omp;
  %min  |X - D*GAMMA|_2     s.t.  |GAMMA|_0 <= T
else
  ompfunc = @omp2;
  %min  |GAMMA|_0           s.t.  |X - D*GAMMA|_2 <= EPSILON
end


% status messages %

printiter = 0;
printreplaced = 0;
printerr = 0;
printgerr = 0;

verbose = 't';
msgdelta = -1;

for i = 1:length(varargin)
  if (ischar(varargin{i}))
    verbose = varargin{i};
  elseif (isnumeric(varargin{i}))
    msgdelta = varargin{i};
  else
    error('Invalid call syntax');
  end
end

for i = 1:length(verbose)
  switch lower(verbose(i))
    case 'i'                  %    i - iteration number
      printiter = 1;
    case 'r'                  %    r - number of replaced atoms
      printiter = 1;
      printreplaced = 1;
    case 't'                  %    t - target function value (and its value on the test data if provided)
      printiter = 1;
      printerr = 1;
      if (isfield(params,'testdata'))
        printgerr = 1;
      end
  end
end

if (msgdelta<=0 || isempty(verbose))
  msgdelta = -1; 
end

ompparams{end+1} = 'messages';
ompparams{end+1} = msgdelta;



% compute error flag %

comperr = (nargout>=3 || printerr);


% validation flag %

testgen = 0;
if (isfield(params,'testdata'))
  testdata = params.testdata;
  if (nargout>=4 || printgerr)
    testgen = 1;
  end
end


% data norms %

XtX = []; XtXg = [];
if (codemode==CODE_ERROR && memusage==MEM_HIGH)
  XtX = colnorms_squared(data);   %
  if (testgen)
    XtXg = colnorms_squared(testdata);
  end
end


% mutual incoherence limit %

if (isfield(params,'muthresh'))
  muthresh = params.muthresh;
else
  muthresh = 0.99;
end
if (muthresh < 0)
  error('invalid muthresh value, must be non-negative');
end


% exact svd computation %

exactsvd = 0;
if (isfield(params,'exact') && params.exact~=0)
  exactsvd = 1;
end
 

% determine dictionary size %
if (isfield(params,'dictsize'))    % this superceedes the size determined by initdict
  dictsize = params.dictsize;
end

if (isfield(params,'initdict'))
    dictsize = size(params.initdict,2);
end


% initialize the dictionary %

if (isempty(params.initdict)==1)
  data_ids = find(colnorms_squared(data) > 1e-6);   % ensure no zero data elements are chosen
  perm = randperm(length(data_ids));
  D = data(:,data_ids(perm(1:IncAtom)));
else
    D=params.initdict;
    G = [];
  if (memusage >= MEM_NORMAL)
    G = D'*D;
  end
  
 Gamma = sparsecode(data,D,XtX,G,thresh); %稀锟斤拷锟斤拷
 
 %%%%%%%%  set the initial value of the new atoms and add to the dictionary %%%%%%%%
 DicLocation=NotSarse(Gamma,dictsize,IncAtom);
  if ~isempty(DicLocation)
   AddD=data(:,DicLocation);
   D=[D,AddD];
  end
end
dictsize = size(D,2);
% normalize the dictionary %

D = normcols(D);   % normalize the dictionary

err = zeros(1,ksvd_iternum);
gerr = zeros(1,ksvd_iternum);

 errstr1 = 'RMSE';
 errstr2 = 'mean atomnum';



%%%%%%%%%%%%%%%%%  main loop  %%%%%%%%%%%%%%%%%


for iter = 1:ksvd_iternum

  G = [];
  if (memusage >= MEM_NORMAL)
    G = D'*D;
  end
  
  
  %%%%%  sparse coding  %%%%%
  
  Gamma = sparsecode(data,D,XtX,G,thresh);
  
  
  %%%%%  dictionary update  %%%%%
  
  replaced_atoms = zeros(1,dictsize);  % mark each atom replaced by optimize_atom
  
  unused_sigs = 1:size(data,2);  % tracks the signals that were used to replace "dead" atoms.
                                 % makes sure the same signal is not selected twice
  
  p = randperm(IncAtom)+(dictsize-IncAtom);
  tid = timerinit('updating atoms', dictsize);
  for j = 1:IncAtom
    [D(:,p(j)),gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom(data,D,p(j),Gamma,unused_sigs,replaced_atoms);
    Gamma(p(j),data_indices) = gamma_j;
    if (msgdelta>0)
      timereta(tid, j, msgdelta);
    end
  end
  if (msgdelta>0)
    printf('updating atoms: iteration %d/%d', dictsize, dictsize);
  end
  
  
  %%%%%  compute error  %%%%%
  
 if (comperr)
    [err1(iter),err2(iter)] = compute_err(D,Gamma,data);
  end
  if (testgen)
    if (memusage >= MEM_NORMAL)
      G = D'*D;
    end
    GammaG = sparsecode(testdata,D,XtXg,G,thresh);
    gerr(iter) = compute_err(D,GammaG,testdata);
  end
   
  
  %%%%%  clear dictionary  %%%%%
  
  [D,cleared_atoms] = cleardict(D,Gamma,data,muthresh,unused_sigs,replaced_atoms,IncAtom);
  
  
  %%%%%  print info  %%%%%
  
  info = sprintf('Iteration %d / %d complete', iter, ksvd_iternum);
  if (printerr)
    info = sprintf('%s, %s = %.4g', info, errstr1, err1(iter));
    info = sprintf('%s, %s = %.4g', info, errstr2, err2(iter));
  end
  if (printgerr)
    info = sprintf('%s, test %s = %.4g', info, errstr, gerr(iter));
  end
  if (printreplaced)
    info = sprintf('%s, replaced %d atoms', info, sum(replaced_atoms) + cleared_atoms);
  end
  
  if (printiter)
    disp(info);
    if (msgdelta>0), disp(' '); end
  end
  
end
%% Dictionary Decorrelation %%

Dictionary=D(:,1:end-IncAtom);
[n,m]=size(D);
if (m>n)
for ii=1:IPR_iternum
G=D'*D;
%%%%%%%% get the structural constraint set K %%%%%%%%%%%%%
K=G-eye(m);
KK=K(:,(end-IncAtom+1):end);
KK=abs(KK);
if (isempty(find(KK>=mu_coherence))==0)
K(find(K>mu_coherence))=mu_coherence;
K(find(K<-mu_coherence))=-mu_coherence;
else
     return
end
K=K+eye(m);

%%%%%%%%% get the spectral constraint set F  %%%%%%%%%%%%%%%%%
[Q,A]       = eig(K);
eigenvalues = flipud(sort(diag(A)));
eig_thresh      = eigenvalues(n);
A(find(A<eig_thresh))=0;
D=(A^0.5)*Q';
correct_part=D((end-n+1):end,(end-IncAtom+1):end);
D=[Dictionary,correct_part];

%%%%%%%%% Dictionary Rotation  %%%%%%%%%%%%
G = [];
  if (memusage >= MEM_NORMAL)
    G = D'*D;
  end
  D = normcols(D);
  D=double(D);
 Gamma = sparsecode(data,D,XtX,G,thresh);

 %C=data*(D*Gamma)';
C=D*Gamma*data'; 
[u,d,v]=svd(C);
 W=v*u';
 D=W*D;
correct_part=D(:,(end-IncAtom+1):end);
D=[Dictionary,correct_part];

end  
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            optimize_atom             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [atom,gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom(X,D,j,Gamma,unused_sigs,replaced_atoms)

global exactsvd

% data samples which use the atom, and the corresponding nonzero
% coefficients in Gamma
[gamma_j, data_indices] = sprow(Gamma, j);  %%获得稀疏系数矩阵中第j行的非零元素及其位置

if (length(data_indices) < 1)
  maxsignals = 5000;
  perm = randperm(length(unused_sigs));
  perm = perm(1:min(maxsignals,end));
  E = sum((X(:,unused_sigs(perm)) - D*Gamma(:,unused_sigs(perm))).^2);
  [d,i] = max(E);
  atom = X(:,unused_sigs(perm(i)));
  atom = atom./norm(atom);
  gamma_j = zeros(size(gamma_j));
  unused_sigs = unused_sigs([1:perm(i)-1,perm(i)+1:end]);
  replaced_atoms(j) = 1;
  return;
end

smallGamma = Gamma(:,data_indices);
Dj = D(:,j);

if (exactsvd)

  [atom,s,gamma_j] = svds(X(:,data_indices) - D*smallGamma + Dj*gamma_j, 1);
  gamma_j = s*gamma_j;
  
else
  
  atom = collincomb(X,data_indices,gamma_j') - D*(smallGamma*gamma_j') + Dj*(gamma_j*gamma_j');
  atom = atom/norm(atom);
  gamma_j = rowlincomb(atom,X,1:size(X,1),data_indices) - (atom'*D)*smallGamma + (atom'*Dj)*gamma_j;

end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             sparsecode               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Gamma = sparsecode(data,D,XtX,G,thresh)

global CODE_SPARSITY codemode
global MEM_HIGH memusage
global ompfunc ompparams

if (memusage < MEM_HIGH)
  Gamma = ompfunc(D,data,G,thresh,ompparams{:});
  
else  % memusage is high
  
  if (codemode == CODE_SPARSITY)
    Gamma = ompfunc(D'*data,G,thresh,ompparams{:});
    
  else
    Gamma = ompfunc(D'*data,XtX,G,thresh,ompparams{:});
  end
  
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             compute_err              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err1,err2] = compute_err(D,Gamma,data)
 
  mid =reperror2(data,D,Gamma);
  err1 = sqrt(sum(mid)/numel(data));

  err2 = nnz(Gamma)/size(data,2);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           cleardict                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [D,cleared_atoms] = cleardict(D,Gamma,X,muthresh,unused_sigs,replaced_atoms,IncAtom)

use_thresh = 4;  % at least this number of samples must use the atom to be kept 

    dictsize = size(D,2);

% compute error in blocks to conserve memory
err = zeros(1,size(X,2));
blocks = [1:3000:size(X,2) size(X,2)+1];
for i = 1:length(blocks)-1
  err(blocks(i):blocks(i+1)-1) = sum((X(:,blocks(i):blocks(i+1)-1)-D*Gamma(:,blocks(i):blocks(i+1)-1)).^2);
  
end

cleared_atoms = 0;
usecount = sum(abs(Gamma)>1e-7, 2);

for j = (dictsize-IncAtom+1):dictsize
  
  % compute G(:,j)
  Gj = D'*D(:,j);
  Gj(j) = 0;
  
  % replace atom
  if ( (max(Gj.^2)>muthresh^2 || usecount(j)<use_thresh) && ~replaced_atoms(j) )
    [y,i] = max(err(unused_sigs));
    D(:,j) = X(:,unused_sigs(i)) / norm(X(:,unused_sigs(i)));
    unused_sigs = unused_sigs([1:i-1,i+1:end]);
    cleared_atoms = cleared_atoms+1;
  end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  find location of the initial value of the new atoms   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DicLocation]=NotSarse(A,dictsize,IncAtom)
d=sum(A ~=0,1);        %
d=full(d);
[not_sparse location]=sort(d,'descend'); 
if (size(location,2)> dictsize)
NoSparseLocation=location(1:dictsize); %
% NoSparseA=full(A(:,NoSparseLocation));
NoSparseA=A(:,NoSparseLocation);
MinCoeff=min(min(NoSparseA)); 
NoSparseA(NoSparseA~=0)=1-MinCoeff+NoSparseA(NoSparseA~=0);
NormA=bsxfun(@rdivide,NoSparseA,sum(NoSparseA));   %
NormA=full(NormA);
NormA=NormA.*log(NormA);              
NormA(find(isnan(NormA)==1))=0;
NormA=-sum(NormA);                       %
% NormA=full(NormA);
[H DicLocation]=sort(NormA,'descend');
DicLocation=NoSparseLocation(DicLocation(1:IncAtom));
% AddD=A(:,DicLocation);
else
    DicLocation=[];
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            misc functions            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err2 = reperror2(X,D,Gamma)

% compute in blocks to conserve memory
err2 = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
  blockids = i : min(i+blocksize-1,size(X,2));
  err2(blockids) = sum((X(:,blockids) - D*Gamma(:,blockids)).^2);
end

end


function Y = colnorms_squared(X)

% compute in blocks to conserve memory
Y = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
  blockids = i : min(i+blocksize-1,size(X,2));
  Y(blockids) = sum(X(:,blockids).^2);
end

end


