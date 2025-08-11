function[varargout]=nkp(M, d, varargin)

% NKP: Computes the nearest Kronecker product of a matrix M.
%
% Mh = NKP(M,d) computes the nearest Kronecker product approximation
% Mh = M1 ⨂ M2 ⨂ ... ⨂ Md of a matrix M.
%
% Mh = NKP(M,d,q) computes the nearest Kronecker rank q approximation.
% Default: Kronecker rank 1 approximation. Leave empty ([]) to use default.
%
% Mh = NKP(M,d,q,tol) attempts to compute an approximation of Kronecker rank
% smaller or equal to q satisfying a tolerance tol. If the tolerance cannot
% be met, the approximation has Kronecker rank q. Default: 1e-14.
% Leave empty ([]) to use default.
%
% Mh = NKP(M,d,q,tol,name,value) specifies name/value pairs for optional
% parameters.
%   'blocksize' - specifies the blocksize for the block
%                 partitioning of the matrix M, defined as
%                 blocksize = {[m1 n1], [m2 n2], ..., [md nd]}.
%                 If this parameter is not specified, the algorithm
%                 automatically tries to find it.
%   'algo'      - specifies the algorithm for the low rank approximation.
%                 Possible choices are 'svd' (default) and 'aca' for d = 2
%                 and 'cp_als' for d > 2.
%   'format'    - output format for the factor matrices. Available formats
%                 are 'tensor' and 'cell'. Default: 'cell'.
%   'singv'     - plots the singular values of the reordered matrix for
%                 d = 2. Default: true.
%
% [M1,M2,...,Md] = NKP(M,d,...) returns instead the factor matrices of
% the approximation, stored along the pages (or elements) of Mk.
%
% [Mh,M1,M2,...,Md] = NKP(M,d,...) returns the approximation as well as
% the factor matrices, stored along the pages (or elements) of Mk.
%
% References:
% [1] C. F. Van Loan and N. Pitsianis. Approximation with Kronecker products.
% In Linear Algebra for Large Scale and Real-Time Applications. Springer, 1993.
% [2] G. H. Golub and C. F. Van Loan. Matrix computations. JHU press, 2013.
% [3] A. N. Langville and W. J. Stewart. A Kronecker product approximate 
% preconditioner for SANs. Numerical Linear Algebra with Applications, 2004.

%% Set algorithm parameters

% Set default parameters
Default{1}=1;
Default{2}=1e-14;

% Replace empty inputs with default parameters
def=cell2mat(cellfun(@isempty, varargin, 'UniformOutput', false));
[varargin{def}]=Default{def};

Param = inputParser;
Param.addRequired('matrix',     @ismatrix);
Param.addRequired('dimension',  @(x) isscalar(x) & x > 0);
Param.addOptional('rank',       Default{1}, @(x) isscalar(x) & x > 0);
Param.addOptional('tolerance',  Default{2}, @(x) isscalar(x) & x > 0);
Param.addParameter('blocksize', 'auto', @iscell);
Param.addParameter('algo',      'svd', @(x) ismember(x,{'svd','aca','cp_als'}));
Param.addParameter('format',    'cell', @(x) ismember(x,{'tensor','cell'}));
Param.addParameter('singv',     true, @islogical);
Param.parse(M, d, varargin{:});

%% Retrieve parameters
q=Param.Results.rank;
tol=Param.Results.tolerance;
blocksize=Param.Results.blocksize;
algo=Param.Results.algo;
fout=Param.Results.format;
singv=Param.Results.singv;

%% Nearest Kronecker product approximation

if ~iscell(blocksize)
    f=kronfact(M, false);

    if isempty(f)
        error('Blocksize partitioning failed. Consider providing it explicitly.')
    else
        if length(f)>1
            warning('The sparsity pattern admits multiple factorizations. Choosing the first one.')
        end

        n=f{1}';
        blocksize=num2cell([n n], 2);
    end
end

if length(blocksize)~=d
    error('The factorization length (%d) does not match the prescribed length (%d).', length(blocksize), d)
end

s=prod(cat(1,blocksize{:}),2); % Size of R
R=Rop(M, blocksize); % Compute R(M)
U=cell(1,d);

if isempty(R.subs)
    I=cell(1,d);
else
    [Rs,I]=squash(R); % Remove zero slices

    if d==2 % Length 2 factorizations

        % Convert the sptensor Rs to a standard Matlab sparse matrix
        Rs=spmatrix(Rs);

        switch algo
            case 'svd' % Truncated SVD
                [U{1},S,U{2}]=svds(Rs, q);
                S=diag(S);

                ind=S>tol;
                S=S(ind);

                U{1}=sqrt(S)'.*U{1}(:,ind);
                U{2}=sqrt(S)'.*U{2}(:,ind);

                q=sum(ind);

            case 'aca' % Adaptive cross approximation
                [U{1},U{2},~,~,q]=aca(Rs, tol, q);
        end


        if singv
            figure
            sigma=svd(Rs);
            sR=size(Rs);
            l=min([sR 50]);
            semilogy(1:l, sigma(1:l), '.-b', 'MarkerSize', 15)
            grid on;
            title('Singular values of reordered matrix')
        end

    else % Factorizations of length > 2

        % CP decomposition
        [T]=cp_als(Rs, q, 'tol', tol, 'printitn', 0);

        U=tocell(T);
        q=ncomponents(T);
    end
end

% Factor matrices
F=cell(1,d);

switch fout
    case 'tensor'

        if issparse(M) && q==1
            for k=1:d
                F{k}=sparse(s(k), q);
                F{k}(I{k},:)=U{k};
                F{k}=reshape(F{k}, [blocksize{k} q]);
            end

        else
            for k=1:d
                F{k}=zeros(s(k), q);
                F{k}(I{k},:)=U{k};
                F{k}=reshape(F{k}, [blocksize{k} q]);
            end
        end

    case 'cell'
        for k=1:d
            F{k}=sparse(s(k), q);
            F{k}(I{k},:)=U{k};
            F{k}=mat2cell(F{k}, s(k), ones(1,q));
            F{k}=cellfun(@(x) reshape(x, blocksize{k}), F{k}, 'UniformOutput', false);
        end
end

if nargout==1 || nargout==d+1
    Mh=kron2mat(F{:});
end

switch nargout
    case 1
        varargout{1}=Mh;

    case d
        [varargout{1:d}]=F{:};

    case d+1
        varargout{1}=Mh;
        [varargout{2:d+1}]=F{:};

    otherwise
        error('Invalid number of output arguments')
end
end