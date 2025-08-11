function[varargout]=Rop(varargin)

% ROP: Computes the rearrangement of the matrix M.
%
% R = ROP(M, blocksize) rearranges the matrix M of size m1m2...md x
% n1n2...nd into a tensor R of size m1n1 x m2n2 x ... x mdnd. For d = 2,
% see [1], for d = 3, see [2] for their definition. The blocksize is
% specified as a cell array blocksize = {[m1 n1], [m2 n2], ..., [md nd]}.
% The output R is an sptensor (requiring the Tensor Toolbox). For converting
% R to a standard Matlab tensor, call double(R). For d = 2, R may be
% converted to a standard sparse Matlab matrix using spmatrix(R).
%
% R = ROP(I, J, V, blocksize) provides instead the position (I,J) of the
% nonzero entries V of M, all stored as nnz(M) x 1 column vectors.
%
% [I, V] = ROP(M, blocksize) returns instead a nnz(M) x d matrix I storing
% the position of nonzero entries V of R(M).
%
% References:
% [1] C. F. Van Loan and N. Pitsianis. Approximation with Kronecker products.
% In Linear Algebra for Large Scale and Real-Time Applications. Springer, 1993.
% [2] A. N. Langville and W. J. Stewart. A Kronecker product approximate
% preconditioner for SANs. Numerical Linear Algebra with Applications, 2004.
%
% See also RINVOP.

if nargin==2
    M=varargin{1};
    [I,J,V]=find(M);
    blocksize=varargin{2};
else
    I=varargin{1};
    J=varargin{2};
    V=varargin{3};
    blocksize=varargin{4};
end

s=cat(1,blocksize{:});
d=length(blocksize);

Is=cell(1,d);
Js=cell(1,d);

% Map i -> (i1,i2,...,id)
% Map j -> (j1,j2,...,jd)
[Is{:}]=ind2sub(flipud(s(:,1))', I);
[Js{:}]=ind2sub(flipud(s(:,2))', J);

Is=fliplr(Is);
Js=fliplr(Js);

Is=cat(2,Is{:});
Js=cat(2,Js{:});

% Position of nonzero entries of R(M)
L=(Js-1).*s(:,1)'+Is;

if nargout==1
    varargout{1}=sptensor(L, V, prod(s,2)');
else
    varargout{1}=L;
    varargout{2}=V;
end

end