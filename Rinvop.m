function[M]=Rinvop(varargin)

% RINVOP: Computes the inverse rearrangement of the tensor R.
%
% M = RINVOP(R, blocksize) rearranges the tensor R of size
% m1n1 x m2n2 x ... x mdnd into a matrix M of size m1m2...md x n1n2...nd.
% The blocksize is specified as a cell array
% blocksize = {[m1 n1], [m2 n2], ..., [md nd]}. For the definition of the
% rearrangement operator, see [1] for d = 2 and [2] for d = 3. The tensor R
% is either an sptensor or tensor from the Tensor Toolbox or a standard
% multidimentional Matlab array.
%
% M = RINVOP(I, V, blocksize) provides instead a nnz(M) x d matrix I
% containing the position of nonzero entries of R as well as their value
% stored in the nnz(M) x 1 column vector V.
%
% References:
% [1] C. F. Van Loan and N. Pitsianis. Approximation with Kronecker products.
% In Linear Algebra for Large Scale and Real-Time Applications. Springer, 1993.
% [2] A. N. Langville and W. J. Stewart. A Kronecker product approximate
% preconditioner for SANs. Numerical Linear Algebra with Applications, 2004.
%
% See also ROP.

if nargin==2
    R=varargin{1};

    if isa(R, 'double') % Convert the multidimentional array to a tensor
        R=tensor(R);
    end

    [S,V]=find(R);
    blocksize=varargin{2};
else
    S=varargin{1};
    V=varargin{2};
    blocksize=varargin{3};
end

d=length(blocksize);
s=cat(1,blocksize{:});
nz=length(V);

Is=cell(1,d);
Js=cell(1,d);

% Convert linear indices to subscripts
for k=1:d
    [Is{k}, Js{k}]=ind2sub(blocksize{k}, S(:,k));
end

Is=fliplr(Is);
Js=fliplr(Js);

I=sub2ind(flipud(s(:,1))', Is{:});
J=sub2ind(flipud(s(:,2))', Js{:});

M=sparse(I, J, V, prod(s(:,1)), prod(s(:,2)), nz);
end