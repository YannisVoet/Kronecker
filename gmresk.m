function[x, varargout] = gmresk(A, b, varargin)

% GMRESK: (Restarted) Generalized Minimal Residual Method (GMRES(k)) for 
% solving linear systems of equations.
%
% x = GMRESK(A, b) iteratively solves the linear system A*x = b where A is
% a n x n matrix and b is a column vector of length n. The matrix A is
% either provided explicitly or as a function handle.
%
% x = GMRESK(A, b, restart) restarts the method every restart iterations.
% If this parameter is n or [] the non-restarted method is used.
% Default: n (non-restarted method).
%
% x = GMRESK(A, b, restart, tol) specifies the tolerance on the absolute 
% residual. For a stopping criterion based on the relative residual, 
% replace tol with tol*norm(b). Default: 1e-8.
%
% x = GMRESK(A, b, restart, tol, maxiter) also specifies the maximum number
% of outer iterations (i.e. the maximum number of restarts). If the method
% is not restarted, the maximum number of total iterations is maxiter
% (instead of maxiter*restart).
% Default: min(n,10).
%
% x = GMRESK(A, b, restart, tol, maxiter, M) also specifies a preconditioning
% matrix. The preconditioning matrix M is either provided explicitly or as
% a function handle. The algorithm uses the right preconditioned version of
% GMRES (see [1]).
%
% x = GMRESK(A, b, restart, tol, maxiter, M, x0) also specifies the initial
% starting vector. Default: all zero vector.
%
% x = GMRESK(A, b, restart, tol, maxiter, M, X0, reorth_tol) also specifies
% the reorthogonalization tolerance for the modified Gram-Schmidt process.
% Default: 0.7.
%
% [x, res] = GMRESK(A, b, ...) returns the vector of absolute residuals for
% each iteration.
%
% [x, res, niter] = GMRESK(A, b, ...) also returns a vector with the number
% of outer, inner and total number of iterations.
% 1 ≤ niter(1) ≤ maxiter
% 1 ≤ niter(2) ≤ restart
% 1 ≤ niter(3) ≤ maxiter*restart (or maxiter for the non-restarted method)
%
% Reference:
% [1] Y. Saad. Iterative methods for sparse linear systems. SIAM, 2003.
%
% See also GLGMRESK, BICGSTB, CG.

%% Set algorithm parameters

% Set default parameters
Default{1}=length(b);
Default{2}=1e-8;
Default{3}=min(length(b),10);
Default{4}=@(x) x;
Default{5}=zeros(size(b));
Default{6}=0.7;

% Replace empty inputs with default parameters
def=cell2mat(cellfun(@isempty, varargin, 'UniformOutput', false));
[varargin{def}]=Default{def};

Param = inputParser;
Param.addRequired('A', @(x) isa(x, 'function_handle') || isa(x, 'double'));
Param.addRequired('b');
Param.addOptional('restart',    Default{1}, @(x) isscalar(x) & x > 0);
Param.addOptional('tol',        Default{2}, @(x) isscalar(x) & x > 0);
Param.addOptional('maxiter',    Default{3}, @(x) isscalar(x) & x > 0);
Param.addOptional('M',          Default{4}, @(x) isa(x, 'function_handle') || isa(x, 'double'));
Param.addOptional('x0',         Default{5}, @(x) all(size(x)==size(b)));
Param.addOptional('reorth_tol', Default{6}, @(x) isscalar(x) & x > 0);
Param.parse(A, b, varargin{:});

%% Retrieve parameters
restart=Param.Results.restart;
tol=Param.Results.tol;
maxiter=Param.Results.maxiter;
M=Param.Results.M;
x0=Param.Results.x0;
reorth_tol=Param.Results.reorth_tol;

if isa(A, 'double')
    A=@(x) A*x;
end

if isa(M, 'double')
    dM=decomposition(M);
    M=@(x) dM\x;
end

if restart==length(b)
    restart=maxiter;
    maxiter=1;
end

%% GMRES method

res = cell(maxiter,1);

for outerit=1:maxiter

    [x,beta,innerit]=innergmres(A,M,b,x0,restart,tol,reorth_tol);
    res{outerit}=beta;

    if beta(end)<tol
        break
    end
    x0 = x;
end

varargout{1}=sort(uniquetol(cat(1,res{:}),1e-14), 'descend');
varargout{2}=[outerit innerit (outerit-1)*restart+innerit];


    function[x,beta,k]=innergmres(A,M,b,x0,m,tol,reorth_tol)

        r0=b-A(x0);
        n=length(r0);
        U=zeros(n, m);
        H=zeros(m+1, m);
        beta=zeros(m+1,1);
        beta(1)=norm(r0);
        U(:,1)=r0/beta(1);


        for k=1:m
            w=M(U(:,k));
            w=A(w);

            h=U(:,1:k)'*w;
            u_tilde=w-U(:,1:k)*h;

            if norm(u_tilde-w)<reorth_tol % Re-orthogonalization
                h_hat=U(:,1:k)'*u_tilde;
                h=h+h_hat;
                u_tilde=u_tilde-U(:,1:k)*h_hat;
            end

            h_tilde=norm(u_tilde);
            U(:,k+1)=u_tilde/h_tilde;
            H(1:(k+1) ,k)=[h; h_tilde];

            % Be careful: parentheses are crucial
            y=beta(1)*(H(1:k+1,1:k)\eye(k+1,1));
            beta(k+1)=norm(beta(1)*eye(k+1,1)-H(1:k+1,1:k)*y);

            if beta(k+1)<tol
                break
            end

        end
        x=x0+M(U(:,1:k)*y);
        beta=beta(1:k+1);

    end
end