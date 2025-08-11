% Application 1: Space-time isogeometric discretization of a 3D problem.

clc
clear variables
close all

% Import matrices
Data=load('Data.mat');
Ks=Data.Ks; Wt=Data.Wt;
Ms=Data.Ms; Mt=Data.Mt;

B = kron(Wt,Ms) + kron(Mt,Ks);

[fact,~,F]=kronfact(B);

%% Validation checks
% Bandwidth
[lt, ut]=bandwidth(F{1}{1});
[ls, us]=bandwidth(F{1}{2});

%% Visualization
% Sparsity pattern
figure
spy(B)

%% Compare errors for various compatible pairs

b=fact{1}';
n=prod(b);
blocksize=num2cell([b b], 2);

b1=[31,144,6,2]';   blocksize1=num2cell([b1 b1], 2);
b2=[31,12,6,24]';   blocksize2=num2cell([b2 b2], 2);
b3=[62,6,6,24]';    blocksize3=num2cell([b3 b3], 2);

[Bh]=nkp(B,4,2,[],'blocksize', blocksize);
[B1]=nkp(B,4,2,[],'blocksize', blocksize1);
[B2]=nkp(B,4,2,[],'blocksize', blocksize2);
[B3]=nkp(B,4,2,[],'blocksize', blocksize3);

err=norm(B-Bh,'fro');
err1=norm(B-B1,'fro');
err2=norm(B-B2,'fro');
err3=norm(B-B3,'fro');

% GMRES solver
rhs=sparse(n,1);
rhs(1)=1;
restart=30;
tol=1e-8;
maxiter=10;

[x, res, niter]=gmresk(B, rhs, restart, tol, maxiter, Bh);
[x1, res1, niter1]=gmresk(B, rhs, restart, tol, maxiter, B1);
[x2, res2, niter2]=gmresk(B, rhs, restart, tol, maxiter, B2);
[x3, res3, niter3]=gmresk(B, rhs, restart, tol, maxiter, B3);

% Cheaper Kronecker rank 2 approximation
bt=b(2:end);
blocksizet=num2cell([bt bt], 2);
Kh=nkp(Ks,3,1,[],'blocksize', blocksizet);
Bt = kron(Wt,Ms) + kron(Mt,Kh);
errt=norm(B-Bt,'fro');

[xt, rest, nitert]=gmresk(B, rhs, restart, tol, maxiter, Bt);






