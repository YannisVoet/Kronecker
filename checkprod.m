function[isprod, S1, S2]=checkprod(S)

% CHECKPROD: Checks whether a set S is a Cartesian product of two sets.
%
% isprod = CHECKPROD(S) returns a Boolean indicator that is true if S is a 
% Cartesian product and is false otherwise.
%
% [isprod, S1, S2] = CHECKPROD(S) if S is a Cartesian product, also returns
% the sets S1 and S2 such that S = S1 x S2.

%% Check if S is a Cartesian product

% If S is a Cartesian product of two sets Sx and Sy, of cardinality nx and
% ny, respectively, the cardinality of S must be nx*ny. This means that
% each distinct left index must appear ny times and each distinct right
% index must appear nx times. Use the first pair in S for guessing nx and
% ny. If nx*ny != nnz(A), S cannot be a Cartesian product.
% If nx*ny = nnz(A), then S is potentially a Cartesian product. Sort it
% (costly step) to confirm it.

ny=sum(S(:,1)==S(1,1));
nx=sum(S(:,2)==S(1,2));

if ny*nx ~= size(S,1)
    isprod=false;
    S1=[];
    S2=[];
else

    S=sortrows(S);

    R1=reshape(S(:,1), ny, nx);
    R2=reshape(S(:,2), ny, nx);

    S1=R1(1,:);
    S2=R2(:,1);

    isprod = all(R1==S1) & all(R2==S2);
end
