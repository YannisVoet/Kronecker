function[f2]=kroncompr(n1, f1, n, f)

% KRONCOMPR (Kronecker Complete Right): Finds the linear indices of the
% right factor of a Kronecker product given the linear indices of the
% product and those of its left factor.
%
% fr = KRONCOMPR(f1, n1, f, n) computes the right linear indices from those
% of the Kronecker product (f) of size n and those of its left factor (f1)
% of size n1.

if isempty(f1) || isempty(f)
    f2=[];
    return
end

% Size of right factor
n2=n/n1;

% Reverse to index pairs
[i1,j1]=ind2sub([n1 n1], f1(1));
[i,j]=ind2sub([n n], f);

% Find the range of indices in f
ir=(i1-1)*n2+1:i1*n2;
jr=(j1-1)*n2+1:j1*n2;

is=ismember(i, ir);
js=ismember(j, jr);

s=is & js;

i2=i(s);
j2=j(s);

i2=mod(i2-1, n2)+1;
j2=mod(j2-1, n2)+1;

f2=sub2ind([n2 n2], i2, j2);
end