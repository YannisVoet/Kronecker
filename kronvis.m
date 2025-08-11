function[varargout]=kronvis(A, varargin)

% KRONVIS: Kronecker graph visualization. Exploits the Kronecker structure
% of the adjancency matrix to produce a visually pleasing and easily
% interpretable graph of the network.
%
% KRONVIS(A) plots the Kronecker graph represented by the adjacency
% matrix A.
%
% KRONVIS(A, s) additionally specifies the (row) vector of sizes for the
% visualization. If this parameter is not provided, KRONVIS identifies
% suitable sizes from the Kronecker structure of the adjacency matrix.
%
% p = KRONVIS(A) returns the plot object for editing its properties.

if nargin==1
    [fact]=kronfact(A, false);

    if isempty(fact)
        s=size(A,1);
    else
        s=fact{1};
    end
else
    s=varargin{1};
end

d=length(s);
pts=cell(d,1); % Points
rot=cell(d,1); % Rotations
r=ones(d+1,1); % Radii

for k=1:d
    pts{k}=r(k)*exp((0:s(k)-1)*2*pi*1i/s(k));
    rot{k}=exp(pi*1i/2+angle(pts{k})*1i); % Additional rotation of pi/2 for esthetics
    r(k+1)=0.5*k/d*r(k)*abs(1-pts{k}(2));
end

J = length(pts{end});
P = reshape(pts{end},[J 1]);

for n = d-1:-1:1
    I = length(pts{n});
    H = reshape(pts{n},[1 I]);
    R = reshape(rot{n},[1 I]);
    P = bsxfun(@times,R,P);
    P = reshape(bsxfun(@plus,H,P),[I*J 1]);
    J = I*J;
end

X=real(P);
Y=imag(P);

figure
G = digraph(A);
p = plot(G, 'XData', X, 'YData', Y);
axis equal

varargout{1}=p;
end