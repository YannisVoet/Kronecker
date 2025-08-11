function[fact, lin, M]=kronfact(A, vis)

% KRONFACT: Seeks a Kronecker product factorization of a square binary matrix.
%
% F = KRONFACT(A) returns a cell array containing the sizes of the factor
% matrices in each prime factorization of A. If A does not admit any
% factorization, F is an empty cell.
%
% F = KRONFACT({L, n}) specifies instead the left indices of all possible
% length 2 factorizations and the size n of the binary matrix. Note: this
% information is insufficient for recovering the sparsity pattern of the
% factor matrices.
%
% [F, l] = KRONFACT(A) also returns the linear indices for the factor
% matrices in the decompositions (only if A is supplied explicitly).
%
% [F, l, M] = KRONFACT(A) also returns the factor matrices in the
% decompositions. M{k}{i} is the ith factor matrix in the kth prime
% decomposition (only if A is supplied explicitly).
%
% F = KRONFACT(A, vis) specifies whether to visualize the decomposition
% graph. Default: true.

arguments
    A
    vis (1,1) {mustBeNumericOrLogical} = true
end

%% Find the set of length 2 factorizations

if iscell(A)
    L=A{1};
    n=A{2};
    nl=length(L);

    Fl=cell(1,nl);
    Fr=cell(1,nl);
    expl=false;

else % Compute the set of left indices from the sparsity pattern of A.

    A=A~=0;
    [n,m]=size(A);
    expl=true;

    if n~=m
        error('This implementation only supports square matrices.')
    end

    if isprime(n)
        fact={};
        return
    end

    d=divisors(n);
    d=setdiff(d, [1, n]);
    n1=d; n2=n./d;
    nd=length(d);

    [I,J]=find(A);

    [Ir, ~, ir]=unique(I);
    [Jr, ~, jr]=unique(J);

    ir2=mod(Ir-1, n2)+1;
    jr2=mod(Jr-1, n2)+1;

    ir1=(Ir-ir2)./n2+1;
    jr1=(Jr-jr2)./n2+1;

    i1=ir1(ir,:); j1=jr1(jr,:);
    i2=ir2(ir,:); j2=jr2(jr,:);

    l1=(j1-1).*n1+i1;
    l2=(j2-1).*n2+i2;

    L=[];
    i=1;

    for k=1:nd

        [isfact, fl, fr]=checkprod([l1(:,k), l2(:,k)]);

        if isfact
            L=[L d(k)];
            Fl{i}=fl;
            Fr{i}=fr;
            i=i+1;
        end
    end
end

%% Find factorizations of greater length
Lt=reduction(L); % Roots (all factorizations start with these indices)
paths=num2cell(Lt);
lt=length(Lt);

for k=1:lt
    [paths]=branch(paths, Lt(k), L);
end

np=length(paths);   % Number of paths = number of branches
E=cell(np,1);       % Path edges
p=cell(np,1);       % p-indices
edgeids=cell(np,1); % Edge IDs
code=cell(np,1); % Edge path code: all edges belonging to the kth path are labeled with an index k
fact=cell(np,1); % Sizes of the factor matrices in the decompositions
lin=cell(np,1);  % Linear indices of each factor matrix in each decomposition
M=cell(np,1);    % Factor matrices in the decompositions
s=1;

for k=1:np
    path=paths{k};
    % Length of the path = nploc >= 1
    % Length of the factorization = nploc + 1 >= 2
    % Number of edges = nploc - 1 >= 0
    nploc=length(path);
    lin{k}=cell(1,nploc+1);       % lin{k} stores the linear indices of each factor matrix in the kth decomposition
    lin{k}{1}=Fl{L==path(1)};     % Initialize leftmost factor
    lin{k}{end}=Fr{L==path(end)}; % Initialize rightmost factor

    if nploc>1
        E{k}=num2cell([path(1:end-1); path(2:end)]', 2);
        p{k}=path(2:end)./path(1:end-1);
        edgeids{k}=s:(s+nploc-2);
        code{k}=k*ones(1,nploc-1);
        s=s+nploc-1;

        for j=1:nploc-1 % Find the linear indices for the nploc - 1 inner factors
            lin{k}{j+1}=kroncompr(path(j), lin{k}{j}, path(j+1), Fl{L==path(j+1)});
        end
    end

    fact{k}=[path(1) p{k} n/path(end)];

    % Recover the factor matrices in the decomposition
    if expl
        M{k}=cell(1,nploc+1);
        for j=1:nploc+1
            [I,J]=ind2sub([fact{k}(j) fact{k}(j)], lin{k}{j});
            M{k}{j}=sparse(I, J, ones(length(lin{k}{j}),1), fact{k}(j), fact{k}(j));
        end
    end
end

E=cell2mat(cat(1,E{:}));
p=cat(2,p{:})';
code=cat(2,code{:})';
edgeids=cat(2,edgeids{:})';

%% Visualization

if vis
    figure

    if size(E,1)>0
        EdgeTable = table(E, p, code, edgeids, 'VariableNames',{'EndNodes', 'Weight', 'Code', 'IDs'});
        G=digraph(EdgeTable);
        H=subgraph(G, L);


        plt=plot(H, 'NodeLabel', L, 'EdgeLabel', H.Edges.Weight);
        plt.NodeFontSize=15;
        plt.EdgeFontSize=15;
        plt.MarkerSize=10;
        plt.LineWidth=2.5;
        plt.EdgeAlpha=1;
        colors=rand(np,3);

        for k=1:np
            highlight(plt, 'edges', find(H.Edges.Code==k), 'EdgeColor', colors(k,:))
        end

    else % The graph only contains isolated nodes

        G=digraph([], [], [], n);
        H=subgraph(G, L);
        plt=plot(H, 'NodeLabel', L);
        plt.NodeFontSize=15;
        plt.EdgeFontSize=15;
        plt.MarkerSize=10;
        plt.LineWidth=2.5;
        plt.EdgeAlpha=1;
    end
end