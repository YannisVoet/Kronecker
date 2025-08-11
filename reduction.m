function[Lr]=reduction(L)

% REDUCTION: Reduces a set of integers by removing all elements which are 
% multiples of other elements.
%
% Lr = REDUCTION(L) computes the reduced set Lr from L.

m=length(L);

if m==0
    Lr=L;
else
    L=sort(L);
    Lr=L(1);
    L(mod(L, Lr(end))==0)=[];

    while ~isempty(L)
        Lr=[Lr L(1)];
        L(mod(L, Lr(end))==0)=[];
    end
end