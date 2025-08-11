function[paths]=branch(paths, l, L)

% BRANCH: Expands the set of paths given a current path l and the set L of 
% left indices of length 2 factorizations.
%
% paths = BRANCH(paths, l, L) updates the set of paths.


el=l(end);
M=L(mod(L,el)==0 & L>el); % Elements in L that are multiples of the last element in l.
Mt=reduction(M); % Remove elements in M that are multiples of other elements.
mt=length(Mt);

ind=cellfun(@(x) x(end)==el, paths, 'UniformOutput', true);

path=paths{ind}; % Path that is updated
old_paths=paths(~ind); % Paths that are not updated
new_paths=cell(1,mt);

for k=1:mt
    new_paths{k}=[path Mt(k)];
end

% Do not update the paths if we have reached a leaf node.
if mt>0 % Update the paths
    paths=[new_paths old_paths];
end

for k=1:mt
    [paths]=branch(paths, [l Mt(k)], L);
end

end