% Examples from the article

%% Example 2.5
clear variables

A{1}=[1 1; 0 0];
A{2}=[1 1; 0 0];
A{3}=[1 1 1; 0 0 0; 1 1 1];
A=kron2mat(A{:});

[fact,l,M]=kronfact(A);

%% Example 2.6
clear variables

A{1}=[1 0; 0 0];
A{2}=[1 0; 1 0];
A{3}=[1 0 0; 0 0 0; 0 0 0];
A=kron2mat(A{:});

[fact,l,M]=kronfact(A);

%% Example 3.15
clear variables

n=24;
L=[2 3 4 6 8 12];

fact=kronfact({L, n});

%% Example 3.16
clear variables

A{1}=[1 0; 1 0];
A{2}=[1 0 1; 1 0 1; 1 0 1];
A=kron2mat(A{:});
A=kron(A,A);

[fact,l,M]=kronfact(A);

%% Test 1
clear variables

A{1}=[1 1 0; 0 1 1; 0 0 1];
A{2}=[1 1; 1 1];
A=kron2mat(A{:});

fact=kronfact(A);











