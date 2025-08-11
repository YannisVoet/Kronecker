% Application 2: Visualization of Kronecker graphs

clc
close all

% Reference:
% [Leskovec et al., 2010] Leskovec, J., Chakrabarti, D., Kleinberg, J., 
% Faloutsos, C., and Ghahramani, Z. (2010). Kronecker graphs: an approach 
% to modeling networks. Journal of Machine Learning Research, 11(2).

%% Example from [Leskovec et al., 2010] (Figure 2)

clear variables
d=4;
A=cell(d,1);
A{1}=[1 1 0; 1 1 1; 0 1 1];
[A{2:end}]=deal(A{1});
A=kron2mat(A{:});

[fact]=kronfact(A);
kronvis(A,fact{1});

%% Example from [Leskovec et al., 2010] (Figure 3, Example 1)
clear variables
d=3;
A=cell(d,1);
A{1}=[1 1 1 1; 1 1 0 0; 1 0 1 0; 1 0 0 1];
[A{2:end}]=deal(A{1});

Ad3=kron2mat(A{1:3});
Ad2=kron2mat(A{1:2});

[factd3]=kronfact(Ad3);
[factd2]=kronfact(Ad2);

pd3=kronvis(Ad3,factd3{1});
pd2=kronvis(Ad2,factd2{1});
pd1=kronvis(A{1});

pd=[pd1 pd2 pd3];
set(pd, 'EdgeColor', [0.65 0.65 0.65], 'EdgeAlpha', 0.2)


%% Test 1
clear variables
A{1}=[1 1; 0 0];
A{2}=[1 1; 0 0];
A{3}=[1 1 1; 0 0 0; 1 1 1];
A=kron2mat(A{:});

[fact]=kronfact(A);

kronvis(A,fact{1});
kronvis(A,fact{2});
kronvis(A,fact{3});

%% Test 2
clear variables
A{1}=[1 0; 1 0];
A{2}=[1 0 1; 1 0 1; 1 0 1];
A=kron2mat(A{:});
A=kron(A,A);

[fact]=kronfact(A);
kronvis(A,fact{1});

%% Test 3
clear variables
A{1}=[1 0; 1 0];
A{2}=[1 0 1; 1 0 1; 1 0 1];
A=kron2mat(A{:});

[fact]=kronfact(A);
kronvis(A,fact{1});

%% Test 4
clear variables
A{1}=[1 1; 0 0];
A{2}=[1 1; 0 0];
A{3}=[1 1 1; 0 0 0; 1 1 1];
A=kron2mat(A{:});

kronvis(A,[4 3]);

%% Test 5
clear variables
A{1}=[1 0 0 0; 1 1 0 0; 1 1 1 0; 0 0 0 1];
A{2}=[1 0 0; 1 1 0; 0 0 1];
A{3}=[0 1; 1 0];
A=kron2mat(A{:});

kronvis(A,[4 3 2]);

%% Visualize graph construction

s=[4 3 2];

d=length(s);
pts=cell(d,1); % Points
rot=cell(d,1); % Rotations
r=ones(d+1,1); % Radii
c=cell(d+1,1); % Centers

for k=1:d
    pts{k}=r(k)*exp((0:s(k)-1)*2*pi*1i/s(k));
    rot{k}=exp(pi*1i/2+angle(pts{k})*1i); % Additional rotation of pi/2 for esthetics
    r(k+1)=0.5*k/d*r(k)*abs(1-pts{k}(2));
end



c{1}=0;
c{2}=pts{1}';

for k=0:d-2

    J = length(pts{d-k});
    P = reshape(pts{d-k},[J 1]);

    for n = d-k-1:-1:1
        I = length(pts{n});
        H = reshape(pts{n},[1 I]);
        R = reshape(rot{n},[1 I]);
        P = bsxfun(@times,R,P);
        P = reshape(bsxfun(@plus,H,P),[I*J 1]);
        J = I*J;
    end
    c{d-k+1}=P;
end


m=1e3;
c=cellfun(@(x) [real(x) imag(x)], c, 'UniformOutput', false);
colors=rand(d,3);
figure
hold on

for k=1:d
       x=c{k}(:,1)+r(k)*cos((1:m)*2*pi/m);
       y=c{k}(:,2)+r(k)*sin((1:m)*2*pi/m);

        plot(x',y', 'Color', colors(k,:), 'LineWidth', 2)
        plot(c{k}(:,1), c{k}(:,2), 'Color', colors(k,:), 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
end

% Plot graph vertices
plot(c{d+1}(:,1), c{d+1}(:,2), 'Color', 'k', 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
axis equal
box on
