function [r0,u0,v0,p0,tEnd,cfl] = Euler_IC2d(x,y)
%fprintf('Setting boundary conditions:\n');

% Pre-Allocate variables
r0 = zeros(size(x));
u0 = zeros(size(x));
v0 = zeros(size(x));
p0 = zeros(size(x));

r0(:,:) = 1.185;
u0(:,:) = 686.47;
v0(:,:) = 0;
p0(:,:) = 99719;

tEnd = 0.3;
cfl = 0.475;

%% Display
% h=surf(x,y,r0);
% set(h,'LineStyle','none');
% view(0,90);
