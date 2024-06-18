function [x0,y0,p_x,p_y,e_x,e_y,jacob] = Mesh0()
AB = 1;         % Distance from left to right corner
AE = 2.4;       % Distance from bottom to top corner
DE = 4;         % Horizontal distance from bottom to top
angle = 15;     % Angle in degrees

% Corner points
x_l = 0;        % Left boundary x-coordinate
y_b = 0;        % Bottom boundary y-coordinate
x_B = AB;       % X-coordinate of bottom right corner
x_r = DE;       % X-coordinate of top right corner
y_t = AE;       % Top boundary y-coordinate
y_C = (DE - AB) * tan(angle * (pi / 180)); % Calculate y-coordinate of bottom right corner

% Grid resolution
nx = 80;        % Number of grid points in x-direction
ny = 50;        % Number of grid points in y-direction
nx1 = round(nx * (AB / DE));
nx2 = nx + 1 - nx1;
nx = nx1 + nx2 - 1;

% Generate domain
X_b1 = linspace(x_l, x_B, nx1);
X_b2 = linspace(x_B, x_r, nx2);
Y_b1 = ones(size(X_b1)) * y_b;
Y_b2 = y_b + (X_b2 - x_B) * tan(angle * (pi / 180));
X_t1 = linspace(x_l, x_B, nx1);
X_t2 = linspace(x_B, x_r, nx2);
X_t = linspace(x_l, x_r, nx);
X_bottom = [X_b1, X_b2(2:end)];
Y_bottom = [Y_b1, Y_b2(2:end)];
X_top = [X_t1, X_t2(2:end)];
Y_top = ones(size(X_top)) * y_t;

% Generate grid
x = zeros(nx, ny); % Preallocate x coordinate matrix
y = zeros(nx, ny); % Preallocate y coordinate matrix
p_x = zeros(nx, ny);
p_y = zeros(nx, ny);
e_x = zeros(nx, ny);
e_y = zeros(nx, ny);
jacob = zeros(nx, ny);
for i = 1:nx
    for j = 1:ny
        x(i, j) = X_bottom(i) + (X_top(i) - X_bottom(i)) * (j - 1) / (ny - 1);
        y(i, j) = Y_bottom(i) + (Y_top(i) - Y_bottom(i)) * (j - 1) / (ny - 1);
    end
end

% disp(y(end,:));
x0 = x;
y0 = y;
p_x(:,:)=1;
p_y(:,:)=0;


%for i = 1:nx
 %   for j = 1:ny
  %      if x(i,j)<=1
   %         e_x()

a = find(x(:,1) <= 1);
e_x(a,:) = 0;
e_y(a,:) = 1;
%disp(a);
a = find(x(:,1) > 1);
e_x(a,:) = (y(a,:)-2.4)*0.26795*2.4 ./(2.4-(x(a,:)-1)*0.26795) ./(2.4-(x(a,:)-1)*0.26795);
e_y(a,:) = 2.4 ./(2.4-(x(a,:)-1)*0.26795);
jacob(:,:)=p_x(:,:) .* e_y(:,:)-p_y(:,:) .* e_x(:,:);
end