% 参数设置
L = 3;          % 空间周期
T = 1;           % 总时间
Nx = 1000;        % 空间网格点数
Nt = 500;        % 时间网格点数
dx = L / Nx;     % 空间步长
dt = T / Nt;     % 时间步长
r = dt / dx;     % Courant数

% 初始化网格
x = linspace(0, L, Nx+1);
x = x(1:end-1);  % 去掉最后一个点以保证周期性
u = zeros(Nx, Nt+1);

% 初始条件：周期性的正弦波
u(:, 1) = sin(2 * pi * x);

% Lax-Wendroff时间迭代
for n = 1:Nt
    un = u(:, n);
    unp1 = u(:, n+1);
    for i = 2:Nx-1
        unp1(i) = un(i) - r/2 * (un(i+1) - un(i-1)) + (r^2/2) * (un(i+1) - 2*un(i) + un(i-1));
    end
    % 周期边界条件
    unp1(1) = unp1(Nx-1);
    unp1(Nx) = unp1(2);
    u(:, n+1) = unp1;
end

acc_u3t = sin(2 * pi * 2);
u3t = u(end,end);
% 可视化结果
for n = 1:Nt+1
    plot(x, u(:, n));
    axis([0 L -1.5 1.5]);
    xlabel('x');
    ylabel('u(x,t)');
    title(sprintf('Time = %.2f', (n-1)*dt));
    pause(0.05);
end

