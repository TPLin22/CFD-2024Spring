% 参数设置
L = 3;          % 空间周期
T = 1;           % 总时间
Nx = 100;        % 空间网格点数
Nt = 500;        % 时间网格点数
dx = L / Nx;     % 空间步长
dt = T / Nt;     % 时间步长
r = dt / dx;     % Courant数

% 初始化网格
x = linspace(0, L, Nx+1);
x = x(1:end-1);  % 去掉最后一个点以保证周期性
u = zeros(Nx, Nt+1);

max_u_values = zeros(1, Nt+1);

% 初始条件：周期性的正弦波
u(:, 1) = sin(2 * pi * x);
max_u_values(1) = max(u(:, 1));  % 存储初始条件下的最大u值

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
    max_u_values(n+1) = max(u(:, n+1));
end

% 观察有无数值耗散
%plot(max_u_values);
%xlabel('Time Step');
%ylabel('Maximum u value');
%title('Maximum Value of u at Each Time Step');

acc_u3t = sin(2 * pi * (3-T));
u3t = u(end,end);
err = abs(acc_u3t - u3t);

% 可视化解的结果
for n = 1:Nt+1
    plot(x, u(:, n));
    axis([0 L -1.5 1.5]);
    xlabel('x');
    ylabel('u(x,t)');
    title(sprintf('Time = %.2f', (n-1)*dt));
    pause(0.05);
end

% 可视化T=100的数值解与精确解以观察相位差
%acc_sol = sin(2 * pi * (x - 100));

%plot(x, acc_sol); hold on;
%plot(x, u(:, end)); hold off;
%legend('accurate solution', 'T=100');
%title('精确解与数值解 相位对比');
%xlabel('x');
%ylabel('u');
%grid on;




