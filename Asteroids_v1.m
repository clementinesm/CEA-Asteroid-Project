clear;
clc;

n = input('number of bodies = '); % Number of bodies
G = 1; % Gravitational constant
dt = 0.1;
T = 10;

for i = 1:n
    m(i) = input('mass = '); % Mass vector
    sx(i) = input('x-pos = '); % x-position vector
    sy(i) = input('y-pos = '); % y-position vector
    sz(i) = 0; % z-position vector
    vx(i) = input('x-vel = '); % x-velocity vector
    vy(i) = input('y-vel = '); % y-velocity vector
    vz(i) = 0; % z-velocity vector
end
x(:,:,1) = [transpose(sx),transpose(sy),transpose(sz)]; % Combine x and y-positions
v_old = [transpose(vx),transpose(vy),transpose(vz)]; % Combine x- and y-velocities

for t = 1:(T/dt)
    for i = 1:3 % Determines which dimension we're looking in (i = 1:x, i = 2:y, etc.)
        for j = 1:n
            for k = 1:n
                if j == k
                    A(j,k) = 0; % i (dimension) comes last
                else
                    A(j,k) = G * m(k) * (x(k,i,t) - x(j,i,t))/ ((x(j,1,t) - x(k,1,t))^2 + (x(j,2,t) - x(k,2,t))^2)^(3/2);
                end
            end
        end
    
        a(:,i) = diag(A*ones(n));
    end
    
    v_new = v_old + a * dt;
    x(:,:,t+1) = x(:,:,t) + v_old * dt;
    v_old = v_new;
end

for i = 1:n
    plot_x(:,:) = x(i,:,:);
    plot(plot_x(1,:),plot_x(2,:))
    hold on
end