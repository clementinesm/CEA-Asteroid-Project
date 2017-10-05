clear;
clc;

n = 2; % Number of bodies
G = 10; % Gravitational constant

for i = 1:n
    m(i) = input('mass = '); % Mass vector
    sx(i) = input('x-pos = '); % x-position vector
    sy(i) = input('y-pos = '); % y-position vector
    sz(i) = 0; % z-position vector
    vx(i) = input('x-vel = '); % x-velocity vector
    vy(i) = input('y-vel = '); % y-velocity vector
    vz(i) = 0; % z-velocity vector
end
x = [transpose(sx),transpose(sy),transpose(sz)]; % Combine x and y-positions
v = [transpose(vx),transpose(vy),transpose(vz)]; % Combine x- and y-velocities

for i = 1:3 % Determines which dimension we're looking in (i = 1:x, i = 2:y, etc.)
    for j = 1:n
        for k = 1:n
            if j == k
                A(j,k) = 0; % i (dimension) comes last
            else
                A(j,k) = G * m(k) * (x(k,i) - x(j,i))/ ((x(j,1) - x(k,1))^2 + (x(j,2) - x(k,2))^2)^(3/2);
            end
        end
    end
    
    a(:,i) = diag(A*ones(n));
end

