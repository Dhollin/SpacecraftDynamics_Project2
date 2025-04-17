% Spacecraft Dynamics Project Part 5-1
close all
clc

% Inertia Matrix (kg*m^2)
J = [500 0 0;
     0 400 -7;
     0 -7 440];

% Wheel Axial Inertia (kg*m^2)
I_w = 0.1;
J_w = diag([I_w, I_w, I_w]);

% Initial Angular Velocity (rpm -> rad/s)
omg_i = [5; 0; 0];
omg_i = omg_i * 2 * pi / 60;

% Initial MRPs <1
sigma_0 = [0.2; 0.2; 0.2];

%Initial Wheel Angular Velocities (start at rest we could also use a small value like 0.01)
omega_w0 = [0; 0; 0];

% Time span
tspan = [0, 60];

% Initial state vector [sigma; omega_s; omega_w]
y0 = [sigma_0; omg_i; omega_w0];

% X2 Gain Higher Control Gains for stronger response
D = 0.316 * diag([2400, 4400, 6200]);
K = diag([120, 220, 310]);

% Run ODE
[t, y] = ode45(@(t, y) MRP_dynamics(t, y, J, J_w, D, K, I_w), tspan, y0);

% Extract states
sigma = y(:, 1:3);
omega = y(:, 4:6);
omega_w = y(:, 7:9);

%% Lyapunov function calculation
V = zeros(length(t),1);
for i = 1:length(t)
    s = sigma(i,:)';
    w = omega(i,:)';
    V(i) = 0.5 * w' * J * w + 2 * trace(K) * log(1 + s' * s);  % Optional: use s' * K * s / (1 + s' * s)
end

%% Plot MRPs
figure
plot(t, sigma)
xlabel('Time [s]')
ylabel('MRPs')
legend('\sigma_1','\sigma_2','\sigma_3')
title('Modified Rodrigues Parameters')
grid on

%% Plot Angular Velocity
figure
plot(t, omega)
xlabel('Time [s]')
ylabel('\omega [rad/s]')
legend('\omega_1','\omega_2','\omega_3')
title('Spacecraft Angular Velocity')
grid on

%% Plot Reaction Wheel Speeds (RPM)
figure
plot(t, omega_w * 60 / (2*pi)) % rad/s to rpm
xlabel('Time [s]')
ylabel('Wheel Speed [rpm]')
legend('RW1','RW2','RW3')
title('Reaction Wheel Speeds')
grid on

%% Plot Lyapunov Function
figure
plot(t, V)
xlabel('Time [s]')
ylabel('V(t)')
title('Lyapunov Function')
grid on

%% Dynamics Function
function dydt = MRP_dynamics(~, y, J, J_w, D, K, I_w)
    sigma = y(1:3);
    omega_s = y(4:6);
    omega_w = y(7:9);

    % MRP kinematics
    sigma_dot = MRP_kinematics(omega_s, sigma);

    % Control torque from Lyapunov-based PD control law
    u = -D * omega_s - K * sigma;

    % Spacecraft angular acceleration (taken from Ben's code)
    %domg = J \ (-skew3(omega_s) * J * omega_s - skew3(omega_s) * I_w * omega_s - I_w * omega_s);

    %Spacecraft angular acceleration(updated to include the actual control
    %torque applied to the space craft)
    h_rw = J_w * omega_w;
    domg = J \ (-skew3(omega_s) * (J * omega_s + h_rw) + u);

    % Reaction wheel acceleration (user's original form using Omega)
    w1 = omega_w(1); w2 = omega_w(2); w3 = omega_w(3);
    Omega = [0 -w3 w2;
             w3 0 -w1;
            -w2 w1 0];
    omega_dot_wheels = inv(J_w) * (Omega * J_w * omega_w + u);

    dydt = [sigma_dot; domg; omega_dot_wheels];
end

%% MRP Kinematics Function
function sigma_dot = MRP_kinematics(omega, sigma)
    sigma_CPM = [0 -sigma(3) sigma(2);
                 sigma(3) 0 -sigma(1);
                -sigma(2) sigma(1) 0];
    A = (1/4) * ((1 - sigma' * sigma) * eye(3) + 2 * sigma_CPM + 2 * sigma * sigma');
    sigma_dot = A * omega;
end

%% Proper Skew-Symmetric Matrix
function skew = skew3(vec)
    skew = [    0     -vec(3)  vec(2);
              vec(3)     0    -vec(1);
             -vec(2)  vec(1)     0 ];
end
