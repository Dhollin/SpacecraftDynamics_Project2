% Clearing Workspace and Closing Open Figures
clear
close all

% Adding function folder to the file path
addpath('C:\Users\Ben\Google Drive\College\Senior\Spring Semester\Spacecraft Dynamics\Functions')

%% Definitions
% Problem 1
syms J J_s J_w 3 matrix
syms omg domg omg_w domg_w [3 1] matrix
% J = Inertia Matrix, 3x3
% J_s = Satellite Inertia Matrix, 3x3
% omg = Angular Velocity of Satellite WRT Inertial Frame, 3x1
% domg = Time Derivative of Angular Velocity of Satellite WRT Inertial Frame, 3x1
% J_w = Reaction Wheel Inertia Matrix, 3x3
% omg_w = Wheel Rotation Speed WRT Body Frame, 3x1
% domg_w = Time Derivative of Wheel Rotation Speed WRT Body Frame, 3x1

% Problem 2
syms I_w
syms b_3 Omg_w dOmg_w [3 1] matrix
% b_3 = Reaction Wheel Axis, 3x1
% I_w = Wheel Axial Inertia, 1x1
% Omg_w = Wheel Speed, 3x1
% dOmg_w = Wheel Acceleration, 3x1

%% Question 1
% Total Angular Momentum About CoM
h = J*omg+J_w*omg_w;
% Euler's Rotational Equations of Motion

%% Question 2
% Euler's Rotational Equations of Motion

%% Question 3-1
% Given
% Inertia Matrix (kg*m^2)
J = [500 0 0;
    0 400 -7;
    0 -7 440];
% Wheel Axial Inertia (kg*m^2)
I_w = .1;
% Initial Angular Velocity of Satellite WRT Inerial Frame (rpm)
omg_i = [5; 0; 0];
% Initial Wheel Momentum (kg*m^2/sec)
h_i = 0;
% Nominal Wheel Momentum (kg*m^2/sec)
h_n = 55;
% Spin-up Duration (sec)
t_s = 5000;
% Wheel Control (kg*m^2/sec^2)
dh_w3 = h_n/t_s;

% Calculations
% Converting omg_i to rad/s
omg_i = omg_i*2*pi/60;
omg_init = [omg_i; 0];
domg_w3 = dh_w3/I_w;
J_w = diag([0,0,I_w]);

domg_w = [0; 0; domg_w3];

[t, Omg] = ode45(@ (t, Omg) ang_vel(t, Omg, J, J_w, domg_w), [0 t_s], omg_init);
b3 = [0;0;1];
beta = zeros(size(t));


for i=1:length(t)
    h_vec = J*Omg(i,1:3).'+J_w*[0; 0; Omg(i,4)];
    beta(i) = acosd(dot(h_vec,b3)/(norm(h_vec)*norm(b3)));
end

figure(1)
fig1 = tiledlayout("vertical");
title(fig1, "Time Responce of the Angular Velocity of the Satellite",'Off-Diagonal Elements Non-Zero')
xlabel(fig1, 'Time (sec)')
ylabel(fig1, 'Angular Velocity, \omega (rad/s)')

nexttile
hold on
title('\omega_1')
grid on
plot(t,Omg(:,1))
hold off

nexttile
hold on
title('\omega_2')
grid on
plot(t,Omg(:,2))
hold off

nexttile
hold on
title('\omega_3')
grid on
plot(t,Omg(:,3))
hold off

figure(2)
hold on
title('Time Responce of the Nutation of the Satellite','Off-Diagonal Elements Non-Zero')
plot(t,beta)
grid on
hold off

figure(3)
plot3(Omg(:,1), Omg(:,2), Omg(:,3))
hold on
grid on
title('3D Plot of Angular Velocity of the Satellite','Off-Diagonal Elements Non-Zero')
hold off

%% Question 3-2
% Given
% Inertia Matrix (kg*m^2)
J = [500 0 0;
    0 400 0;
    0 0 440];
% Wheel Axial Inertia (kg*m^2)
I_w = .1;
% Initial Angular Velocity of Satellite WRT Inerial Frame (rpm)
omg_i = [5; 0; 0];
% Initial Wheel Momentum (kg*m^2/sec)
h_i = 0;
% Nominal Wheel Momentum (kg*m^2/sec)
h_n = 55;
% Spin-up Duration (sec)
t_s = 5000;
% Wheel Control (kg*m^2/sec^2)
dh_w3 = h_n/t_s;

% Calculations
% Converting omg_i to rad/s
omg_i = omg_i*2*pi/60;
omg_init = [omg_i; 0];
domg_w3 = dh_w3/I_w;
J_w = diag([0,0,I_w]);

domg_w = [0; 0; domg_w3];

[t, Omg] = ode45(@ (t, Omg) ang_vel(t, Omg, J, J_w, domg_w), [0 t_s], omg_init);
b3 = [0;0;1];
beta = zeros(size(t));

for i=1:length(t)
    h_vec = J*Omg(i,1:3).'+J_w*[0; 0; Omg(i,4)];
    beta(i) = acosd(dot(h_vec,b3)/(norm(h_vec)*norm(b3)));
end

figure(4)
fig1 = tiledlayout("vertical");
title(fig1, "Time Responce of the Angular Velocity of the Satellite",'Off-Diagonal Elements Zero')
xlabel(fig1, 'Time (sec)')
ylabel(fig1, 'Angular Velocity, \omega (rad/s)')

nexttile
hold on
title('\omega_1')
grid on
plot(t,Omg(:,1))
hold off

nexttile
hold on
title('\omega_2')
grid on
plot(t,Omg(:,2))
hold off

nexttile
hold on
title('\omega_3')
grid on
plot(t,Omg(:,3))
hold off

figure(5)
hold on
title('Time Responce of the Nutation of the Satellite','Off-Diagonal Elements Zero')
plot(t,beta)
grid on
hold off

figure(6)
plot3(Omg(:,1), Omg(:,2), Omg(:,3))
hold on
grid on
title('3D Plot of Angular Velocity of the Satellite','Off-Diagonal Elements Zero')
hold off

%% Functions
function [skew] = skew3(vec)
% Returns the 3x3 skew matrix for a given vector
    skew = cross([vec vec vec],eye(3));
end

function [dOmg] = ang_vel(t, Omg, J, I_w, domg_w)
    omg = Omg(1:3);
    omg_w = [0; 0; Omg(4)];
    domg = J\(-skew3(omg)*J*omg-skew3(omg)*I_w*omg_w-I_w*domg_w);
    dOmg = [domg; domg_w(3)];
end