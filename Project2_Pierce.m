%% Project 2 - Pierce Elliott
clear; clc; close all;



%% Question 3.1

Iw = 0.1; %kg*m^2
initial = [5;0;0;0;0;0]./60.*2.*pi;
opts = odeset('MaxStep',0.5);
ts = 5000;
hw = 55;
Iww = [0,0,0;0,0,0;0,0,Iw];
J = [500,0,0;0,400,-7;0,-7,440];
[t,y] = ode45(@(t,w) angular_rates(w,J,hw,Iw,ts),[0,5000],initial,opts);
y2 = y;

nutation_angle = nutation(J,Iww,y,[0,0,1]);

figure(1);
tcl = tiledlayout(4,1);

nexttile(tcl)
hold on;
plot(t,y(:,1));
ylabel('W1 (Rad/s)',"Interpreter","latex","FontWeight","Bold")

nexttile(tcl)
hold on;
plot(t,y(:,2));
ylabel('W2 (Rad/s)',"Interpreter","latex","FontWeight","Bold")

nexttile(tcl)
hold on;
plot(t,y(:,3));
ylabel('W3 (Rad/s)',"Interpreter","latex","FontWeight","Bold")
% xlabel('Time (s)',"FontWeight","Bold")

nexttile(tcl)
hold on;
plot(t,nutation_angle);
ylabel('Nutation Angle (deg)',"FontWeight","Bold")
xlabel('Time (s)',"FontWeight","Bold")

figure(2)
plot3(y(:,1),y(:,2),y(:,3));
xlabel('W1 (rad/s)'); ylabel('W2 (rad/s)'); zlabel('W3 (rad/s)')
grid on;

wi = [5;0;0]./60.*2.*pi;
wf = [y(end,1); y(end,2); y(end,3)];
wwf = [y(end,4); y(end,5); y(end,6)];
H_i = norm(J*wi);
H_f = norm(J*wf+ Iww*wwf);

%% Question 3.2

Iw = 0.1; %kg*m^2
initial = [5;0;0;0;0;0]./60.*2.*pi;
opts = odeset('MaxStep',0.5);
ts = 5000;
hw = 55;
Iww = [0,0,0;0,0,0;0,0,Iw];
J = [500,0,0;0,400,0;0,0,440];
[t,y] = ode45(@(t,w) angular_rates(w,J,hw,Iw,ts),[0,5000],initial,opts);

nutation_angle = nutation(J,Iww,y,[0,0,1]);

figure(3);
tcl = tiledlayout(4,1);

nexttile(tcl)
hold on;
plot(t,y(:,1));
ylabel('W1 (Rad/s)',"Interpreter","latex","FontWeight","Bold")

nexttile(tcl)
hold on;
plot(t,y(:,2));
ylabel('W2 (Rad/s)',"Interpreter","latex","FontWeight","Bold")

nexttile(tcl)
hold on;
plot(t,y(:,3));
ylabel('W3 (Rad/s)',"Interpreter","latex","FontWeight","Bold")
% xlabel('Time (s)',"FontWeight","Bold")

nexttile(tcl)
hold on;
plot(t,nutation_angle);
ylabel('Nutation Angle (deg)',"FontWeight","Bold")
xlabel('Time (s)',"FontWeight","Bold")

figure(4)
plot3(y(:,1),y(:,2),y(:,3));
xlabel('W1 (rad/s)'); ylabel('W2 (rad/s)'); zlabel('W3 (rad/s)')
grid on;

wi = [5;0;0]./60.*2.*pi;
wf = [y(end,1); y(end,2); y(end,3)];
wwf = [y(end,4); y(end,5); y(end,6)];
H_i = norm(J*wi);
H_f = norm(J*wf+ Iww*wwf);

%% Question 4

% Initialize values for calculation
Bv1 = [-0.3;-0.1;0.9];
Bv2 = [0.8;-0.5;0.2];
Nv1 = [0;0;1];
Nv2 = [1;0;0];

C_BN_triad = TRIAD(Bv1,Bv2,Nv1,Nv2);

Bv1 = Bv1./norm(Bv1);
Bv2 = Bv2./norm(Bv2);
Nv1 = Nv1./norm(Nv1);
Nv2 = Nv2./norm(Nv2);
Bv = [Bv1,Bv2];
Nv = [Nv1,Nv2];
w = [2,1];

C_BN_quest = QUEST(Bv,Nv,w);

angle_diff = PR(C_BN_quest*C_BN_triad');

%% Question 5
Iw = 0.1; %kg*m^2
mrp = DCM2MRPs(C_BN_quest);
initial = [mrp;y2(end,:)'];
opts = odeset('MaxStep',2);
J = [500,0,0;0,400,-7;0,-7,440];
D = diag([700,700,700]);
K = 1000;
[t,s] = ode45(@(t,w) full_dyn(w,J,Iw,D,K,mrp),[0,50],initial,opts);

figure(5);
tcl = tiledlayout(3,1);

nexttile(tcl)
hold on;
plot(t,s(:,1));
ylabel('$\sigma _1$',"Interpreter","latex","FontWeight","Bold")
xlabel('Time (s)',"Interpreter","latex","FontWeight","Bold")

nexttile(tcl)
hold on;
plot(t,s(:,2));
ylabel('$\sigma _2$',"Interpreter","latex","FontWeight","Bold")
xlabel('Time (s)',"Interpreter","latex","FontWeight","Bold")

nexttile(tcl)
hold on;
plot(t,s(:,3));
ylabel('$\sigma _3$',"Interpreter","latex","FontWeight","Bold")
xlabel('Time (s)',"Interpreter","latex","FontWeight","Bold")

figure(6);
tcl = tiledlayout(3,1);

nexttile(tcl)
hold on;
plot(t,s(:,4));
ylabel('w1',"Interpreter","latex","FontWeight","Bold")
xlabel('Time (s)',"Interpreter","latex","FontWeight","Bold")

nexttile(tcl)
hold on;
plot(t,s(:,5));
ylabel('w2',"Interpreter","latex","FontWeight","Bold")
xlabel('Time (s)',"Interpreter","latex","FontWeight","Bold")

nexttile(tcl)
hold on;
plot(t,s(:,6));
ylabel('w3',"Interpreter","latex","FontWeight","Bold")
xlabel('Time (s)',"Interpreter","latex","FontWeight","Bold")

%% Question 5.3
Iw = 0.1;
J = [500,400,440];
A = [0,0,0,0.25,0,0,0,0,0;
     0,0,0,0,0.25,0,0,0,0;
     0,0,0,0,0,0.25,0,0,0;
     0,0,0,0,0,0,-Iw/J(1),0,0;
     0,0,0,0,0,0,0,-Iw/J(2),0;
     0,0,0,0,0,0,0,0,-Iw/J(3);
     0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0];

B = zeros(9,9);
B(7,7) = 1/Iw;B(8,8) = 1/Iw;B(9,9) = 1/Iw;

Q = diag([10,10,10,5,5,5,1,1,1]);
R = diag([0,0,0,0,0,0,1,1,1]);

% [K,S,P] = dlqr(A,B,Q,R);



%% Functions

% Full dynamics including control law
function [prime] = full_dyn(x,J,Iw,D,K,goal)
    
    % hw_dot = hw/ts;
    % ww_dot = hw_dot/Iw;
    % ww_dot = [0;0;ww_dot];
    sig = x(1:3);
    sig2 = sig'*sig;
    if sig2 > 1
        sig = -sig./sig2;
    end
    sig_e = ((1-norm(sig)^2)*goal - (1-norm(goal)^2)*sig + cross(2*goal,sig))...
            /(1+norm(sig)^2*norm(goal)^2 - 2*dot(sig,goal));

    w = x(4:6);
    ww = x(7:9);
    
    Iww = [Iw,0,0;0,Iw,0;0,0,Iw];

    u = -K*sig_e + D*w + skew(w)*J*w;

    ww_dot = inv(Iww)*(skew(ww)*Iww*ww+u);
    
    wdot = J\(-skew(w)*J*w-skew(w)*Iww*ww-Iww*ww_dot);

    sig_dot = 0.25*((1-sig2)*eye(3) + 2*skew(sig) + 2*sig*sig')*w;

    prime = [sig_dot;wdot;ww_dot];   
end

% Full dynamics without control law
function [prime] = lin_dyn(x,J,Iw,D,K,goal)
    
    % hw_dot = hw/ts;
    % ww_dot = hw_dot/Iw;
    % ww_dot = [0;0;ww_dot];
    sig = x(1:3);
    sig2 = sig'*sig;
    if sig2 > 1
        sig = -sig./sig2;
    end
    sig_e = ((1-norm(sig)^2)*goal - (1-norm(goal)^2)*sig + cross(2*goal,sig))...
            /(1+norm(sig)^2*norm(goal)^2 - 2*dot(sig,goal));

    w = x(4:6);
    ww = x(7:9);
    
    Iww = [Iw,0,0;0,Iw,0;0,0,Iw];

    % u = -K*sig_e + D*w + skew(w)*J*w;

    ww_dot = inv(Iww)*(skew(ww)*Iww*ww+u);
    
    wdot = J\(-skew(w)*J*w-skew(w)*Iww*ww-Iww*ww_dot);

    sig_dot = 0.25*((1-sig2)*eye(3) + 2*skew(sig) + 2*sig*sig')*w;

    prime = [sig_dot;wdot;ww_dot];   
end

function [xp,vs] = EP_dyn(x,J,D,K)
    qe = x(2:4);
    w = x(5:7);
    omega = skew(w);
    u = -omega*J*w - D*w - K*qe;
    wdot = inv(J)*(omega*J*w+u);
    Bdot = 0.5*[x(1), -x(2), -x(3), -x(4);
                x(2), x(1), -x(4), x(3);
                x(3), x(4), x(1), -x(2);
                x(4), -x(3), x(2), x(1)] * [0;w];

    xp = [Bdot;wdot];  
    vs = u;
end

% Define the differential kinematic equation for 3-2-1 Euler angles
function [rotation_rates] = angular_rates(w,J,hw,Iw,ts)
    
    hw_dot = hw/ts;
    ww_dot = hw_dot/Iw;
    ww_dot = [0;0;ww_dot];

    Iww = [0,0,0;0,0,0;0,0,Iw];
    
    wdot = J\(-skew(w(1:3))*J*w(1:3)-skew(w(1:3))*Iww*w(4:6)-Iww*ww_dot);

    rotation_rates = [wdot;ww_dot];
   
end

% Function to assemble skew matrix from a vector
function M = skew(s)
    M = [0,-s(3),s(2); s(3),0,-s(1); -s(2),s(1),0];
end

function angle = nutation(J,Iw,w,control)
    for i = 1:length(w(:,1))
        H(i,:) = (J*w(i,1:3)' + Iw*w(i,4:6)')';
        angle(i) = acosd(dot(H(i,:),control)/(norm(H(i,:))*norm(control)));
    end    
end

% Calculate attitude DCM using TRIAD
function C_BN = TRIAD(Bv1,Bv2,Nv1,Nv2)
    % Normalize all vectors
    Bv1 = Bv1./norm(Bv1);
    Bv2 = Bv2./norm(Bv2);
    Nv1 = Nv1./norm(Nv1);
    Nv2 = Nv2./norm(Nv2);
    % Calculate vectors in intermediate frame
    Bt1 = Bv1;
    Bt2 = cross(Bv1,Bv2)/norm(cross(Bv1,Bv2));
    Bt3 = cross(Bt1,Bt2);
    Nt1 = Nv1;
    Nt2 = cross(Nv1,Nv2)/norm(cross(Nv1,Nv2));
    Nt3 = cross(Nt1,Nt2);
    % Compute DCMs for both B and N frames
    C_BT = [Bt1,Bt2,Bt3];
    C_NT = [Nt1,Nt2,Nt3];
    % Compute attitude DCM
    C_BN = C_BT*C_NT';
end

% Calculate attitude DCM using QUEST
function C_BN = QUEST(Bv,Nv,w)    
    % Calculate intermediate values
    B = zeros(3,3);
    for i = 1:length(w)
        B = B + w(i)*Bv(:,i)*Nv(:,i)';
    end
    S = B + B';
    sig = trace(B);
    Z = [B(2,3)-B(3,2);B(3,1)-B(1,3);B(1,2)-B(2,1)];
    % Assemble K matrix
    K = [sig,Z';Z,S-sig*eye(3,3)];
    % Calculate max eigenvalue using Newton Raphson
    eig_val = NR(sum(w),1e-10);
    
    % Calculate Rodrigues parameters (eq. 3.238)
    q = inv((eig_val+sig)*eye(3,3)-S)*Z;
    
    % Assemble DCM from Rodrigues parameters
    C_BN = (1/(1+q'*q))*[1+q(1)^2-q(2)^2-q(3)^2, 2*(q(1)*q(2)+q(3)), 2*(q(1)*q(3)-q(2))
                      2*(q(1)*q(2)-q(3)), 1-q(1)^2+q(2)^2-q(3)^2, 2*(q(2)*q(3)+q(1))
                      2*(q(1)*q(3)+q(2)), 2*(q(2)*q(3)-q(1)), 1-q(1)^2-q(2)^2+q(3)^2];
    
    % Eq 3.236 for Newton Raphson root solving
    function y = f(s)
        y = det(K-s*eye(4,4));
    end
    % 1st derivative of Eq 3.236
    function dy = fp(s)
        dy = -det(K-s*eye(4,4))*trace(inv(K-s*eye(4,4)));
    end
    % Root solving using Newton Raphson for the eigenvalue
    function eig_val = NR(initial,tol)
        e = initial;
        err = f(e);
        while err > tol
            e = e - f(e)/fp(e);
            err = f(e);
        end
        eig_val = e;
    end
end

% Find principle rotation angle of a DCM
function angle = PR(C)
    angle = acosd(0.5*(trace(C)-1));
end

function mrp = DCM2MRPs(DCM)
    c = sqrt(trace(DCM)+1);
    mrp = 1/(c*(c+2)) *[DCM(2,3)-DCM(3,2);DCM(3,1)-DCM(1,3);DCM(1,2)-DCM(2,1)];
end




