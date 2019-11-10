clear; clc;
% 3 Dofs manipularot
%Given:
syms theta1;        
syms theta2;
syms d3;

syms q1;
syms q2;
syms q3;

syms d1;
syms a2;


%% Forward Kinematics
theta1 = deg2rad(0) ; theta2 = deg2rad(0) ; d3 = 10; 


T = simplify(trotz(q1)*transl(0, 0, d1)*troty(-q2)*transl(a2, 0, 0)*transl(q3,0, 0))

px = T(1,4);
py = T(2,4);
pz = T(3,4);

T1 = trotz(theta1)*transl(0, 0, d1)*troty(-theta2)*transl(a2, 0, 0)*transl(d3,0, 0)


%% Inverse Kinematics
px_i = 50;
py_i = 40;
pz_i = 50;

a2i = 10;
d1i = 20;

r = sqrt(px_i^2 + py_i^2);
s = pz - d1i;


theta1_i = rad2deg(atan2(py_i, px_i))
d3_i = sqrt(r^2 + s^2) - a2i
theta2_i = rad2deg(atan2(s, r))

%% Jacobians

% Classic approach
dpx = diff(px, q1) + diff (px, q2) + diff(px, q3)  
dpy = diff(py, q1) + diff (py, q2) + diff(py, q3)
dpz = diff(pz, q1) + diff (pz, q2) + diff(pz, q3)

JV = [diff(px, q1), diff(px, q2), diff(px, q3);
      diff(py, q1), diff(py, q2), diff(py, q3);
      diff(pz, q1), diff(pz, q2), diff(pz, q3); ]
J_om = [0, sin(q1), 0;
        0, -cos(q1), 0;
        1, 0, 0];
Jacob = simplify([JV; J_om])
Jacob1 = subs(Jacob, {a2, d1}, {10,20})
%% Geometric approach

R_00 = rotz(0);
R_10 = rotz(q1)*rotx(pi/2);
R_20 = R_10*rotz(q2)*roty(pi/2)*rotz(pi/2);
R_30 = R_20;

O3 = trotz(q1)*transl(0, 0, d1)*troty(-q2)*transl(a2, 0, 0)*transl(q3,0, 0);
O2 = trotz(q1)*transl(0, 0, d1)*troty(-q2)*transl(a2, 0, 0);
O1 = trotz(q1)*transl(0, 0, d1);

d_30 = [px; py;pz];
d_20 = [O2(1,4); O2(2,4); O2(3,4)];
d_10 = [0 ; 0 ; d1];
d_00 = [0;0;0];

uz = [0;0;1];
Jv_g = simplify([cross(R_00*uz, (d_30 - d_00)), cross(R_10*uz, (d_30 - d_10)), R_20*uz]);
J_omg = simplify([R_00*uz, R_10*uz, d_00])
Jacobian_G = simplify([Jv_g;J_omg])
%% Singularities
Sing = Jacob((1:3),:);
DS = simplify(det(Sing));

%% Velocity of the tool frame
syms theta1 theta2 d3 t 
theta1 = sin(t);
theta2 = cos(2*t);
d3 = sin (3*t);

dtheta1 = cos(t);
dtheta2 = -2*sin(2*t);
dd3 = 3*cos(3*t);

J = simplify(subs(Jacob, {d1, a2, q1, q2, q3}, {20, 10, theta1, theta2, d3}));

q = [theta1; theta2; d3];
qdot = [dtheta1; dtheta2; dd3];

xi = simplify(J * qdot);

time = 0 : 0.1 :10;
xi_t = subs(xi, {t}, {time});
xi_t = double(xi_t);
%%
%Plot
figure(1);
plot(time, xi_t(1:3, :), 'LineWidth', 2)
xlabel('Time of moving')
ylabel('Linear velocity value')
title('Linear velocities plots')
legend('v_x', 'v_y', 'v_z', 'Fontsize', 15)
grid on

figure(2);
plot(time, xi_t(4:6, :), 'LineWidth', 2)
xlabel('Time of moving')
ylabel('Angular velocity value')
title('Angular velocities plots')
grid on
legend('\omega_x', '\omega_y', '\omega_z', 'Fontsize', 14)

%% Velocity of the tool frame
% Inverse kinematics
sym t 
a2n = 10;
d1n = 20;
px_d = 2*a2n*cos(t);
py_d = 2*a2n*cos(2*t);
pz_d = d1n*sin(3*t);

r_d = sqrt(px_d^2 + py_d^2);
s_d = pz_d - d1n;

%Case1
theta1_d = atan2(py_d, px_d);
d3_d = sqrt(r_d^2 + s_d^2) - a2n;
theta2_d =atan2(s_d, r_d);

time = 0 : 0.1 :10;
theta1_d_time = subs(theta1_d, {t}, {time});
theta1_d_t = double(theta1_d_time);

theta2_d_time = subs(theta2_d, {t}, {time});
theta2_d_t = double(theta2_d_time);

d3_d_time = subs(d3_d, {t}, {time});
d3_d_t = double(d3_d_time);

%Plot
figure(3);
plot(time, [theta1_d_t; theta2_d_t; d3_d_t], 'LineWidth', 2)
xlabel('Time ')
ylabel('Changing of joint variables')
title('Case 1')
legend('\theta_1', '\theta_2', 'd_3', 'Fontsize', 15)
grid on

%Case 2
theta1_d_2 = atan2(py_d, px_d)+pi;
d3_d_2 = sqrt(r_d^2 + s_d^2) - a2n;
theta2_d_2 =pi - atan2(s_d, r_d);

time = 0.1 : 0.1 :10;
theta1_d_time_2 = subs(theta1_d_2, {t}, {time});
theta1_d_t_2 = double(theta1_d_time_2);

theta2_d_time_2 = subs(theta2_d_2, {t}, {time});
theta2_d_t_2 = double(theta2_d_time_2);

d3_d_time_2 = subs(d3_d_2, {t}, {time});
d3_d_t_2 = double(d3_d_time_2);

%Plot
figure(4);
plot(time, [theta1_d_t_2; theta2_d_t_2; d3_d_t_2], 'LineWidth', 2)
xlabel('Time ')
ylabel('Changing of joint variables')
title('Case 2')
legend('\theta_1', '\theta_2', 'd_3', 'Fontsize', 15)
grid on

%% Differential kinematic approach
px_d = 2*a2n*cos(t);
py_d = 2*a2n*cos(2*t);
pz_d = d1*sin(3*t);

dpx_d = -2*a2*sin(t);
dpy_d = -4*a2*sin(2*t);
dpz_d = 3*d1*sin(2*t);
J_v = Jacob1(1:3,:); 
J_inv = simplify(inv(J_v))

v = [dpx_d; dpy_d; dpz_d];

dotq = J_inv*v;
dq = subs(dotq, {a2, d1}, {10, 20})

%%numerical integration
% Euler
theta1_0 = 0;
theta2_0 = 0;
d3_0 = 0;
dq = matlabFunction(dq, 'Vars', {t, [q1;q2;q3]});
q_IK_0 = [theta1_d_t(1); theta2_d_t(1); d3_d_t(1)];
q_IDK = q_IK_0(:, 1);

time = 0 : 0.1 :10;
for i = 1:(length(time) - 1)
    timenow = time(1, i);
    q_IDK(:, i+1) = q_IDK(:, i) + dq(timenow, q_IDK(:, i)) * timenow;
end

figure(5);
plot(time, q_IDK, 'LineWidth', 2)
xlabel('Time ')
ylabel('Changing of joint variables')
title('Euler method')
legend('\theta_1', '\theta_2', 'd_3', 'Fontsize', 15)
grid on

%% ODE45
tspan = [0 : 0.01 : 10];

q0 = [theta1_d_t(1); theta2_d_t(1); d3_d_t(1)]; 

[t,q] = ode45(@(t,q) dq(t,q), tspan, q0);
figure(6);
plot(tspan, q, 'LineWidth', 2)
xlabel('Time ')
ylabel('Changing of joint variables')
title('Runge - Kutta method')
legend('\theta_1', '\theta_2', 'd_3', 'Fontsize', 15)
grid on














