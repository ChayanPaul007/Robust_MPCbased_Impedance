function Working_Adaptive_MPC_Impedance()
% Adaptive_MPC_Impedance_Fresh
% One fresh MATLAB file, no duplicate function names.
% Includes:
% 1) MPC on x-axis (double integrator) with force-limit induced penetration bound
% 2) Adaptive impedance filter (v1, v3)
% 3) Your full regressor filtering + integral filtering + parameter update law
% 4) Numerics: smooth wall contact + avoid inv(J) + saturations

clc; close all;

% -------------------- User parameters --------------------
p.xe      = 1.0;        % wall location (m)
p.K_env   = 5000;       % wall stiffness (N/m)
p.F_limit = 50;        % wall force limit (N)

% MPC
p.dt      = 0.01;
p.N       = 10;
p.F_margin = 20;        % safety margin in force (N)
p.ax_max   = 8;         % accel bound (m/s^2)
p.Qx       = 200;
p.Ru       = 0.5;
p.Ws       = 1e6;

% Contact smoothing
p.contact_eps = 2e-3;   % meters, smooth max(0,.) for force

% Impedance filter
p.lambda  = 80;
p.alpha   = 35;
p.k_gain  = eye(2);     % your "k" in k\( ... )
p.D       = 65*eye(2);
p.K       = 100*eye(2);

% Parameter adaptation gains (Gamma = diag(gain1))
p.gain1   = [20 10 20]; % you can tune
p.theta_dot_sat = 1500;

% Numerics / conditioning
p.pinv_damp = 1e-3;     % damped J pseudo-inverse
p.tau_max   = 300;      % torque saturation (Nm)
p.pow_hi    = 1.67;
p.pow_lo    = 0.33;

% Geometry
p.l1 = 0.75; p.l2 = 0.75; p.l = p.l1;

% Build MPC optimizer once
fprintf('Designing MPC Optimizer...\n');
p.P_opt = design_mpc_optimizer_fresh(p);

% -------------------- Initial conditions (36 states) --------------------
% State layout (36):
% 1:2   q
% 3:4   q_dot
% 5:7   thetahat (3)
% 8:9   reserved (unused)
% 10:11 Yf2 (2)
% 12:15 h1 (2x2 packed col-wise)
% 16:17 h2 (2)
% 18:19 tauf (2)
% 20:28 YIF (3x3 packed col-wise)
% 29:31 TauIF (3)
% 32:33 v1 (2)
% 34:35 v3 (2)
% 36    ax_filt state (filtered MPC accel)

q0 = [pi/4; pi/2];
dq0 = [0; 0];
thetahat0 = [3.0; 0.1; 0.4];

x0 = zeros(36,1);
x0(1:2) = q0;
x0(3:4) = dq0;
x0(5:7) = thetahat0;
x0(36)  = 0; % filtered ax initial

% Simulate
tspan = 0:p.dt:10;
opts = odeset('RelTol',2e-3,'AbsTol',1e-5,'MaxStep',0.005);
[t,Y] = ode15s(@(t,Y) robot_ode_fresh(t,Y,p), tspan, x0, opts);

% -------------------- Plots --------------------
l1=p.l1; l2=p.l2;
xaxis = l1*cos(Y(:,1)) + l2*cos(Y(:,1)+Y(:,2));
pen   = softplus_fresh(xaxis - p.xe, p.contact_eps);
F     = p.K_env*pen;

figure;
subplot(3,1,1);
plot(t,xaxis,'LineWidth',1.5); hold on;
yline(p.xe,'r--','Wall');
ylabel('x (m)'); title('End-effector x and wall');

subplot(3,1,2);
plot(t,F,'LineWidth',1.5); hold on;
yline(p.F_limit,'k--','Force limit');
ylabel('Force (N)'); title('Wall force');

subplot(3,1,3);
plot(t,Y(:,5),'LineWidth',1.2); hold on;
plot(t,Y(:,6),'LineWidth',1.2);
plot(t,Y(:,7),'LineWidth',1.2);
xlabel('Time (s)'); ylabel('\theta hat');
legend('\theta_1','\theta_2','\theta_3');
title('Parameter estimates');

end

% =====================================================================
% MPC OPTIMIZER (unique name to avoid collisions)
% =====================================================================
function P_opt = design_mpc_optimizer_fresh(p)
N  = p.N;
dt = p.dt;

u = sdpvar(N-1,1);    % accel
x = sdpvar(N,1);      % position
v = sdpvar(N,1);      % velocity
s = sdpvar(N,1);      % slack

x0 = sdpvar(1,1);
v0 = sdpvar(1,1);
xr = sdpvar(N,1);

x_limit = p.xe + (p.F_limit - p.F_margin)/p.K_env;

Constraints = [x(1)==x0, v(1)==v0, s>=0];
for k=1:N-1
    Constraints = [Constraints, ...
        x(k+1) == x(k) + v(k)*dt, ...
        v(k+1) == v(k) + u(k)*dt, ...
        x(k+1) <= x_limit + s(k+1), ...
        -p.ax_max <= u(k) <= p.ax_max];
end

Objective = p.Qx*sum((x-xr).^2) + p.Ru*sum(u.^2) + p.Ws*sum(s.^2);

P_opt = optimizer(Constraints,Objective, ...
    sdpsettings('solver','quadprog','verbose',0), {x0,v0,xr}, u);
end

% =====================================================================
% ODE (unique name)
% =====================================================================
function dYdt = robot_ode_fresh(t, Y, p)

% Unpack
q     = Y(1:2);
qdot  = Y(3:4);
th    = Y(5:7);

% Geometry and trig
l1=p.l1; l2=p.l2; l=p.l;
c1=cos(q(1)); s1=sin(q(1));
c2=cos(q(2)); s2=sin(q(2));
c12=cos(q(1)+q(2)); s12=sin(q(1)+q(2));

% Jacobian
J = [-l1*s1-l2*s12, -l2*s12;
      l1*c1+l2*c12,  l2*c12];

% End-effector position and velocity
x = [l1*c1+l2*c12;
     l1*s1+l2*s12];
xdot = J*qdot;

% Jdot matrix used in your expressions (Jacd)
Jd = [-l1*c1*qdot(1)-l2*c12*(qdot(1)+qdot(2)), -l2*c12*(qdot(1)+qdot(2));
      -l1*s1*qdot(1)-l2*s12*(qdot(1)+qdot(2)), -l2*s12*(qdot(1)+qdot(2))];

% Desired trajectory
xd    = [0.8+0.3*sin(t); 0.8+0.3*cos(t)];
xdotd = [0.3*cos(t);    -0.3*sin(t)];

e    = x - xd;
edot = xdot - xdotd;

% Smooth wall force
pen  = softplus_fresh(x(1) - p.xe, p.contact_eps);
Fext = p.K_env*pen;

% Force limit channel (F1)
F1 = min(Fext, p.F_limit);

% ---------------- MPC x accel ----------------
N=p.N; dt=p.dt;
future_t  = t + (0:N-1)'*dt;
xr_future = 0.8 + 0.3*sin(future_t);
[u_sol, err] = p.P_opt({x(1), xdot(1), xr_future});
if err~=0 || isempty(u_sol)
    ax_raw = 0;
else
    ax_raw = u_sol(1);
end

% Filter ax to prevent violent switching
ax_f = Y(36);
tau_ax = 0.03;
dax_f = (ax_raw - ax_f)/tau_ax;
ax = ax_f;

% PD y accel
yd     = 0.8 + 0.3*cos(t);
ydot_d = -0.3*sin(t);
ay = -0.3*cos(t) - 100*(x(2)-yd) - 20*(xdot(2)-ydot_d);

xdddot = [ax; ay];

% ---------------- Adaptive impedance filter ----------------
v1 = Y(32:33);
v3 = Y(34:35);

dv1 = -p.lambda*v1 - p.k_gain\(p.D*edot + p.K*e + [F1;0]);
dv3 = -p.lambda*v3 + edot;
v2  = edot - p.lambda*v3;
v   = v1 + v2;

% ---------------- Regressor Z ----------------
qr = xdddot - (p.D*edot + p.K*e + [Fext;0]) - Jd*qdot;

% ---------------- Regressor yfilter (2x3) FIXED ----------------
term13 = 2*c2*qr(1) + c2*qr(2) ...
       - s2*Y(3)*Y(4) ...
       - c2*(Y(4)+Y(3))*Y(4) ...
       - (c12*Y(4) - s12*Y(3))*(v(1)/l) ...
       - (c12*(Y(3)+Y(4)))*(v(2)/l);

term23 = c2*qr(1) + s2*(Y(3)^2) ...
       + (c1*Y(4) + c12*Y(4) - s1*Y(3) - s12*Y(3))*(v(1)/l) ...
       + (c12*(Y(3)+Y(4)))*(v(2)/l);

yfilter = [ qr(1),         qr(2),         term13;
            0,     (qr(1)+qr(2)),         term23 ];


% ---------------- Filtered regressor terms ----------------
% y2 and Yf2
y2   = [-s2*(Y(4)^2 + 2*Y(3)*Y(4));  s2*Y(3)^2];
Yf2  = Y(10:11);
dotYf2 = -p.alpha*Yf2 + y2;

% h1 (2x2 packed)
h1 = reshape(Y(12:15),2,2);
f1 = [1 1; 0 1];
doth1 = -p.alpha*h1 + p.alpha*f1*[Y(3) 0; 0 Y(4)];
Yf11  = f1*[Y(3) 0; 0 Y(4)] - h1;

% h2 (2x1)
h2 = Y(16:17);
f2 = [2*c2 c2; c2 0];
f2dot = Y(4)*[-2*s2 -s2; -s2 0];
doth2 = -p.alpha*h2 + (f2dot + p.alpha*f2)*qdot;
Yf22  = f2*qdot - h2;

% final Yf (2x3)
Yf = [Yf11, Yf22 + Yf2];

% ---------------- Estimated dynamics ----------------
hatM = [th(1) + 2*th(3)*c2, th(2) + th(3)*c2;
        th(2) + th(3)*c2,   th(2)];
hatV = [-th(3)*s2*qdot(2), -th(3)*s2*(qdot(1)+qdot(2));
         th(3)*s2*qdot(1), 0];

% Damped J inverse mapping: J \ z is safer than inv(J)*z
tau = hatM*(J\(xdddot - (p.D*edot + p.K*e + [F1;0]) - Jd*qdot)) ...
    + hatV*(qdot - (J\v));

% Torque sat
tau = sat_vec_fresh(tau, -p.tau_max, p.tau_max);

% Torque filtering
tauf = Y(18:19);
dottauf = -p.alpha*tauf + tau;

% ---------------- Integral filtering ----------------
YF   = Yf'*Yf;
TauF = Yf'*tauf;

YIF   = reshape(Y(20:28),3,3);
TauIF = Y(29:31);

dotYIF   = YF;
dotTauIF = TauF;

% ---------------- Parameter update law (your exact structure) ----------------
Gamma = diag(p.gain1);

e1 = tauf - Yf*th;       % 2x1
e2 = TauIF - YIF*th;     % 3x1

thetahatdot = Gamma*Yf'*((abs(e1).^p.pow_hi).*sign(e1)) ...
           + Gamma*Yf'*((abs(e1).^p.pow_lo).*sign(e1)) ...
           + 2*Gamma*((abs(e2).^p.pow_hi).*sign(e2)) ...
           + 2*Gamma*((abs(e2).^p.pow_lo).*sign(e2));

thetahatdot = sat_vec_fresh(thetahatdot, -p.theta_dot_sat, p.theta_dot_sat);

% ---------------- True plant dynamics ----------------
p1 = 3.4; p2 = 0.2; p3 = 0.5;
M  = [p1 + 2*p3*c2, p2 + p3*c2;
      p2 + p3*c2,   p2];
Vm = [-p3*s2*qdot(2), -p3*s2*(qdot(1)+qdot(2));
       p3*s2*qdot(1), 0];

qdd = M \ (tau - Vm*qdot - J'*[Fext;0]);

% ---------------- Pack derivatives ----------------
dYdt = zeros(36,1);
dYdt(1:2)   = qdot;
dYdt(3:4)   = qdd;
dYdt(5:7)   = thetahatdot;

dYdt(10:11) = dotYf2;
dYdt(12:15) = doth1(:);
dYdt(16:17) = doth2;

dYdt(18:19) = dottauf;

dYdt(20:28) = dotYIF(:);
dYdt(29:31) = dotTauIF;

dYdt(32:33) = dv1;
dYdt(34:35) = dv3;

dYdt(36)    = dax_f;

end

% =====================================================================
% Helpers (unique names)
% =====================================================================
function y = softplus_fresh(x, epsval)
% smooth approx of max(0,x)
z = x/epsval;
if z > 50
    y = x;
elseif z < -50
    y = 0;
else
    y = epsval*log1p(exp(z));
end
end

function y = sat_vec_fresh(x, xmin, xmax)
y = min(max(x, xmin), xmax);
end
