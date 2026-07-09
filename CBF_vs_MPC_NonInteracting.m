function CBF_vs_MPC_NonInteracting()
%% CBF_vs_MPC_NonInteracting.m
%
%  COMPARISON: MPC safety vs HOCBF safety on the NON-INTERACTING plane (Y-axis)
%
%  SYSTEM DECOMPOSITION
%  --------------------
%  X-axis (interacting plane)   : adaptive impedance control handles contact
%                                 force at the wall.  IDENTICAL in both controllers.
%  Y-axis (non-interacting plane): the EE must remain inside a workspace box
%                                  [y_min, y_max].  This is where the two
%                                  approaches differ:
%    Proposed : MPC (N-step predictive QP) on Y ensures y in [y_min, y_max]
%    Baseline : HOCBF filter (order 2)    on Y ensures y in [y_min, y_max]
%
%  FAIR COMPARISON DESIGN
%  ----------------------
%  Both controllers use:
%    - Identical 2-DOF planar robot model and true parameters
%    - Identical fixed-time adaptive parameter law (same Gamma, pow_hi, pow_lo)
%    - Identical adaptive impedance structure on X (lambda, D, K, v1/v3 filter)
%    - Identical trajectory xd(t), yd(t) and initial conditions
%    - Identical contact model on X (wall at xe = 1.0 m, K_env = 5000 N/m)
%  Only the Y-axis safety mechanism differs.
%
%  TRAJECTORY: xd = 0.8 + 0.3*sin(t),  yd = 0.8 + 0.3*cos(t)
%  Y-workspace bound: [y_min, y_max] = [0.60, 1.00] m
%  -> yd oscillates between 0.5 and 1.1 m, both bounds are exceeded by the
%     nominal reference, making the safety mechanism active every cycle.
%
%  CBF REFERENCE
%  -------------
%  Singletary A., Kolathaya S., Ames A.D. (2022).
%  "Safety-Critical Kinematic Control of Robotic Systems."
%  IEEE Control Systems Letters, 6, pp. 139-144.
%  DOI: 10.1109/LCSYS.2021.3050609
%
%  HOCBF FOR Y-AXIS (upper bound  h_u = y_max - y):
%    psi_u = -ydot + alpha1*(y_max - y)
%    HOCBF: ay <= -(alpha1+alpha2)*ydot + alpha1*alpha2*(y_max - y)
%
%  HOCBF FOR Y-AXIS (lower bound  h_l = y - y_min):
%    psi_l =  ydot + alpha1*(y - y_min)
%    HOCBF: ay >= -(alpha1+alpha2)*ydot - alpha1*alpha2*(y - y_min)
%
%  Combined: ay_safe = clip(ay_nom, lower_cbf, upper_cbf)

clear; clc; close all;

%% ============================================================
%  SHARED PARAMETERS
%% ============================================================
p.l1 = 0.75;  p.l2 = 0.75;
p.p_true = [3.4; 0.2; 0.5];

% X-axis contact
p.xe        = 1.0;
p.K_env     = 5000;
p.F_limit   = 50;
p.contact_eps = 2e-3;

% Y-axis workspace safety bound
p.y_min = 0.60;     % lower safe limit (m)
p.y_max = 1.00;     % upper safe limit (m)
% Note: yd = 0.8 + 0.3*cos(t) ranges [0.5, 1.1] -> both bounds exceeded

% Simulation
p.dt = 0.01;
p.Tf = 16;
tspan = 0 : p.dt : p.Tf;

% Adaptive impedance (shared, X-axis contact handling)
p.lambda  = 80;
p.alpha   = 35;
p.k_gain  = eye(2);
p.D       = 65*eye(2);
p.K       = 100*eye(2);

% Fixed-time adaptation (IDENTICAL in both controllers)
p.gain1         = [20 10 20];
p.theta_dot_sat = 1500;
p.pow_hi        = 1.67;
p.pow_lo        = 0.33;

% Numerics
p.tau_max   = 300;
p.pinv_damp = 1e-3;

% ------- MPC settings (Y-axis) -------
p.N       = 10;       % prediction horizon (steps)
p.Qy      = 300;      % tracking weight
p.Ru      = 0.5;      % effort weight
p.Ws      = 1e6;      % slack penalty (soft safety)
p.ay_max  = 8;        % accel bound (m/s^2)
p.y_margin = 0.02;    % inner safety margin so constraint is F_limit - margin

% ------- HOCBF settings (Y-axis, Singletary et al. 2022) -------
p.cbf_a1 = 20;        % class-K coefficient level 0->1
p.cbf_a2 = 20;        % class-K coefficient level 1->dot
p.cbf_ay_max = 8;     % same accel bound

%% ============================================================
%  INITIAL CONDITIONS (36-state vector)
%% ============================================================
q0        = [pi/4; pi/2];
dq0       = [0;    0];
thetahat0 = [3.0;  0.1;  0.4];   % off from p_true = [3.4; 0.2; 0.5]

x0 = zeros(36,1);
x0(1:2) = q0;  x0(3:4) = dq0;  x0(5:7) = thetahat0;

%% ============================================================
%  BUILD MPC OPTIMIZER FOR Y-AXIS (once, YALMIP + quadprog)
%% ============================================================
fprintf('=== Building Y-axis MPC optimizer ===\n');
p.Popt_y = build_y_mpc(p);
fprintf('    Done.\n\n');

%% ============================================================
%  SIMULATION A : ADAPTIVE IMPEDANCE (X) + MPC SAFETY (Y)  [Proposed]
%% ============================================================
fprintf('=== Simulating: Adaptive Impedance + MPC on Y (Proposed) ===\n');
tic;
opts = odeset('RelTol',2e-3,'AbsTol',1e-5,'MaxStep',0.005);
[t_mpc, Y_mpc] = ode15s(@(t,Y) ode_mpc(t,Y,p), tspan, x0, opts);
T_mpc = toc;
fprintf('    Done in %.2f s wall-time.\n\n', T_mpc);

%% ============================================================
%  SIMULATION B : ADAPTIVE IMPEDANCE (X) + HOCBF SAFETY (Y) [Baseline]
%% ============================================================
fprintf('=== Simulating: Adaptive Impedance + HOCBF on Y (Baseline) ===\n');
tic;
[t_cbf, Y_cbf] = ode15s(@(t,Y) ode_cbf(t,Y,p), tspan, x0, opts);
T_cbf = toc;
fprintf('    Done in %.2f s wall-time.\n\n', T_cbf);

%% ============================================================
%  POST-PROCESSING
%% ============================================================
ee_mpc = fwd_kin(Y_mpc, p);   % [Nx2]: x,y positions
ee_cbf = fwd_kin(Y_cbf, p);

x_mpc = ee_mpc(:,1);   y_mpc = ee_mpc(:,2);
x_cbf = ee_cbf(:,1);   y_cbf = ee_cbf(:,2);

% Contact forces on X
F_mpc = p.K_env * arrayfun(@(v) sp(v, p.contact_eps), x_mpc - p.xe);
F_cbf = p.K_env * arrayfun(@(v) sp(v, p.contact_eps), x_cbf - p.xe);

% Desired trajectories
yd_mpc = 0.8 + 0.3*cos(t_mpc);
yd_cbf = 0.8 + 0.3*cos(t_cbf);
xd_mpc = 0.8 + 0.3*sin(t_mpc);
xd_cbf = 0.8 + 0.3*sin(t_cbf);

% Y-axis tracking errors
ey_mpc = y_mpc - yd_mpc;
ey_cbf = y_cbf - yd_cbf;

% Y-axis constraint violation: distance outside [y_min, y_max]
viol_mpc = max(0, y_mpc - p.y_max) + max(0, p.y_min - y_mpc);
viol_cbf = max(0, y_cbf - p.y_max) + max(0, p.y_min - y_cbf);

% Number of timesteps violating the bound
n_viol_mpc = sum(y_mpc > p.y_max | y_mpc < p.y_min);
n_viol_cbf = sum(y_cbf > p.y_max | y_cbf < p.y_min);

%% ============================================================
%  METRICS TABLE
%% ============================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('  SAFETY & PERFORMANCE METRICS  (Y-axis non-interacting plane)\n');
fprintf('================================================================\n');
fprintf('%-40s | %10s | %10s\n', 'Metric', 'MPC (Prop.)', 'HOCBF (CBF)');
fprintf('%s\n', repmat('-', 65, 1));
fprintf('%-40s | %10.4f | %10.4f\n', 'Y-axis tracking RMSE (m)',          rms(ey_mpc),         rms(ey_cbf));
fprintf('%-40s | %10.4f | %10.4f\n', 'Max Y constraint violation (m)',     max(viol_mpc),       max(viol_cbf));
fprintf('%-40s | %10.0f | %10.0f\n', 'Steps outside [y_min,y_max]',       n_viol_mpc,          n_viol_cbf);
fprintf('%-40s | %10.4f | %10.4f\n', 'Integrated violation (m*s)',         trapz(t_mpc,viol_mpc),  trapz(t_cbf,viol_cbf));
fprintf('%-40s | %10.4f | %10.4f\n', 'Y std during boundary contact',      std(ey_mpc(viol_mpc>0|[false; diff(viol_mpc(:))~=0])), std(ey_cbf(viol_cbf>0|[false; diff(viol_cbf(:))~=0])));
fprintf('%-40s | %10.4f | %10.4f\n', 'X-tracking RMSE (contact plane, m)',rms(x_mpc-xd_mpc),   rms(x_cbf-xd_cbf));
fprintf('%-40s | %10.2f | %10.2f\n', 'Peak contact force on X (N)',       max(F_mpc),           max(F_cbf));
fprintf('%-40s | %10.3f | %10.3f\n', 'theta1 final est. error',           abs(Y_mpc(end,5)-p.p_true(1)), abs(Y_cbf(end,5)-p.p_true(1)));
fprintf('%-40s | %10.3f | %10.3f\n', 'theta2 final est. error',           abs(Y_mpc(end,6)-p.p_true(2)), abs(Y_cbf(end,6)-p.p_true(2)));
fprintf('%-40s | %10.3f | %10.3f\n', 'theta3 final est. error',           abs(Y_mpc(end,7)-p.p_true(3)), abs(Y_cbf(end,7)-p.p_true(3)));
fprintf('%s\n', repmat('-', 65, 1));
fprintf('%-40s | %10s | %10s\n', 'Safety type',      'Predictive', 'Reactive');
fprintf('%-40s | %10s | %10s\n', 'Y-safety computation', 'N-step QP', 'Closed-form');
fprintf('================================================================\n\n');

%% ============================================================
%  FIGURES
%% ============================================================
c_mpc = [0.10 0.40 0.82];
c_cbf = [0.82 0.25 0.10];
lw = 1.9;

figure('Color','w','Position',[60 50 1180 860], ...
       'Name','CBF vs MPC on Non-Interacting (Y) Plane');

%--- 1. Y position (non-interacting plane - main comparison)
ax1 = subplot(3,2,1);
hold on; grid on;
plot(t_mpc, y_mpc, '-',  'Color',c_mpc, 'LineWidth',lw);
plot(t_cbf, y_cbf, '--', 'Color',c_cbf, 'LineWidth',lw);
plot(t_mpc, yd_mpc, 'k:', 'LineWidth',1.2);
yline(p.y_max,'r-','y_{max}=1.0 m','LineWidth',1.3,'LabelHorizontalAlignment','left');
yline(p.y_min,'b-','y_{min}=0.6 m','LineWidth',1.3,'LabelHorizontalAlignment','left');
ylabel('y_{EE} (m)');
title('Y-axis Position (Non-Interacting Plane) — safety-critical axis');
legend('Adap-MPC (Proposed)','HOCBF (Singletary 2022)','y_d(t)','Location','best');
xlim([0 p.Tf]);

%--- 2. Y constraint violation (positive = outside safe set)
ax2 = subplot(3,2,2);
hold on; grid on;
plot(t_mpc, viol_mpc*1000, '-',  'Color',c_mpc, 'LineWidth',lw);
plot(t_cbf, viol_cbf*1000, '--', 'Color',c_cbf, 'LineWidth',lw);
ylabel('Constraint violation (mm)');
title('Y-axis Safety Violation: dist outside [y_{min}, y_{max}]');
legend('Adap-MPC','HOCBF','Location','best');
xlim([0 p.Tf]);

%--- 3. Y tracking error
ax3 = subplot(3,2,3);
hold on; grid on;
plot(t_mpc, ey_mpc, '-',  'Color',c_mpc, 'LineWidth',lw);
plot(t_cbf, ey_cbf, '--', 'Color',c_cbf, 'LineWidth',lw);
yline(0,'k-','LineWidth',0.8);
ylabel('e_y = y_{EE} - y_d (m)');
title('Y-axis Tracking Error');
legend('Adap-MPC','HOCBF','Location','best');
xlim([0 p.Tf]);

%--- 4. X position and contact force (should match between both — shared impedance)
ax4 = subplot(3,2,4);
yyaxis left;
plot(t_mpc, x_mpc, '-',  'Color',c_mpc, 'LineWidth',lw); hold on;
plot(t_cbf, x_cbf, '--', 'Color',c_cbf, 'LineWidth',lw);
yline(p.xe, 'k-', 'Wall', 'LineWidth',1.2);
ylabel('x_{EE} (m)');
yyaxis right;
plot(t_mpc, F_mpc, '-.', 'Color',c_mpc*0.7, 'LineWidth',1.2);
plot(t_cbf, F_cbf, ':', 'Color',c_cbf*0.7, 'LineWidth',1.2);
ylabel('F_{contact} (N)');
title('X-axis (Interacting Plane) — adaptive impedance, should be similar');
legend('x MPC','x CBF','F MPC','F CBF','Location','best');
xlim([0 p.Tf]);

%--- 5. Parameter estimates
ax5 = subplot(3,2,5);
hold on; grid on;
lh(1) = plot(t_mpc, Y_mpc(:,5), '-',  'Color',c_mpc,      'LineWidth',lw);
lh(2) = plot(t_mpc, Y_mpc(:,6), '-',  'Color',c_mpc*0.6,  'LineWidth',lw);
lh(3) = plot(t_mpc, Y_mpc(:,7), '-',  'Color',c_mpc*0.3,  'LineWidth',lw);
lh(4) = plot(t_cbf, Y_cbf(:,5), '--', 'Color',c_cbf,      'LineWidth',lw);
lh(5) = plot(t_cbf, Y_cbf(:,6), '--', 'Color',c_cbf*0.6,  'LineWidth',lw);
lh(6) = plot(t_cbf, Y_cbf(:,7), '--', 'Color',c_cbf*0.3,  'LineWidth',lw);
yline(p.p_true(1),'k-','LineWidth',0.8);
yline(p.p_true(2),'k-','LineWidth',0.8);
yline(p.p_true(3),'k-','LineWidth',0.8);
xlabel('Time (s)'); ylabel('\hat{\theta}');
title('Parameter Estimates (black = true values, adaptation identical in both)');
legend(lh([1 4]),{'MPC: \theta_{1,2,3}','CBF: \theta_{1,2,3}'},'Location','best');
xlim([0 p.Tf]);

%--- 6. Quantitative bar chart
ax6 = subplot(3,2,6);
hold on; grid on;
bar_labels = {'Y-RMSE \times100 (m)', 'Max viol \times1000 (m)', ...
              'Int. viol \times100', 'Steps viol / 100'};
mpc_b = [rms(ey_mpc)*100,  max(viol_mpc)*1000, ...
         trapz(t_mpc,viol_mpc)*100,  n_viol_mpc/100];
cbf_b = [rms(ey_cbf)*100,  max(viol_cbf)*1000, ...
         trapz(t_cbf,viol_cbf)*100,  n_viol_cbf/100];
bh = bar([mpc_b; cbf_b]', 0.75);
bh(1).FaceColor = c_mpc;
bh(2).FaceColor = c_cbf;
set(gca,'XTickLabel', bar_labels);
xtickangle(20);
ylabel('Scaled metric value');
title('Y-axis Safety Metrics Summary');
legend('Adap-MPC (Proposed)','HOCBF (Singletary 2022)','Location','best');

sgtitle({'MPC vs HOCBF Safety on Non-Interacting (Y) Plane', ...
    'X-axis: adaptive impedance (identical in both)  |  Y-axis: safety mechanism differs'}, ...
    'FontSize',12,'FontWeight','bold');

out_fig = fullfile(fileparts(mfilename('fullpath')), 'fig_CBF_vs_MPC_NonInteracting.png');
saveas(gcf, out_fig);
fprintf('Figure saved: %s\n', out_fig);
end

%% ============================================================
%%  ODE — ADAPTIVE IMPEDANCE (X)  +  MPC SAFETY (Y)  [Proposed]
%% ============================================================
function dYdt = ode_mpc(t, Y, p)

[q, qdot, th] = unpack3(Y);
[J, pos, vel, Jd] = kinematics(q, qdot, p);
[xd, yd, xdotd, ydotd, xddotd, yddotd] = traj(t);
e    = pos - [xd; yd];
edot = vel - [xdotd; ydotd];

% X contact force
Fext = p.K_env * sp(pos(1)-p.xe, p.contact_eps);
F1x  = min(Fext, p.F_limit);

% ---- X: PD impedance (same as CBF controller, no MPC on X) ----
ax = xddotd - 100*(pos(1)-xd) - 20*(vel(1)-xdotd);

% ---- Y: MPC (N-step QP with workspace safety) ----
N  = p.N;  dt = p.dt;
future_t = t + (0:N-1)'*dt;
yr_future = 0.8 + 0.3*cos(future_t);

[u_sol, err] = p.Popt_y({pos(2), vel(2), yr_future});
if err ~= 0 || isempty(u_sol)
    ay = yddotd - 100*(pos(2)-yd) - 20*(vel(2)-ydotd);   % fallback PD
else
    ay = u_sol(1);
end

a_cmd = [ax; ay];

dYdt = core_ode(Y, q, qdot, th, J, pos, vel, Jd, e, edot, Fext, F1x, a_cmd, p);
end

%% ============================================================
%%  ODE — ADAPTIVE IMPEDANCE (X)  +  HOCBF SAFETY (Y)  [CBF Baseline]
%% ============================================================
function dYdt = ode_cbf(t, Y, p)
%
% Y-axis safety via HOCBF order 2 (Singletary et al. 2022):
%
%  Upper bound:  h_u = y_max - y,    HOCBF: ay <= -(a1+a2)*ydot + a1*a2*(y_max-y)
%  Lower bound:  h_l = y - y_min,    HOCBF: ay >= -(a1+a2)*ydot - a1*a2*(y-y_min)
%
%  Combined:  ay_safe = clip(ay_nom, lower_bound, upper_bound)
%  This is a closed-form filter — no QP required at runtime.

[q, qdot, th] = unpack3(Y);
[J, pos, vel, Jd] = kinematics(q, qdot, p);
[xd, yd, xdotd, ydotd, xddotd, yddotd] = traj(t);
e    = pos - [xd; yd];
edot = vel - [xdotd; ydotd];

Fext = p.K_env * sp(pos(1)-p.xe, p.contact_eps);
F1x  = min(Fext, p.F_limit);

% ---- X: PD impedance (identical to MPC controller) ----
ax = xddotd - 100*(pos(1)-xd) - 20*(vel(1)-xdotd);

% ---- Y: Nominal PD reference ----
ay_nom = yddotd - 100*(pos(2)-yd) - 20*(vel(2)-ydotd);

% ---- HOCBF Y safety filter (Singletary et al. 2022, Eq. 12) ----
a1 = p.cbf_a1;  a2 = p.cbf_a2;
y  = pos(2);  ydot = vel(2);

% Upper bound: y <= y_max
cbf_upper = -(a1+a2)*ydot + a1*a2*(p.y_max - y);

% Lower bound: y >= y_min
cbf_lower = -(a1+a2)*ydot - a1*a2*(y - p.y_min);

ay = min(ay_nom, cbf_upper);      % enforce upper bound
ay = max(ay,     cbf_lower);      % enforce lower bound
ay = max(min(ay, p.cbf_ay_max), -p.cbf_ay_max);

a_cmd = [ax; ay];

dYdt = core_ode(Y, q, qdot, th, J, pos, vel, Jd, e, edot, Fext, F1x, a_cmd, p);
end

%% ============================================================
%%  SHARED ODE CORE  (adaptive impedance + fixed-time adaptation)
%%  Identical for both controllers
%% ============================================================
function dYdt = core_ode(Y, q, qdot, th, J, ~, vel, Jd, e, edot, Fext, F1x, a_cmd, p)

c2 = cos(q(2));  s2 = sin(q(2));

% ---- Impedance filter v1, v3 (contact compensation on X) ----
v1 = Y(32:33);  v3 = Y(34:35);
dv1 = -p.lambda*v1 - p.k_gain\(p.D*edot + p.K*e + [F1x; 0]);
dv3 = -p.lambda*v3 + edot;
v_imp = v1 + (edot - p.lambda*v3);

% ---- Regressor filter: h1, h2, Yf2 ----
h1   = reshape(Y(12:15), 2, 2);
f1   = [1 1; 0 1];
Yf11 = f1*[qdot(1) 0; 0 qdot(2)] - h1;
doth1 = -p.alpha*h1 + p.alpha*f1*[qdot(1) 0; 0 qdot(2)];

h2   = Y(16:17);
f2   = [2*c2 c2; c2 0];
f2d  = qdot(2)*[-2*s2 -s2; -s2 0];
Yf22 = f2*qdot - h2;
doth2 = -p.alpha*h2 + (f2d + p.alpha*f2)*qdot;

Yf2    = Y(10:11);
y2     = [-s2*(qdot(2)^2 + 2*qdot(1)*qdot(2));  s2*qdot(1)^2];
dotYf2 = -p.alpha*Yf2 + y2;

Yf = [Yf11,  Yf22 + Yf2];   % 2x3

% ---- Estimated dynamics ----
hatM = [th(1)+2*th(3)*c2,  th(2)+th(3)*c2;
        th(2)+th(3)*c2,    th(2)];
hatV = [-th(3)*s2*qdot(2), -th(3)*s2*(qdot(1)+qdot(2));
         th(3)*s2*qdot(1),  0];

% ---- Control torque ----
tau = hatM*(J\(a_cmd - (p.D*edot + p.K*e + [F1x; 0]) - Jd*qdot)) ...
    + hatV*(qdot - J\v_imp);
tau = min(max(tau, -p.tau_max), p.tau_max);

% ---- Filtered torque ----
tauf    = Y(18:19);
dottauf = -p.alpha*tauf + tau;

% ---- Integral filtering ----
YIF    = reshape(Y(20:28), 3, 3);
TauIF  = Y(29:31);
dotYIF  = Yf'*Yf;
dotTauIF = Yf'*tauf;

% ---- Fixed-time parameter adaptation ----
Gamma = diag(p.gain1);
e1 = tauf  - Yf*th;
e2 = TauIF - YIF*th;
thdot = Gamma*Yf'*( (abs(e1).^p.pow_hi).*sign(e1) ...
                  + (abs(e1).^p.pow_lo).*sign(e1) ) ...
      + 2*Gamma*(   (abs(e2).^p.pow_hi).*sign(e2) ...
                  + (abs(e2).^p.pow_lo).*sign(e2) );
thdot = min(max(thdot, -p.theta_dot_sat), p.theta_dot_sat);

% ---- True plant (uncertain) ----
p1=p.p_true(1); p2=p.p_true(2); p3=p.p_true(3);
M_true = [p1+2*p3*c2, p2+p3*c2; p2+p3*c2, p2];
V_true = [-p3*s2*qdot(2), -p3*s2*(qdot(1)+qdot(2)); p3*s2*qdot(1), 0];
qdd = M_true \ (tau - V_true*qdot - J'*[Fext; 0]);

% ---- Pack ----
dYdt = zeros(36,1);
dYdt(1:2)   = qdot;
dYdt(3:4)   = qdd;
dYdt(5:7)   = thdot;
dYdt(10:11) = dotYf2;
dYdt(12:15) = doth1(:);
dYdt(16:17) = doth2;
dYdt(18:19) = dottauf;
dYdt(20:28) = dotYIF(:);
dYdt(29:31) = dotTauIF;
dYdt(32:33) = dv1;
dYdt(34:35) = dv3;
dYdt(36)    = 0;
end

%% ============================================================
%%  BUILD Y-AXIS MPC OPTIMIZER (YALMIP + quadprog)
%% ============================================================
function Popt = build_y_mpc(p)
N  = p.N;
dt = p.dt;

u  = sdpvar(N-1, 1);   % Y accelerations
y  = sdpvar(N,   1);   % predicted Y positions
vy = sdpvar(N,   1);   % predicted Y velocities
su = sdpvar(N,   1);   % slack for upper bound
sl = sdpvar(N,   1);   % slack for lower bound

y0_v  = sdpvar(1,1);
vy0_v = sdpvar(1,1);
yr    = sdpvar(N,1);   % reference Y trajectory

y_u = p.y_max - p.y_margin;   % inner upper bound (leave margin)
y_l = p.y_min + p.y_margin;   % inner lower bound

C = [y(1)==y0_v, vy(1)==vy0_v, su>=0, sl>=0];
for k = 1:N-1
    C = [C, ...
        y(k+1)  == y(k)  + vy(k)*dt,   ...
        vy(k+1) == vy(k) + u(k)*dt,    ...
        y(k+1)  <= y_u + su(k+1),      ...   % soft upper bound
        y(k+1)  >= y_l - sl(k+1),      ...   % soft lower bound
        -p.ay_max <= u(k) <= p.ay_max];
end

Obj = p.Qy*sum((y-yr).^2) + p.Ru*sum(u.^2) + p.Ws*(sum(su.^2) + sum(sl.^2));

Popt = optimizer(C, Obj, ...
    sdpsettings('solver','quadprog','verbose',0), ...
    {y0_v, vy0_v, yr}, u);
end

%% ============================================================
%%  HELPERS
%% ============================================================
function [J, pos, vel, Jd] = kinematics(q, qdot, p)
l1=p.l1; l2=p.l2;
c1=cos(q(1)); s1=sin(q(1));
c12=cos(q(1)+q(2)); s12=sin(q(1)+q(2));
J  = [-l1*s1-l2*s12, -l2*s12;  l1*c1+l2*c12,  l2*c12];
pos = [l1*cos(q(1))+l2*c12;  l1*sin(q(1))+l2*s12];
vel = J*qdot;
Jd = [-l1*c1*qdot(1)-l2*c12*(qdot(1)+qdot(2)),  -l2*c12*(qdot(1)+qdot(2));
      -l1*s1*qdot(1)-l2*s12*(qdot(1)+qdot(2)),  -l2*s12*(qdot(1)+qdot(2))];
end

function [xd,yd,xdotd,ydotd,xddotd,yddotd] = traj(t)
xd     =  0.8 + 0.3*sin(t);
yd     =  0.8 + 0.3*cos(t);
xdotd  =  0.3*cos(t);
ydotd  = -0.3*sin(t);
xddotd = -0.3*sin(t);
yddotd = -0.3*cos(t);
end

function [q, qdot, th] = unpack3(Y)
q=Y(1:2); qdot=Y(3:4); th=Y(5:7);
end

function pos = fwd_kin(Y, p)
% Compute EE positions from state matrix Y (rows = time)
pos = zeros(size(Y,1), 2);
for i = 1:size(Y,1)
    q = Y(i,1:2);
    pos(i,1) = p.l1*cos(q(1)) + p.l2*cos(q(1)+q(2));
    pos(i,2) = p.l1*sin(q(1)) + p.l2*sin(q(1)+q(2));
end
end

function y = sp(x, eps_val)
% Scalar softplus: smooth max(0, x)
z = x/eps_val;
if z > 50; y = x; elseif z < -50; y = 0;
else;      y = eps_val*log1p(exp(z)); end
end
