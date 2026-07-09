function CBF_vs_AdaptiveMPC_Admittance()
%% CBF_vs_AdaptiveMPC_Admittance.m
%
%  COMPARISON: Fixed-Time Adaptive MPC Admittance Control (Proposed)
%              vs HOCBF-Filtered Admittance Control (CBF Baseline)
%
%  SCENARIO
%  --------
%  A 2-DOF planar robot performs a task under a periodic external force
%  F_ext (simulated human interaction).  An admittance model converts F_ext
%  into a compliant motion x_c so the robot yields to the human naturally.
%  The compliance can push the end-effector outside the safe workspace on
%  the non-interacting (Y) plane.  The safety enforcement mechanism is what
%  differs between the two controllers:
%
%    Proposed  — MPC (N-step predictive QP) on Y keeps y in [y_min, y_max]
%    Baseline  — HOCBF filter (order 2, closed-form) on Y, same bounds
%
%  WHAT IS IDENTICAL IN BOTH CONTROLLERS
%  ---------------------------------------
%   * Admittance model (Ma, Da, Ka) and external force F_ext(t)
%   * Adaptive impedance on X (v1/v3 filter, regressor, integral filter)
%   * Fixed-time parameter update law (Gamma, pow_hi, pow_lo)
%   * Modified reference: x_ref = x_d + x_c,  xdot_ref = xdot_d + xcdot
%   * True robot model and initial conditions / parameter estimates
%   * Sampling time, trajectory, geometry
%
%  WHAT DIFFERS
%  -------------
%   * Y-axis safety: MPC predictive QP  vs  HOCBF reactive filter
%
%  BASE CODE:  Adaptive_ImpedanceX_MPCY_trial.m  (provided by user)
%
%  CBF REFERENCE
%  -------------
%    Singletary A., Kolathaya S., Ames A.D. (2022).
%    "Safety-Critical Kinematic Control of Robotic Systems."
%    IEEE Control Systems Letters, 6, 139-144.
%    DOI: 10.1109/LCSYS.2021.3050609
%
%  ADMITTANCE DYNAMICS
%  -------------------
%    Ma * xc_ddot + Da * xcdot + Ka * xc = F_ext(t)
%    x_ref = x_d + x_c   (admittance-modified reference)
%
%  HOCBF FORMULATION (order 2, Singletary 2022) for Y upper bound h=y_max-y:
%    psi_1  = -ydot + alpha1*(y_max - y)
%    HOCBF: ay <= -(alpha1+alpha2)*ydot + alpha1*alpha2*(y_max - y)
%  For lower bound h = y - y_min:
%    HOCBF: ay >= -(alpha1+alpha2)*ydot - alpha1*alpha2*(y - y_min)
%  Combined: ay_safe = clip(ay_nom, lower_cbf, upper_cbf)
%
%  STATE VECTOR (40 states):
%  1:2   q          3:4  qdot        5:7  thetahat
%  8:9   reserved   10:11 Yf2        12:15 h1(2x2)
%  16:17 h2         18:19 tauf       20:28 YIF(3x3)
%  29:31 TauIF      32:33 v1         34:35 v3
%  36    ay_filt    37:38 xc(2)      39:40 xcdot(2)

% =====================================================================
%  SHARED PARAMETERS  (identical for both simulations)
% =====================================================================

% Geometry
p.l1 = 0.75;  p.l2 = 0.75;  p.l = p.l1;

% Sampling / simulation
p.dt    = 0.01;
p.Tsim  = 15;

% Contact smoothing (no wall, but kept for numerical robustness)
p.contact_eps = 2e-3;

% Impedance filter (same as base code)
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

% Torque saturation
p.tau_max = 300;

% ---- Admittance model ----
% Ma * xc_ddot + Da * xcdot + Ka * xc = F_ext
% Parameters tuned so compliance amplitude ~0.15 m under 8 N at 0.2 Hz
p.Ma = [5;  3];     % virtual mass   [x; y]  (kg)
p.Da = [20; 20];    % virtual damping [x; y]  (N s/m)
p.Ka = [80; 50];    % virtual spring  [x; y]  (N/m)

% ---- External force (simulated human push in Y) ----
% F_ext = [0; F_amp * sin(omega_f * t)]
% omega_f = 0.4*pi ~ 1.26 rad/s (0.2 Hz), slow enough to be human-like
% Admittance compliance amplitude at this freq ~ F_amp / |Z_y|
% |Z_y| at omega=1.26: sqrt((50-3*1.26^2)^2 + (20*1.26)^2) ~ 47 N/m
% -> amplitude ~ 8/47 ~ 0.17 m, so y_ref can reach 1.1+0.17 = 1.27 m
p.F_amp   = 8;        % external force amplitude (N)
p.omega_f = 0.4*pi;   % external force frequency (rad/s)

% ---- Y workspace safety bounds (NON-INTERACTING PLANE) ----
% yd = 0.8 + 0.3*cos(t) ranges [0.5, 1.1]
% With compliance +/-0.17 m, y_ref ranges [0.33, 1.27]
% Safety bounds well-inside -> mechanism needed every cycle
p.y_min = 0.60;     % lower safe limit (m)
p.y_max = 1.00;     % upper safe limit (m)

% ---- Y MPC settings ----
pNy.N      = 40;          % horizon (same as base code)
pNy.dt     = p.dt;
pNy.ay_max = 8;
pNy.y_max  = p.y_max;
pNy.y_min  = p.y_min;
pNy.e_max  = [];          % disabled: admittance ref may exceed bounds
pNy.eT_max = [];          % disabled: same reason
pNy.Qy     = 200;
pNy.Ru     = 0.5;
pNy.Ws     = 1e6;
p.pNy = pNy;

% ---- HOCBF settings ----
p.cbf_a1 = 20;    % class-K coefficient alpha_1
p.cbf_a2 = 20;    % class-K coefficient alpha_2

% =====================================================================
%  BUILD Y-MPC OPTIMIZER ONCE
% =====================================================================
fprintf('=== Building Y-MPC Optimizer (YALMIP + quadprog) ===\n');
p.Py_opt = mpc_y_admittance_opt(pNy);
fprintf('    Done.\n\n');

% =====================================================================
%  INITIAL CONDITIONS  (40 states)
% =====================================================================
q0        = [pi/4; pi/2];
dq0       = [0;    0];
thetahat0 = [3.0;  0.1;  0.4];   % off from true = [3.4; 0.2; 0.5]

x0 = zeros(40, 1);
x0(1:2)  = q0;
x0(3:4)  = dq0;
x0(5:7)  = thetahat0;
% States 37:40 (xc, xcdot) initialised to zero (no initial compliance)

% =====================================================================
%  SIMULATE A: FIXED-TIME ADAPTIVE MPC ADMITTANCE  (Proposed)
% =====================================================================
fprintf('=== Simulating: Fixed-Time Adaptive MPC on Y (Proposed) ===\n');
tspan = 0 : p.dt : p.Tsim;
opts  = odeset('RelTol',2e-3,'AbsTol',1e-5,'MaxStep',0.005);
tic;
[t_mpc, Y_mpc] = ode15s(@(t,Y) ode_adap_mpc(t, Y, p), tspan, x0, opts);
fprintf('    Done in %.2f s wall-time.\n\n', toc);

% =====================================================================
%  SIMULATE B: HOCBF ADMITTANCE  (Singletary et al. 2022 baseline)
% =====================================================================
fprintf('=== Simulating: HOCBF Admittance on Y (CBF Baseline) ===\n');
tic;
[t_cbf, Y_cbf] = ode15s(@(t,Y) ode_adap_cbf(t, Y, p), tspan, x0, opts);
fprintf('    Done in %.2f s wall-time.\n\n', toc);

% =====================================================================
%  POST-PROCESSING
% =====================================================================
l1 = p.l1;  l2 = p.l2;

% End-effector positions
x_ee_mpc = l1*cos(Y_mpc(:,1)) + l2*cos(Y_mpc(:,1)+Y_mpc(:,2));
y_ee_mpc = l1*sin(Y_mpc(:,1)) + l2*sin(Y_mpc(:,1)+Y_mpc(:,2));
x_ee_cbf = l1*cos(Y_cbf(:,1)) + l2*cos(Y_cbf(:,1)+Y_cbf(:,2));
y_ee_cbf = l1*sin(Y_cbf(:,1)) + l2*sin(Y_cbf(:,1)+Y_cbf(:,2));

% Admittance compliance
xc_mpc = Y_mpc(:,37:38);   % [xc_x, xc_y] for MPC run
xc_cbf = Y_cbf(:,37:38);

% Modified Y reference (what the robot tries to track)
yd_nom_mpc = 0.8 + 0.3*cos(t_mpc);
yd_nom_cbf = 0.8 + 0.3*cos(t_cbf);
yref_mpc   = yd_nom_mpc + xc_mpc(:,2);   % admittance-modified reference
yref_cbf   = yd_nom_cbf + xc_cbf(:,2);

% External force trajectory
Fext_y_mpc = p.F_amp * sin(p.omega_f * t_mpc);
Fext_y_cbf = p.F_amp * sin(p.omega_f * t_cbf);

% Y tracking errors (against admittance-modified reference)
ey_mpc = y_ee_mpc - yref_mpc;
ey_cbf = y_ee_cbf - yref_cbf;

% Y safety violations: distance outside [y_min, y_max]
viol_mpc = max(0, y_ee_mpc - p.y_max) + max(0, p.y_min - y_ee_mpc);
viol_cbf = max(0, y_ee_cbf - p.y_max) + max(0, p.y_min - y_ee_cbf);

n_viol_mpc = sum(y_ee_mpc > p.y_max | y_ee_mpc < p.y_min);
n_viol_cbf = sum(y_ee_cbf > p.y_max | y_ee_cbf < p.y_min);

% =====================================================================
%  METRICS TABLE
% =====================================================================
fprintf('\n');
fprintf('=================================================================\n');
fprintf(' Y-AXIS SAFETY & TRACKING METRICS  (Non-Interacting Plane)\n');
fprintf('=================================================================\n');
fprintf('%-42s | %8s | %8s\n', 'Metric', 'MPC(Prop)', 'HOCBF');
fprintf('%s\n', repmat('-',65,1));
fprintf('%-42s | %8.4f | %8.4f\n', 'Y tracking RMSE (m)',               rms(ey_mpc),              rms(ey_cbf));
fprintf('%-42s | %8.4f | %8.4f\n', 'Max Y safety violation (m)',         max(viol_mpc),             max(viol_cbf));
fprintf('%-42s | %8.0f | %8.0f\n', 'Steps outside [y_min, y_max]',      n_viol_mpc,                n_viol_cbf);
fprintf('%-42s | %8.4f | %8.4f\n', 'Integrated violation (m*s)',         trapz(t_mpc,viol_mpc),    trapz(t_cbf,viol_cbf));
fprintf('%-42s | %8.4f | %8.4f\n', 'Max compliance amplitude (m)',       max(abs(xc_mpc(:,2))),    max(abs(xc_cbf(:,2))));
fprintf('%-42s | %8.3f | %8.3f\n', 'theta1 final error',                abs(Y_mpc(end,5)-3.4),    abs(Y_cbf(end,5)-3.4));
fprintf('%-42s | %8.3f | %8.3f\n', 'theta2 final error',                abs(Y_mpc(end,6)-0.2),    abs(Y_cbf(end,6)-0.2));
fprintf('%-42s | %8.3f | %8.3f\n', 'theta3 final error',                abs(Y_mpc(end,7)-0.5),    abs(Y_cbf(end,7)-0.5));
fprintf('%s\n', repmat('-',65,1));
fprintf('%-42s | %8s | %8s\n', 'Y-safety type',  'Predictive',  'Reactive');
fprintf('%-42s | %8s | %8s\n', 'Y-safety computation', 'N-step QP','Closed-fm');
fprintf('=================================================================\n\n');

% =====================================================================
%  PLOTS
% =====================================================================
c_mpc = [0.10  0.40  0.82];
c_cbf = [0.82  0.25  0.10];
lw = 1.9;

figure('Color','w','Position',[60 40 1200 880], ...
    'Name','Fixed-Time Adaptive MPC vs HOCBF: Admittance Safety on Y');

%---- 1. Y end-effector position + safety bounds + admittance reference
ax1 = subplot(3,2,1);
hold on; grid on;
fill([t_mpc(1); t_mpc; t_mpc(end); t_mpc(1)], ...
     [p.y_max; p.y_max*ones(length(t_mpc),1); p.y_min; p.y_min], ...
     [0.9 1.0 0.9], 'FaceAlpha',0.25, 'EdgeColor','none');   % safe region (green band)
plot(t_mpc, yref_mpc, 'k:',  'LineWidth',1.2);    % admittance reference
plot(t_mpc, y_ee_mpc, '-',   'Color',c_mpc, 'LineWidth',lw);
plot(t_cbf, y_ee_cbf, '--',  'Color',c_cbf, 'LineWidth',lw);
yline(p.y_max,'r-','y_{max}','LineWidth',1.3,'LabelHorizontalAlignment','left');
yline(p.y_min,'b-','y_{min}','LineWidth',1.3,'LabelHorizontalAlignment','left');
ylabel('y_{EE}  (m)');
title('Y-axis EE Position (Non-Interacting Plane)');
legend('Safe region','y_d + x_c  (admittance ref)','MPC (Proposed)','HOCBF (Singletary 2022)','Location','best');
xlim([0 p.Tsim]);

%---- 2. External force and admittance compliance
ax2 = subplot(3,2,2);
yyaxis left
plot(t_mpc, Fext_y_mpc, 'k-', 'LineWidth',1.2);
ylabel('F_{ext,y}  (N)');
yyaxis right
plot(t_mpc, xc_mpc(:,2), '-',  'Color',c_mpc, 'LineWidth',lw); hold on;
plot(t_cbf, xc_cbf(:,2), '--', 'Color',c_cbf, 'LineWidth',lw);
ylabel('x_{c,y}  (m)  compliance');
title('External Force and Admittance Compliance in Y');
legend('F_{ext,y}(t)','x_{c,y} MPC','x_{c,y} CBF','Location','best');
grid on; xlim([0 p.Tsim]);

%---- 3. Y safety violation (positive = outside safe set)
ax3 = subplot(3,2,3);
hold on; grid on;
plot(t_mpc, viol_mpc*1000, '-',  'Color',c_mpc, 'LineWidth',lw);
plot(t_cbf, viol_cbf*1000, '--', 'Color',c_cbf, 'LineWidth',lw);
ylabel('Constraint violation  (mm)');
title('Y Safety Violation: dist outside [y_{min}, y_{max}]  (lower = safer)');
legend('Adap-MPC (Proposed)','HOCBF (Singletary 2022)','Location','best');
xlim([0 p.Tsim]);

%---- 4. Y tracking error (against admittance-modified reference)
ax4 = subplot(3,2,4);
hold on; grid on;
plot(t_mpc, ey_mpc, '-',  'Color',c_mpc, 'LineWidth',lw);
plot(t_cbf, ey_cbf, '--', 'Color',c_cbf, 'LineWidth',lw);
yline(0,'k-','LineWidth',0.8);
ylabel('e_y = y_{EE} - y_{ref}  (m)');
title('Y Tracking Error w.r.t. Admittance-Modified Reference');
legend('Adap-MPC','HOCBF','Location','best');
xlim([0 p.Tsim]);

%---- 5. Parameter estimates (adaptation identical - should converge similarly)
ax5 = subplot(3,2,5);
hold on; grid on;
ph1 = plot(t_mpc, Y_mpc(:,5), '-',  'Color',c_mpc,     'LineWidth',lw);
ph2 = plot(t_mpc, Y_mpc(:,6), '-',  'Color',c_mpc*0.6, 'LineWidth',lw);
ph3 = plot(t_mpc, Y_mpc(:,7), '-',  'Color',c_mpc*0.3, 'LineWidth',lw);
ph4 = plot(t_cbf, Y_cbf(:,5), '--', 'Color',c_cbf,     'LineWidth',lw);
ph5 = plot(t_cbf, Y_cbf(:,6), '--', 'Color',c_cbf*0.6, 'LineWidth',lw);
ph6 = plot(t_cbf, Y_cbf(:,7), '--', 'Color',c_cbf*0.3, 'LineWidth',lw);
yline(3.4,'k-','LineWidth',0.8);
yline(0.2,'k-','LineWidth',0.8);
yline(0.5,'k-','LineWidth',0.8);
xlabel('Time (s)'); ylabel('\hat{\theta}');
title('Parameter Estimates (black lines = true values; adaptation identical in both)');
legend([ph1 ph4],{'MPC \theta_{1,2,3}','CBF \theta_{1,2,3}'},'Location','east');
xlim([0 p.Tsim]);

%---- 6. Quantitative bar chart
ax6 = subplot(3,2,6);
hold on; grid on;
bar_vals_mpc = [rms(ey_mpc)*100,    max(viol_mpc)*1000, ...
                trapz(t_mpc,viol_mpc)*10,  n_viol_mpc/10];
bar_vals_cbf = [rms(ey_cbf)*100,    max(viol_cbf)*1000, ...
                trapz(t_cbf,viol_cbf)*10,  n_viol_cbf/10];
bh = bar([bar_vals_mpc; bar_vals_cbf]', 0.75);
bh(1).FaceColor = c_mpc;
bh(2).FaceColor = c_cbf;
set(gca,'XTickLabel', {'RMSE\times100 (m)', 'Max viol\times1000 (m)', ...
    'Int.viol\times10', 'Steps/10'});
xtickangle(15);
ylabel('Scaled metric value');
title('Y-axis Safety Summary');
legend('Adap-MPC (Proposed)','HOCBF (Singletary 2022)','Location','best');

sgtitle({'Fixed-Time Adaptive MPC Admittance vs HOCBF Admittance', ...
    'Non-interacting (Y) plane safety  |  F_{ext,y} = 8 sin(0.4\pi t) N', ...
    'X: adaptive impedance identical in both  |  Adaptation law identical in both'}, ...
    'FontSize',11,'FontWeight','bold');

saveas(gcf, fullfile(fileparts(mfilename('fullpath')), 'fig_Admittance_MPC_vs_CBF.png'));
fprintf('Figure saved: fig_Admittance_MPC_vs_CBF.png\n');
end

%% ===================================================================
%%  ODE A: FIXED-TIME ADAPTIVE MPC ADMITTANCE  (Proposed)
%% ===================================================================
function dYdt = ode_adap_mpc(t, Y, p)

[q, qdot, th, xc, xcdot] = unpack(Y);
[J, x_ee, xdot_ee, Jd]   = jacobian_kin(q, qdot, p);

% Admittance-modified reference
[xd, xdotd, xddotd] = nominal_traj(t);
xref    = xd    + xc;
xdotref = xdotd + xcdot;

% Tracking error against modified reference
e    = x_ee  - xref;
edot = xdot_ee - xdotref;

% Admittance model: Ma*xc_ddot + Da*xcdot + Ka*xc = F_ext
F_ext   = ext_force(t, p);
xc_ddot = (F_ext - p.Da.*xcdot - p.Ka.*xc) ./ p.Ma;

% Impedance filter (no wall force)
v1 = Y(32:33);  v3 = Y(34:35);
dv1 = -p.lambda*v1 - p.k_gain\(p.D*edot + p.K*e);
dv3 = -p.lambda*v3 + edot;
v2  = edot - p.lambda*v3;
v   = v1 + v2;

% ---- X: feedforward of admittance-modified reference (no MPC, no wall) ----
ax = xddotd(1) + xc_ddot(1);

% ---- Y: MPC predictive safety (N-step QP) ----
Ny = p.pNy;
N  = Ny.N;  dt = Ny.dt;
future_t  = t + (0:N-1)' * dt;
yr_future = (0.8 + 0.3*cos(future_t)) + xc(2);   % zero-order hold on compliance

d_mpc = v(2);    % impedance disturbance in Y (same as base code)

[u_sol, err] = p.Py_opt({x_ee(2), xdot_ee(2), yr_future, d_mpc});
if err ~= 0 || isempty(u_sol)
    ay_raw = xddotd(2) + xc_ddot(2);   % fallback: feedforward
else
    ay_raw = u_sol(1);
end

% First-order filter on ay (same as base code, tau = 0.03 s)
ay_f  = Y(36);
day_f = (ay_raw - ay_f) / 0.03;
ay    = ay_f;

a_cmd = [ax; ay];

% Build full derivative vector
dYdt = core_dynamics(Y, t, q, qdot, th, J, x_ee, xdot_ee, Jd, ...
                     e, edot, v, a_cmd, xc_ddot, p);
dYdt(36) = day_f;
end

%% ===================================================================
%%  ODE B: HOCBF ADMITTANCE  (Singletary et al. 2022, baseline)
%% ===================================================================
function dYdt = ode_adap_cbf(t, Y, p)
%
%  CBF SAFETY FILTER (HOCBF order 2, Singletary et al. 2022):
%
%  Upper bound  h_u = y_max - y  (y must stay below y_max):
%    psi_u = -ydot + alpha1*(y_max - y)
%    HOCBF:  ay  <=  -(alpha1+alpha2)*ydot  +  alpha1*alpha2*(y_max - y)
%
%  Lower bound  h_l = y - y_min  (y must stay above y_min):
%    psi_l =  ydot + alpha1*(y - y_min)
%    HOCBF:  ay  >=  -(alpha1+alpha2)*ydot  -  alpha1*alpha2*(y - y_min)
%
%  Closed-form filter (no QP needed for 1D two-sided constraint):
%    ay_safe = clip( ay_nominal,  lower_cbf,  upper_cbf )

[q, qdot, th, xc, xcdot] = unpack(Y);
[J, x_ee, xdot_ee, Jd]   = jacobian_kin(q, qdot, p);

[xd, xdotd, xddotd] = nominal_traj(t);
xref    = xd    + xc;
xdotref = xdotd + xcdot;

e    = x_ee  - xref;
edot = xdot_ee - xdotref;

F_ext   = ext_force(t, p);
xc_ddot = (F_ext - p.Da.*xcdot - p.Ka.*xc) ./ p.Ma;

v1 = Y(32:33);  v3 = Y(34:35);
dv1 = -p.lambda*v1 - p.k_gain\(p.D*edot + p.K*e);
dv3 = -p.lambda*v3 + edot;
v2  = edot - p.lambda*v3;
v   = v1 + v2;

% ---- X: same as MPC controller ----
ax = xddotd(1) + xc_ddot(1);

% ---- Y: nominal feedforward of admittance-modified reference ----
ay_nom = xddotd(2) + xc_ddot(2);

% ---- HOCBF filter (Singletary et al. 2022) ----
a1  = p.cbf_a1;   a2  = p.cbf_a2;
y   = x_ee(2);    ydot = xdot_ee(2);

% Upper bound HOCBF
cbf_upper = -(a1+a2)*ydot + a1*a2*(p.y_max - y);

% Lower bound HOCBF
cbf_lower = -(a1+a2)*ydot - a1*a2*(y - p.y_min);

% Reactive clip
ay = min(ay_nom, cbf_upper);
ay = max(ay,     cbf_lower);
ay = max(min(ay, 8), -8);      % accel saturation (same bound as MPC)

a_cmd = [ax; ay];

dYdt = core_dynamics(Y, t, q, qdot, th, J, x_ee, xdot_ee, Jd, ...
                     e, edot, v, a_cmd, xc_ddot, p);
dYdt(36) = 0;   % ay_filt not used in CBF (no filter state update)
end

%% ===================================================================
%%  CORE DYNAMICS  (shared: adaptive law, regressor, plant)
%%  IDENTICAL for both controllers — only a_cmd and xc_ddot differ
%% ===================================================================
function dYdt = core_dynamics(Y, t, q, qdot, th, J, x_ee, xdot_ee, Jd, ...
                               e, edot, v, a_cmd, xc_ddot, p)

c1=cos(q(1)); s1=sin(q(1));
c2=cos(q(2)); s2=sin(q(2));
c12=cos(q(1)+q(2)); s12=sin(q(1)+q(2));
l  = p.l;

% Regressor reference acceleration
qr = a_cmd - (p.D*edot + p.K*e) - Jd*qdot;

% Filtered regressor terms (same as base code)
y2    = [-s2*(qdot(2)^2 + 2*qdot(1)*qdot(2));  s2*qdot(1)^2];
Yf2   = Y(10:11);
dotYf2 = -p.alpha*Yf2 + y2;

h1 = reshape(Y(12:15),2,2);
f1 = [1 1; 0 1];
doth1 = -p.alpha*h1 + p.alpha*f1*[qdot(1) 0; 0 qdot(2)];
Yf11  = f1*[qdot(1) 0; 0 qdot(2)] - h1;

h2 = Y(16:17);
f2 = [2*c2 c2; c2 0];
f2dot = qdot(2)*[-2*s2 -s2; -s2 0];
doth2 = -p.alpha*h2 + (f2dot + p.alpha*f2)*qdot;
Yf22  = f2*qdot - h2;

Yf = [Yf11, Yf22 + Yf2];   % 2x3

% Estimated dynamics
hatM = [th(1)+2*th(3)*c2,  th(2)+th(3)*c2;
        th(2)+th(3)*c2,    th(2)];
hatV = [-th(3)*s2*qdot(2), -th(3)*s2*(qdot(1)+qdot(2));
         th(3)*s2*qdot(1), 0];

% Control torque (no wall force term)
tau = hatM*(J\(a_cmd - (p.D*edot + p.K*e) - Jd*qdot)) ...
    + hatV*(qdot - (J\v));
tau = sat_vec_fresh(tau, -p.tau_max, p.tau_max);

% Torque filtering
tauf    = Y(18:19);
dottauf = -p.alpha*tauf + tau;

% Integral filtering
YIF    = reshape(Y(20:28),3,3);
TauIF  = Y(29:31);
dotYIF  = Yf'*Yf;
dotTauIF = Yf'*tauf;

% Fixed-time parameter update law (same as base code — IDENTICAL in both)
Gamma = diag(p.gain1);
e1 = tauf  - Yf*th;
e2 = TauIF - YIF*th;
thetahatdot = Gamma*Yf'*((abs(e1).^p.pow_hi).*sign(e1) ...
                        + (abs(e1).^p.pow_lo).*sign(e1)) ...
            + 2*Gamma*  ((abs(e2).^p.pow_hi).*sign(e2) ...
                        + (abs(e2).^p.pow_lo).*sign(e2));
thetahatdot = sat_vec_fresh(thetahatdot, -p.theta_dot_sat, p.theta_dot_sat);

% True plant (external force applied BY human TO robot, + sign in Y)
F_ext    = ext_force(t, p);
p1=3.4; p2=0.2; p3=0.5;
M_true = [p1+2*p3*c2,  p2+p3*c2;
          p2+p3*c2,    p2];
V_true = [-p3*s2*qdot(2), -p3*s2*(qdot(1)+qdot(2));
           p3*s2*qdot(1), 0];
qdd = M_true \ (tau - V_true*qdot + J'*F_ext);   % + F_ext: human pushes robot
t
% Admittance states
v1_state = Y(32:33);  v3_state = Y(34:35);
dv1_new  = -p.lambda*v1_state - p.k_gain\(p.D*edot + p.K*e);
dv3_new  = -p.lambda*v3_state + edot;

% Pack derivatives (40 states)
dYdt = zeros(40,1);
dYdt(1:2)   = qdot;
dYdt(3:4)   = qdd;
dYdt(5:7)   = thetahatdot;
dYdt(10:11) = dotYf2;
dYdt(12:15) = doth1(:);
dYdt(16:17) = doth2;
dYdt(18:19) = dottauf;
dYdt(20:28) = dotYIF(:);
dYdt(29:31) = dotTauIF;
dYdt(32:33) = dv1_new;
dYdt(34:35) = dv3_new;
% state 36 (ay_filt): set by each ODE separately
dYdt(37:38) = Y(39:40);     % d/dt xc  = xcdot
dYdt(39:40) = xc_ddot;      % d/dt xcdot = admittance model
end

%% ===================================================================
%%  Y-MPC OPTIMIZER  (same structure as mpc_y_opt_fresh in base code)
%%  Adds lower bound for two-sided safety
%% ===================================================================
function P = mpc_y_admittance_opt(py)
N  = py.N;
dt = py.dt;

u   = sdpvar(N-1,1);   % ay command
y   = sdpvar(N,1);     % y position
vy  = sdpvar(N,1);     % y velocity
s   = sdpvar(N,1);     % slack

y0_v  = sdpvar(1,1);
vy0_v = sdpvar(1,1);
yr    = sdpvar(N,1);   % admittance-modified reference
d     = sdpvar(1,1);   % disturbance: v(2) from impedance filter

Constraints = [y(1)==y0_v, vy(1)==vy0_v, s>=0];
for k = 1:N-1
    Constraints = [Constraints, ...
        y(k+1)  == y(k)  + vy(k)*dt, ...
        vy(k+1) == vy(k) + (u(k) - d)*dt, ...
        -py.ay_max <= u(k) <= py.ay_max];
end

% Two-sided workspace safety (soft via slack)
Constraints = [Constraints, ...
    y <= py.y_max + s, ...    % upper safe bound
    y >= py.y_min - s];       % lower safe bound

% Tracking error constraint (if set)
if ~isempty(py.e_max)
    Constraints = [Constraints, (y-yr) <= py.e_max + s, -(y-yr) <= py.e_max + s];
end

% Terminal error constraint (if set)
if ~isempty(py.eT_max)
    Constraints = [Constraints, ...
        (y(N)-yr(N)) <= py.eT_max + s(N), ...
        -(y(N)-yr(N)) <= py.eT_max + s(N)];
end

Objective = py.Qy*sum((y-yr).^2) + py.Ru*sum(u.^2) + py.Ws*sum(s.^2);
ops = sdpsettings('solver','quadprog','verbose',0,'cachesolvers',1);
P   = optimizer(Constraints, Objective, ops, {y0_v, vy0_v, yr, d}, u);
end

%% ===================================================================
%%  HELPERS
%% ===================================================================
function [q, qdot, th, xc, xcdot] = unpack(Y)
q     = Y(1:2);
qdot  = Y(3:4);
th    = Y(5:7);
xc    = Y(37:38);
xcdot = Y(39:40);
end

function [J, x_ee, xdot_ee, Jd] = jacobian_kin(q, qdot, p)
l1=p.l1; l2=p.l2;
c1=cos(q(1)); s1=sin(q(1));
c12=cos(q(1)+q(2)); s12=sin(q(1)+q(2));
J  = [-l1*s1-l2*s12, -l2*s12;
       l1*cos(q(1))+l2*c12,  l2*c12];
x_ee   = [l1*cos(q(1))+l2*c12;  l1*sin(q(1))+l2*s12];
xdot_ee = J*qdot;
Jd = [-l1*cos(q(1))*qdot(1)-l2*c12*(qdot(1)+qdot(2)),  -l2*c12*(qdot(1)+qdot(2));
      -l1*s1*qdot(1)-l2*s12*(qdot(1)+qdot(2)),          -l2*s12*(qdot(1)+qdot(2))];
end

function [xd, xdotd, xddotd] = nominal_traj(t)
xd     = [0.8 + 0.3*sin(t);   0.8 + 0.3*cos(t)];
xdotd  = [0.3*cos(t);        -0.3*sin(t)];
xddotd = [-0.3*sin(t);       -0.3*cos(t)];
end

function F = ext_force(t, p)
% External force applied BY human TO robot (task-space, [x; y])
% No force in X (non-contact direction assumed free)
% Periodic push in Y simulates human interaction during task
F = [0;  p.F_amp * sin(p.omega_f * t)];
end

function y = softplus_fresh(x, epsval)
z = x/epsval;
if z > 50;    y = x;
elseif z < -50; y = 0;
else;          y = epsval*log1p(exp(z));
end
end

function y = sat_vec_fresh(x, xmin, xmax)
y = min(max(x, xmin), xmax);
end
