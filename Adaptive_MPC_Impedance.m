function Adaptive_MPC_Impedance()
    % Full single-file script: Adaptive impedance filter + MPC force limiting
    % Requires: YALMIP + quadprog

    clear; clc; close all;

    % ------------------ Parameters ------------------
    params.xe = 1.0;          % wall location in x
    params.K_env = 500;       % wall stiffness
    params.F_limit = 50;      % max allowable wall force

    % MPC settings
    params.dt = 0.01;
    params.N  = 12;
    params.F_margin = 5;      % robustness margin (N)
    params.ax_max   = 5.0;    % bound on x-accel (m/s^2)

    % Smooth contact settings (key for ode15s stability)
    params.contact_eps = 2e-3;   % smoothing length in meters (try 1e-3 to 5e-3)

    % Low pass filter for MPC acceleration (key for smoothness)
    params.ax_filt_tau = 0.03;   % seconds (20 to 50 ms typical)

    % Impedance filter settings
    params.lambda = 80;
    params.k_gain = diag([1, 1]);
    params.K0 = diag([100, 100]);
    params.D0 = diag([65,  65]);

    % Adaptive impedance bounds and rates (adapting x only)
    params.Kx_min = 20;   params.Kx_max = 200;
    params.Dx_min = 20;   params.Dx_max = 160;
    params.gammaK = 8;    % stiffness reduction rate near force limit
    params.gammaD = 10;   % damping increase rate near force limit
    params.epsF   = 2;    % start adapting within epsF of limit

    % Damped pseudoinverse and torque clamp
    params.pinv_damp = 1e-3;      % damping for J^T(JJ^T + mu I)^{-1}
    params.tau_max   = 150;       % torque saturation (Nm) adjust if needed

    % ------------------ Build MPC Optimizer ------------------
    fprintf('Designing MPC Optimizer...\n');
    params.P_opt = design_mpc_optimizer(params);

    % ------------------ Initial Conditions (36 states now) ------------------
    % State layout:
    % Y(1:2)=q, Y(3:4)=dq, Y(5:7)=thetahat,
    % Y(8)=Kx, Y(9)=Dx,
    % Y(10)=ax_filt (filtered x accel state),
    % Y(11:32)=unused (22 states),
    % Y(33:34)=v1(2), Y(35:36)=v3(2)

    q0 = [pi/4; pi/2];
    dq0 = [0; 0];
    thetahat0 = [3.0; 0.1; 0.4];

    rest0 = zeros(29,1);          % was 28, now 29 because we add ax_filt
    rest0(1) = params.K0(1,1);    % Y(8)  = Kx
    rest0(2) = params.D0(1,1);    % Y(9)  = Dx
    rest0(3) = 0;                 % Y(10) = ax_filt initial

    x0 = [q0; dq0; thetahat0; rest0];

    % ------------------ Simulate ------------------
    tspan = 0:params.dt:10;

    % Relax tolerances a bit and limit stiffness shock
    options = odeset('RelTol',2e-3,'AbsTol',1e-5,'MaxStep',0.005);
    [t, Y] = ode15s(@(t,Y) robot_ode(t,Y,params), tspan, x0, options);

    % ------------------ Post-processing ------------------
    l1 = 0.75; l2 = 0.75;
    xaxis = l1*cos(Y(:,1)) + l2*cos(Y(:,1)+Y(:,2));

    % use the same smooth contact model for plotting
    pen = softplus(xaxis - params.xe, params.contact_eps);
    F_actual = params.K_env * pen;

    figure;
    subplot(3,1,1);
    plot(t, xaxis, 'LineWidth',1.5); hold on;
    yline(params.xe,'r--','Wall');
    ylabel('X (m)'); title('Task-space X with MPC');

    subplot(3,1,2);
    plot(t, F_actual, 'LineWidth',1.5); hold on;
    yline(params.F_limit,'k--','Force limit');
    ylabel('Force (N)'); title('Wall force');

    subplot(3,1,3);
    plot(t, Y(:,8), 'LineWidth',1.5); hold on;
    plot(t, Y(:,9), 'LineWidth',1.5);
    ylabel('Adaptive gains'); xlabel('Time (s)');
    legend('Kx','Dx'); title('Adaptive impedance gains (x direction)');
end

% =====================================================================
% MPC OPTIMIZER
% =====================================================================
function P_opt = design_mpc_optimizer(p)
    N  = p.N;
    dt = p.dt;

    % Decision variables
    u = sdpvar(N-1,1);     % x-accel
    x = sdpvar(N,1);       % x position
    v = sdpvar(N,1);       % x velocity
    s = sdpvar(N,1);       % slack for soft constraint

    % Parameters to optimizer
    x0 = sdpvar(1,1);
    v0 = sdpvar(1,1);
    xr = sdpvar(N,1);

    % Robust penetration bound
    x_limit = p.xe + (p.F_limit - p.F_margin)/p.K_env;

    Constraints = [x(1)==x0, v(1)==v0, s>=0];

    for k = 1:N-1
        Constraints = [Constraints, ...
            x(k+1) == x(k) + v(k)*dt, ...
            v(k+1) == v(k) + u(k)*dt];

        % Soft constraint: x <= x_limit + slack
        Constraints = [Constraints, x(k+1) <= x_limit + s(k+1)];

        % Input bounds
        Constraints = [Constraints, -p.ax_max <= u(k) <= p.ax_max];
    end

    % Objective
    Qx = 20;
    Ru = 0.5;
    Ws = 1e6;
    Objective = Qx*sum((x - xr).^2) + Ru*sum(u.^2) + Ws*sum(s.^2);

    P_opt = optimizer(Constraints, Objective, ...
        sdpsettings('solver','quadprog','verbose',0), {x0,v0,xr}, u);
end

% =====================================================================
% MAIN ODE
% =====================================================================
function dYdt = robot_ode(t, Y, p)
    % Extract states
    q = Y(1:2);
    dq = Y(3:4);
    th = Y(5:7);

    % Adaptive impedance gains (x only)
    Kx = Y(8);
    Dx = Y(9);

    % Filtered x accel state
    ax_f = Y(10);

    % Filter states
    v1 = Y(33:34);
    v3 = Y(35:36);

    % Kinematics
    l1 = 0.75; l2 = 0.75;
    c1 = cos(q(1)); s1 = sin(q(1));
    c2 = cos(q(2)); s2 = sin(q(2));
    c12 = cos(q(1)+q(2)); s12 = sin(q(1)+q(2));

    Jac = [-l1*s1-l2*s12, -l2*s12;
            l1*c1+l2*c12,  l2*c12];

    x_pos = [l1*c1+l2*c12;
             l1*s1+l2*s12];

    x_dot = Jac*dq;

    % Smooth wall force
    pen = softplus(x_pos(1) - p.xe, p.contact_eps);
    F_actual = p.K_env * pen;

    % Desired trajectory
    xd_vec = [0.8+0.3*sin(t); 0.8+0.3*cos(t)];
    xdot_d_vec = [0.3*cos(t); -0.3*sin(t)];

    e = x_pos - xd_vec;
    edot = x_dot - xdot_d_vec;

    % ---------------- MPC for raw x accel ----------------
    N  = p.N;
    dt = p.dt;
    future_t = t + (0:N-1)'*dt;
    xr_future = 0.8 + 0.3*sin(future_t);

    [u_sol, err] = p.P_opt({x_pos(1), x_dot(1), xr_future});
    if err ~= 0 || isempty(u_sol)
        ax_raw = 0;
    else
        ax_raw = u_sol(1);
    end

    % Low-pass filter the MPC accel: dax = (ax_raw - ax_f)/tau
    dax_f = (ax_raw - ax_f)/p.ax_filt_tau;
    ax = ax_f;

    % ---------------- PD for y accel ----------------
    yd = 0.8 + 0.3*cos(t);
    ydot_d = -0.3*sin(t);
    ay = -0.3*cos(t) - 100*(x_pos(2)-yd) - 20*(x_dot(2)-ydot_d);

    accel_cmd = [ax; ay];

    % ---------------- Adaptive impedance update (x only) ----------------
    sigma = max(0, F_actual - (p.F_limit - p.epsF));
    dKx_raw = -p.gammaK * sigma;
    dDx_raw =  p.gammaD * sigma;

    dKx = proj_scalar(Kx, dKx_raw, p.Kx_min, p.Kx_max);
    dDx = proj_scalar(Dx, dDx_raw, p.Dx_min, p.Dx_max);

    K_ad = p.K0;  D_ad = p.D0;
    K_ad(1,1) = Kx;
    D_ad(1,1) = Dx;

    % ---------------- Adaptive impedance filter ----------------
    F_vec = [min(F_actual, 450); 0];

    dv1 = -p.lambda*v1 - (p.k_gain \ (D_ad*edot + K_ad*e + F_vec));
    dv3 = -p.lambda*v3 + edot;
    v   = v1 + (edot - p.lambda*v3);

    % ---------------- Control torque (computed torque style) ----------------
    % Damped pseudoinverse for better conditioning: J^T (J J^T + mu I)^{-1}
    JJt = Jac*Jac';
    iJJt = (JJt + p.pinv_damp*eye(2)) \ eye(2);
    iJac = Jac' * iJJt;

    Jacd = [-l1*c1*dq(1)-l2*c12*(dq(1)+dq(2)), -l2*c12*(dq(1)+dq(2));
            -l1*s1*dq(1)-l2*s12*(dq(1)+dq(2)), -l2*s12*(dq(1)+dq(2))];

    % Estimated model
    hatM = [th(1) + 2*th(3)*c2, th(2) + th(3)*c2;
            th(2) + th(3)*c2,   th(2)];
    hatVm = [-th(3)*s2*dq(2), -th(3)*s2*(dq(1)+dq(2));
              th(3)*s2*dq(1), 0];

    tau = hatM * iJac * (accel_cmd - v - Jacd*dq) + hatVm*dq + Jac'*[F_actual; 0];

    % Torque saturation to avoid violent spikes
    tau = max(min(tau, p.tau_max), -p.tau_max);

    % ---------------- True plant dynamics ----------------
    p1 = 3.4; p2 = 0.2; p3 = 0.5;
    M = [p1 + 2*p3*c2, p2 + p3*c2;
         p2 + p3*c2,   p2];
    Vm = [-p3*s2*dq(2), -p3*s2*(dq(1)+dq(2));
           p3*s2*dq(1), 0];

    qdd = M \ (tau - Vm*dq - Jac'*[F_actual; 0]);

    % Parameter adaptation placeholder
    thetahatdot = zeros(3,1);

    % ---------------- Assemble derivatives ----------------
    dYdt = zeros(36,1);
    dYdt(1:2)   = dq;
    dYdt(3:4)   = qdd;
    dYdt(5:7)   = thetahatdot;

    dYdt(8)     = dKx;
    dYdt(9)     = dDx;
    dYdt(10)    = dax_f;

    dYdt(33:34) = dv1;
    dYdt(35:36) = dv3;
end

% =====================================================================
% Helpers
% =====================================================================
function y = softplus(x, epsval)
    % Smooth approximation of max(0,x):
    % y = eps * log(1 + exp(x/eps))
    % stable for large negative/positive x
    z = x/epsval;
    if z > 50
        y = x;                 % exp huge => max approx x
    elseif z < -50
        y = 0;                 % exp tiny
    else
        y = epsval*log1p(exp(z));
    end
end

function dx = proj_scalar(x, dx_raw, xmin, xmax)
    if (x <= xmin && dx_raw < 0)
        dx = 0;
    elseif (x >= xmax && dx_raw > 0)
        dx = 0;
    else
        dx = dx_raw;
    end
end

