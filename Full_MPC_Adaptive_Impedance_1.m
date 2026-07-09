function Full_MPC_Adaptive_Impedance()
    %% 1. INITIALIZATION
    clear; clc; close all;
    
    % Constants & Robot Parameters
    p.l = [0.75; 0.75]; 
    p.p_true = [3.4; 0.2; 0.5]; 
    p.xe = 1.5; 
    p.K_env = 500; 
    p.F_limit = 50; % Loosened to 10N for better numerical "breathing room"
    
    % Gains
    p.lambda = 60; % Slightly reduced for stability
    p.k_gain = diag([1, 1]);
    p.D_imp = 65*eye(2); 
    p.K_imp = 100*eye(2);
    p.alpha = 25; 
    p.Gamma = diag([15, 8, 15]); 

    % Design YALMIP MPC Optimizer
    fprintf('Designing Stable MPC Optimizer...\n');
    p.P_opt = design_soft_mpc(p.K_env, p.xe, p.F_limit);

    % Initial Conditions (35 states)
    % Starting at a safe, non-singular configuration
    q_start = [0.45; 1.2]; 
    dq0 = [0; 0];
    thetahat0 = [3.2; 0.15; 0.45]; % Close but not exact estimates
    x0 = [q_start; dq0; thetahat0; zeros(28, 1)];

    %% 2. RUN SIMULATION
    tspan = 0:0.02:10; % Increased step for faster solving
    options = odeset('RelTol', 1e-3, 'AbsTol', 1e-4, 'MaxStep', 0.02);
    
    fprintf('Running simulation... (t=0 to 10)\n');
    [t, Y] = ode15s(@(t, Y) robot_dynamics(t, Y, p), tspan, x0, options);

    %% 3. POST-PROCESSING & PLOTTING
    x_act = p.l(1)*cos(Y(:,1)) + p.l(2)*cos(Y(:,1)+Y(:,2));
    % Calculate force using the same smooth logic as the controller
    dist = x_act - p.xe;
    F_act = p.K_env * (dist + sqrt(dist.^2 + 0.01^2)) / 2; 

    figure('Color', 'w', 'Name', 'Adaptive MPC Impedance Results');
    subplot(2,1,1);
    plot(t, x_act, 'b', 'LineWidth', 2); hold on;
    yline(p.xe, 'r--', 'Wall Position');
    ylabel('X Position (m)'); title('Task Space Tracking'); grid on;

    subplot(2,1,2);
    plot(t, F_act, 'r', 'LineWidth', 2); hold on;
    yline(p.F_limit, 'k--', 'Soft Limit');
    ylabel('Force (N)'); xlabel('Time (s)'); title('Interaction Force'); grid on;
end

%% --- MPC DESIGN (SOFT CONSTRAINTS) ---
function P_opt = design_soft_mpc(K_env, xe, F_max)
    N = 10; dt = 0.02;
    u = sdpvar(N, 1); x = sdpvar(N, 1); v = sdpvar(N, 1); slack = sdpvar(N, 1);
    x0 = sdpvar(1, 1); v0 = sdpvar(1, 1); xr = sdpvar(N, 1);
    x_limit = (F_max / K_env) + xe;
    
    Constraints = [x(1) == x0, v(1) == v0, slack >= 0];
    for k = 1:N-1
        Constraints = [Constraints, ...
            x(k+1) == x(k) + v(k)*dt + 0.5*u(k)*dt^2, ...
            v(k+1) == v(k) + u(k)*dt, ...
            x(k+1) <= x_limit + slack(k)]; 
    end
    % Penalty on slack is high (10^6) to keep it near zero but allow "breathing"
    Objective = 500*norm(x - xr, 2)^2 + 0.1*norm(u, 2)^2 + 1e6*sum(slack);
    P_opt = optimizer(Constraints, Objective, sdpsettings('solver','quadprog','verbose',0), {x0, v0, xr}, u);
end

%% --- ROBOT DYNAMICS & 35-STATE LOGIC ---
function dYdt = robot_dynamics(t, Y, p)
    % State Mapping
    q = Y(1:2); dq = Y(3:4); thetahat = Y(5:7);
    Yf2 = Y(10:11); Yf11 = [Y(12:13), Y(14:15)]; Yf22 = Y(16:17);
    tauf = Y(18:19); YIF = reshape(Y(20:28),3,3); TauIF = Y(29:31);
    v1 = Y(32:33); v3 = Y(34:35);
    
    % Kinematics
    c2 = cos(q(2)); s2 = sin(q(2)); c12 = cos(q(1)+q(2)); s12 = sin(q(1)+q(2));
    Jac = [-p.l(1)*sin(q(1))-p.l(2)*s12, -p.l(2)*s12; p.l(1)*cos(q(1))+p.l(2)*c12, p.l(2)*c12];
    x_pos = [p.l(1)*cos(q(1))+p.l(2)*c12; p.l(1)*sin(q(1))+p.l(2)*s12];
    x_dot = Jac * dq;

    % 1. MPC Command (X-axis)
    N = 10; dt = 0.02;
    xr_future = 0.8 + 0.2*sin(t + (0:N-1)'*dt); % Reduced amplitude for stability
    [u_sol, err] = p.P_opt({x_pos(1), x_dot(1), xr_future});
    ax = u_sol(1); if err ~= 0, ax = 0; end
    
    % Y-axis Tracking
    yd = 0.8 + 0.2*cos(t); yd_dot = -0.2*sin(t);
    ay = -0.2*cos(t) - 100*(x_pos(2)-yd) - 20*(x_dot(2)-yd_dot);
    accel_cmd = [ax; ay];

    % 2. Adaptive Impedance Filter (Smooth Force)
    dist = x_pos(1) - p.xe;
    F_smooth = p.K_env * (dist + sqrt(dist^2 + 0.01^2)) / 2; % Smooth ramp logic
    F_filt = [min(F_smooth, 450); 0]; 
    
    dv1 = -p.lambda*v1 - p.k_gain\(p.D_imp*(x_dot-[0.2*cos(t);-0.2*sin(t)]) + p.K_imp*(x_pos-[0.8+0.2*sin(t);0.8+0.2*cos(t)]) + F_filt);
    dv3 = -p.lambda*v3 + (x_dot-[0.2*cos(t);-0.2*sin(t)]);
    v_corr = v1 + ((x_dot-[0.2*cos(t);-0.2*sin(t)]) - p.lambda*v3);

    % 3. Control Torque
    hatM = [thetahat(1)+2*thetahat(3)*c2, thetahat(2)+thetahat(3)*c2; thetahat(2)+thetahat(3)*c2, thetahat(2)];
    hatVm = [-thetahat(3)*s2*dq(2), -thetahat(3)*s2*(dq(1)+dq(2)); thetahat(3)*s2*dq(1), 0];
    iJac = pinv(Jac); 
    Jacd = [-p.l(1)*cos(q(1))*dq(1)-p.l(2)*c12*sum(dq), -p.l(2)*c12*sum(dq); 
            -p.l(1)*sin(q(1))*dq(1)-p.l(2)*s12*sum(dq), -p.l(2)*s12*sum(dq)];
    
    tau = hatM*iJac*(accel_cmd - v_corr - Jacd*dq) + hatVm*dq + Jac'*[F_smooth; 0];

    % 4. Plant Physics
    M = [p.p_true(1)+2*p.p_true(3)*c2, p.p_true(2)+p.p_true(3)*c2; p.p_true(2)+p.p_true(3)*c2, p.p_true(2)];
    Vm = [-p.p_true(3)*s2*dq(2), -p.p_true(3)*s2*(dq(1)+dq(2)); p.p_true(3)*s2*dq(1), 0];
    qdd = M \ (tau - Vm*dq - Jac'*[F_smooth; 0]);

    % 5. Parameter Estimation Logic (Simplified for stability)
    dotYf2 = -p.alpha*Yf2 + [-s2*(dq(2)^2+2*dq(1)*dq(2)); s2*dq(1)^2];
    doth1 = -p.alpha*Yf11 + p.alpha*[1 1; 0 1]*[dq(1) 0; 0 dq(2)];
    dotYf22 = -p.alpha*Yf22 + (dq(2)*[-2*s2 -s2; -s2 0] + p.alpha*[2*c2 c2; c2 0])*dq;
    
    Yf = [Yf11, Yf22 + Yf2];
    err_a = tauf - Yf*thetahat; err_i = TauIF - YIF*thetahat;
    thetahatdot = p.Gamma*Yf'*(err_a) + 2*p.Gamma*(err_i); % Standard linear for stability
    thetahatdot = max(min(thetahatdot, 100), -100);
    
    % Assemble derivatives
    dYdt = zeros(35, 1);
    dYdt(1:4) = [dq; qdd]; dYdt(5:7) = thetahatdot;
    dYdt(10:11) = dotYf2; dYdt(12:15) = doth1(:); dYdt(16:17) = dotYf22;
    dYdt(18:19) = -p.alpha*tauf + tau; 
    dYF = Yf'*Yf; dTauF = Yf'*tauf;
    dYdt(20:28) = dYF(:); dYdt(29:31) = dTauF;
    dYdt(32:33) = dv1; dYdt(34:35) = dv3;
end