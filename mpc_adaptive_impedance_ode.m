function dYdt = mpc_adaptive_impedance_ode(t, Y, P_opt, params)
    % Extract States
    q = Y(1:2); dq = Y(3:4);
    thetahat = Y(5:7);
    v1 = Y(32:33); v3 = Y(34:35);
    
    % 1. Kinematics
    [x_pos, J] = forward_kin(q);
    x_dot = J * dq;
    
    % 2. MPC Outer Loop (X-axis only for force limiting)
    % Predict reference for the horizon
    dt = 0.01; N = 20;
    future_t = t + (0:N-1)'*dt;
    xr_future = 0.8 + 0.3*sin(future_t);
    
    % Solve MPC for optimal X-acceleration
    [u_mpc_x, error_code] = P_opt({x_pos(1), x_dot(1), xr_future});
    if error_code ~= 0, ax_cmd = 0; else, ax_cmd = u_mpc_x(1); end
    
    % Y-axis remains pure tracking (or another MPC)
    ay_cmd = -0.3*cos(t) - 20*(x_pos(2) - (0.8+0.3*cos(t))) - 10*(x_dot(2) - (-0.3*sin(t)));
    accel_task = [ax_cmd; ay_cmd];

    % 3. Force and Adaptive Filter
    F_raw = params.K_env * max(0, x_pos(1) - params.xe);
    F1 = min(F_raw, 450); % Filter still sees capped force
    F_vec = [F1; 0];
    
    % Your original dv1, dv3 equations
    e = x_pos - [0.8+0.3*sin(t); 0.8+0.3*cos(t)];
    edot = x_dot - [0.3*cos(t); -0.3*sin(t)];
    dv1 = -params.lambda*v1 - params.k_gain\(params.D*edot + params.K*e + F_vec);
    dv3 = -params.lambda*v3 + edot;
    v = v1 + (edot - params.lambda*v3);

    % 4. Feedback Linearization (Inner Loop)
    % Linearize system to double integrator using estimates
    [hatM, hatVm] = get_estimated_dynamics(q, dq, thetahat);
    Jacd = get_jac_dot(q, dq);
    iJac = pinv(Jac);
    
    % Command Torque
    tau = hatM * iJac * (accel_task - v - Jacd*dq) + hatVm*dq + Jac'*F_vec;

    % 5. Real Robot Physics
    [M, Vm] = get_real_dynamics(q, dq);
    qdd = M \ (tau - Vm*dq - Jac'*F_vec);

    % ... (thetahatdot and other filter derivatives same as previous code) ...
    dYdt = [dq; qdd; thetahatdot; ... (rest of states) ... ; dv1; dv3];
end