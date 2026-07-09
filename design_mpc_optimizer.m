function [P_opt] = design_mpc_optimizer(K_env, xe, F_max)
    % N: Prediction Horizon, dt: Sampling time
    N = 20; 
    dt = 0.01;
    
    % Decision Variables: Acceleration (u), Position (x), Velocity (v)
    % We solve for 1D (X-axis) as that is where the wall constraint exists
    u = sdpvar(N, 1); 
    x = sdpvar(N, 1);
    v = sdpvar(N, 1);
    
    % Parameters (Current state and Desired target)
    x0 = sdpvar(1, 1);
    v0 = sdpvar(1, 1);
    xr = sdpvar(N, 1); % Reference trajectory
    
    % Constraints
    Constraints = [x(1) == x0, v(1) == v0];
    for k = 1:N-1
        Constraints = [Constraints, ...
            x(k+1) == x(k) + v(k)*dt + 0.5*u(k)*dt^2, ...
            v(k+1) == v(k) + u(k)*dt];
    end
    
    % FORCE CONSTRAINT: K_env * (x - xe) <= F_max
    % Rearranged: x <= (F_max / K_env) + xe
    x_limit = (F_max / K_env) + xe;
    Constraints = [Constraints, x <= x_limit]; 
    Constraints = [Constraints, -15 <= u <= 15]; % Physical acceleration limits

    % Objective: Minimize tracking error and control effort
    Objective = 100 * norm(x - xr, 2)^2 + 1 * norm(u, 2)^2;
    
    options = sdpsettings('verbose', 0, 'solver', 'quadprog');
    P_opt = optimizer(Constraints, Objective, options, {x0, v0, xr}, u);
end