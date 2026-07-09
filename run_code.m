% Init YALMIP MPC
K_env = 5000; xe = 1.0; F_limit = 5; 
P_opt = design_mpc_optimizer(K_env, xe, F_limit);

% Initial conditions
x0 = [pi/4; pi/2; zeros(33,1)]; 
params.K_env = K_env; params.xe = xe; params.lambda = 80;
params.k_gain = eye(2); params.D = 65*eye(2); params.K = 100*eye(2);

% Run
[t, Y] = ode15s(@(t,Y) mpc_adaptive_impedance_ode(t, Y, P_opt, params), [0 10], x0);

% Plot Force
x_actual = 0.75*cos(Y(:,1)) + 0.75*cos(Y(:,1)+Y(:,2));
force = 5000 * max(0, x_actual - 1.0);
plot(t, force); title('MPC Controlled Contact Force'); ylabel('Force (N)');