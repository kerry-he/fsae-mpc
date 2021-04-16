function [u_opt, x_opt] = ltvmpc_kinetmatic_curvilinear(x0, x_ref, kappa, dt, x_lin, u_lin)
%MPC_KINETMATIC_CURVILINEAR Computes a LTV-MPC step for a kinematic bicycle
%model using a curvilinear coordinate frame
%   INPUTS:
%       x0 - Initial state
%       x_ref - Reference trajectory [x_ref_1, x_ref_2, ...]
%       kappa - Spline function
%       dt - Time step
%       x_lin - Linearisation trajectory for state
%       u_lin - Linearisation trajectory for controls
%   OUTPUTS:
%       u_opt - Optimised control trajectory
%       x_opt - Optimised state trajectory

    % Define time horizon
    N_steps = 40;

    % Define constraints
    state_idx = [4, 5];
    soft_idx = [2];
    x_lb = repmat([0, -0.4, -0.75], N_steps, 1);
    x_ub = repmat([1e10, 0.4, 0.75], N_steps, 1);
    x_lb = x_lb(:); x_ub = x_ub(:);

    u_lb = [repmat([-10; -0.4], N_steps, 1); 0];
    u_ub = [repmat([10; 0.4], N_steps, 1); 1e10];

    % Define cost weights
    Q = [5; 50; 10; 0; 0];
    Q_terminal = Q * 10;
    R = [10, 10];
    R_soft = 1e8;
    
    % Define QP problem
    [A, B, d] = linearise_kinematic_curvilinear(x_lin, u_lin, kappa);
    [A_bar, B_bar, d_bar] = sequential_integration(A, B, d, dt);
    [B_bar, xA, lbA, ubA] = state_constraints(A_bar, B_bar, d_bar, x0, x_lb, x_ub, state_idx, soft_idx);
    [H, f] = generate_qp(A_bar, B_bar, d_bar, x0, x_ref, Q, Q_terminal, R, R_soft);
    
    % Solve QP problem
    [u_opt, ~] = qpOASES(H, f, xA, u_lb, u_ub, lbA, ubA);
    x_opt = A_bar*x0 + B_bar*u_opt + d_bar;

end

