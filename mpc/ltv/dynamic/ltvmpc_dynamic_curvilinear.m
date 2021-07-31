function [u_opt, x_opt, QP, exitflag, fval, slack_opt] = ltvmpc_dynamic_curvilinear(x0, x_ref, kappa, dt, x_lin, u_lin, QP)
%MPC_KINETMATIC_CURVILINEAR Computes a LTV-MPC step for a kinematic bicycle
%model using a curvilinear coordinate frame
%   INPUTS:
%       x0 - Initial state
%       x_ref - Reference trajectory [x_ref_1, x_ref_2, ...]
%       kappa - Curvature profile
%       kappa_d - Derivative of curvature profile
%       dt - Time step
%       x_lin - Linearisation trajectory for state
%       u_lin - Linearisation trajectory for controls
%   OUTPUTS:
%       u_opt - Optimised control trajectory
%       x_opt - Optimised state trajectory

    % Define time horizon
    N_steps = length(x_ref);

    % Define constraints
    state_idx = [4, 7];
    soft_idx = [2];
    
    x_lb = repmat([0, -0.4, -0.75], N_steps, 1);    
    x_ub = repmat([inf, 0.4, 0.75], N_steps, 1);

    x_lb = x_lb(:); x_ub = x_ub(:);

    u_lb = [repmat([0; -0.4], N_steps, 1); 0; 0; 0; 0];
    u_ub = [repmat([1e19; 0.4], N_steps, 1); 1e19; 1e19; 1e19; 1e19];

    % Define cost weights
    Q = [5; 250; 10; 10; 10; 10; 10];
    Q_terminal = Q * 10;
    R = [0.1, 500];
    R_soft = [1e7; 1e6; 1e6; 1e5];
        
    % Define QP problem
    [A, B, d] = rk4_dynamic_curvilinear(x_lin, u_lin, kappa, dt);
    [A_bar, B_bar, d_bar] = sequential_integration(A, B, d, dt);
    [B_bar, xA, lbA, ubA] = dynamic_state_constraints(A_bar, B_bar, d_bar, x0, x_lb, x_ub, state_idx, soft_idx, x_lin, u_lin, kappa);
    [H, f, const] = generate_qp(A_bar, B_bar, d_bar, x0, x_ref, Q, Q_terminal, R, R_soft);
    
    % Solve QP problem
%     options = qpOASES_options('MPC');
    
%     if QP == 0
%         [QP, u_opt, fval, exitflag, iter, lambda, auxOutput] = qpOASES_sequence('i', H, f, xA, u_lb, u_ub, lbA, ubA, options);
%     else
%         [u_opt, fval, exitflag, iter, lambda, auxOutput] = qpOASES_sequence('m', QP, H, f, xA, u_lb, u_ub, lbA, ubA, options);
%     end
    
    [u_opt, fval, exitflag, iter, lambda, auxOutput] = qpOASES(H, f, xA, u_lb, u_ub, lbA, ubA);
    if exitflag
        display(exitflag)
    end
    
    slack_opt = u_opt(end-length(R_soft)+1:end);
    x_opt = A_bar*x0 + B_bar*u_opt + d_bar;
    u_opt = u_opt(1:end-length(R_soft));
    fval = fval + const;

end

