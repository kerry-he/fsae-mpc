function x_traj = kinematic_bicycle_horizon(x0, u_traj, dt)
%KINETMATIC_BICYCLE_HORIZON Sequentially integrates a state forwards given
%a control trajectory
%   INPTUS:
%       x0 - Initial state
%       u_traj - Control trajectory
%       dt - Time step
%   OUTPUTS:
%       x_traj - State trajectory

    N_steps = length(u_traj);
    
    x_traj = zeros(length(x0), N_steps + 1);
    x_traj(:, 1) = x0;
    
    for i = 1:N_steps
        x_traj(:, i + 1) = kinematic_bicycle(x_traj(:, i), u_traj(:, i), dt);
    end
    
end

