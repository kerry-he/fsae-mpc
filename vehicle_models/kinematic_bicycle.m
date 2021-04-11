function x = kinematic_bicycle(x0, u, dt)
%KINEMATIC_BICYCLE Performs single Euler integration step of the kinematic
%bicycle model
%   INPUTS:
%       x0 - State initial conditions [x; y; theta; v]
%       u - Input controls [a, delta]
%   OUTPUTS:
%       x - Updated state vector

    % Define vehicle constants
    lr = 1.;
    lf = 1.;
    lr_ratio = lr / (lr + lf);
    
    % Calculate derivatives
    beta = atan(lr_ratio * tan(u(2)));
    x_d = x0(4) * cos(x0(3) + beta);
    y_d = x0(4) * sin(x0(3) + beta);
    theta_d = x0(4) / lr * sin(beta);
    v_d = u(1);
    
    % Integrate state
    x = [x0(1) + x_d * dt;
         x0(2) + y_d * dt;
         x0(3) + theta_d * dt;
         x0(4) + v_d * dt];

end
