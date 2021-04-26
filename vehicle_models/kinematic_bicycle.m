function x = kinematic_bicycle(x0, u, dt)
%KINEMATIC_BICYCLE Performs single Euler integration step of the kinematic
%bicycle model
%   INPUTS:
%       x0 - State initial conditions [x; y; theta; v]
%       u - Input controls [a, delta]
%       dt - Time step
%   OUTPUTS:
%       x - Updated state vector

    % Define vehicle constants
    k1 = f_cart_kin(x0, u);
    k2 = f_cart_kin(x0 + k1*dt / 2, u);
    k3 = f_cart_kin(x0 + k2*dt / 2, u);
    k4 = f_cart_kin(x0 + k3*dt, u);

    f = (k1 + 2*k2 + 2*k3 + k4) / 6;
    
    % Integrate state
    x = x0 + dt*f;
     
end

