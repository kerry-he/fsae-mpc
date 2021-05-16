function x = integrate_cart_dyn(x0, u, dt)
%KINEMATIC_BICYCLE Performs single Euler integration step of the kinematic
%bicycle model
%   INPUTS:
%       x0 - State initial conditions [x; y; theta; v]
%       u - Input controls [a, delta]
%       dt - Time step
%   OUTPUTS:
%       x - Updated state vector

    % Perform RK6 integration
    k1 = f_cart_dyn(x0, u);
    k2 = f_cart_dyn(x0 + k1*dt / 2, u);
    k3 = f_cart_dyn(x0 + k1*dt / 4 + k2*dt / 8, u);
    k4 = f_cart_dyn(x0 - k2*dt + 2*k3*dt, u);
    k5 = f_cart_dyn(x0 + 7/27 * k2*dt + 10/27 * k2*dt + k4*dt / 27, u);
    k6 = f_cart_dyn(x0 + 28/625 * k1*dt - k2*dt / 5 + 546/625 * k3*dt + 54/625 * k4*dt - 378/625 * k5*dt, u);

    f = k1/24 + 5/48*k4 + 27/56*k5 + 125/336*k6;
    
    % Integrate state
    x = x0 + dt*f;
     
end

