function x = curvilinear_kinematic_bicycle(x0, u, dt, kappa)
%CURVILINEAR_KINEMATIC_BICYCLE 
%   INPUTS:
%       x0 - State initial conditions [x; y; theta; v]
%       u - Input controls [a, delta]
%       dt - Time step
%       kappa - Curvature profile 
%   OUTPUTS:
%       x - Updated state vector

    f = f_curv_kin(x0, u, kappa);
    x = x0 + f * dt;
     

end

