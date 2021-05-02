function [output, status] = pid_controller(target, current, settings, status)
%PID_CONTROLLER Simple implementation of a PID controller
%   Detailed explanation goes here

    kp = settings{1};
    ki = settings{2};
    kd = settings{3};
    max_output = settings{4};
    
    error = target - current;
    integral_error = status{1} + error;
    derivative_error = error - status{2};

    output = kp*error + ki*integral_error + kd*derivative_error;
    output = max(min(output, max_output), -max_output);
    
    status{1} = integral_error;
    status{2} = error;
    
end

