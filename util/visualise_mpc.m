function visualise_mpc(x, x_opt, u_opt, x_P, y_P, dl, dt)
%VISUALISE_MPC Visualise inner workings of MPC
%   INPTUS:
%       x - State of car in Cartesian coordinates
%       x_opt - Optimised state trajectory
%       u_opt - Optimised copntrol trajectory
%       x_P - Spline x coefficients
%       y_P - Spline y coefficients
%       dl - Scale factor (for arclength parameterised spline, otherwise =1)
%       dt - Time step

    persistent f1;
    persistent f2;
    
    if isempty(f1) || isempty(f2)
        f1 = figure(1);
        f2 = figure(2);
    end

    N_x = 7;
    N_steps = length(x_opt) / N_x;
    N_u = round(length(u_opt) / N_steps);

    % Plot trajectory
    set(0, 'CurrentFigure', f1);
    car_marker = plot(x(1), x(2), "ko");
    [x_pred, y_pred, ~] = curvilinear_to_cartesian(x_opt(1:N_x:end), ...
        x_opt(2:N_x:end), x_opt(3:N_x:end), x_P, y_P, dl);
    [x_mid, y_mid, ~] = curvilinear_to_cartesian(x_opt(1:N_x:end), ...
        zeros(N_steps, 1), x_opt(3:N_x:end), x_P, y_P, dl);

%     x_cart_pred = kinematic_bicycle_horizon(x, [u_opt(1:2:N_u*N_steps), ...
%         u_opt(2:2:N_u*N_steps)]', dt);

    car_opt_marker = plot(x_pred, y_pred, "r.");
    mid_opt_marker = plot(x_mid, y_mid, "k-");
%     cart_marker = plot(x_cart_pred(1, :), x_cart_pred(2, :), "b.");

    % Plot controls
    set(0, 'CurrentFigure', f2);
    subplot(2, 1, 1)
    plot(x_opt(4:7:end))
    ylim([0 25])
    ylabel('Velocity (m/s)')
    subplot(2, 1, 2)
    plot(x_opt(7:7:end))
    ylim([-0.4 0.4])
    ylabel('Steering angle (rad)')

    % Animate
    pause(0.1)
    delete(car_marker);
    delete(car_opt_marker);
    delete(mid_opt_marker);
%     delete(cart_marker);

end

