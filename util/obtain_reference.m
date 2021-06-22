function x_ref = obtain_reference(x, ds, N_s, t, s0, dt, N_t)
%REPARAMETERISE_REFERENCE_S_TO_T Summary of this function goes here
%   Detailed explanation goes here

    L = ds * N_s;

    n = x(1:8:end);
    mu = x(2:8:end);
    x_d     = x(3:8:end);
    y_d     = x(4:8:end);
    theta_d = x(5:8:end);
    delta   = x(6:8:end);
    
    a   = x(7:8:end);
    delta_d = x(8:8:end);
    
    
    idx = zeros(N_t+1, 1);
    rto = zeros(N_t+1, 1);

    idx(1) = floor(mod(s0, L) / ds) + 1;
    rto(1) = mod(mod(s0, L) / ds, 1);
    
    for i = 2:N_t+1        
        t_remaining = dt;
        idx(i) = idx(i - 1);
        rto(i) = rto(i - 1) + t_remaining / t(idx(i));
        t_remaining = t_remaining - t(idx(i - 1)) * (1 - rto(i - 1));
        
        while rto(i) > 1
            idx(i) = nxt(idx(i), N_s);
            rto(i) = t_remaining / t(idx(i));
            t_remaining = t_remaining - t(idx(i));
        end
    end
    

    x_ref = zeros(7, N_t);
    
    for i = 2:N_t+1
        x_ref(1, i-1) = s0 + mod(idx(i) + rto(i) - idx(1) - rto(1), N_s) * ds;
        x_ref(2, i-1) = n(idx(i)) + (n(nxt(idx(i), N_s)) - n(idx(i)))*rto(i);
        x_ref(3, i-1) = mu(idx(i)) + (mu(nxt(idx(i), N_s)) - mu(idx(i)))*rto(i);
        x_ref(4, i-1) = x_d(idx(i)) + (x_d(nxt(idx(i), N_s)) - x_d(idx(i)))*rto(i);
        x_ref(5, i-1) = y_d(idx(i)) + (y_d(nxt(idx(i), N_s)) - y_d(idx(i)))*rto(i);
        x_ref(6, i-1) = theta_d(idx(i)) + (theta_d(nxt(idx(i), N_s)) - theta_d(idx(i)))*rto(i);
        x_ref(7, i-1) = delta(idx(i)) + (delta(nxt(idx(i), N_s)) - delta(idx(i)))*rto(i);
    end
    
    
end

function i_p = prv(i, N)
    i_p = mod(i - 2, N) + 1;
end

function i_n = nxt(i, N)
    i_n = mod(i, N) + 1;
end
