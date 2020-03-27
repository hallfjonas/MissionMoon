function x_next = rk4_shuttle(x, u, p_M, h)
% one rk4 step
% inputs:
%  x             initial state of integration
%  u             control, kept constant over integration
%  p_S, p_M      position of shuttle and moon, respectively
%  h             time step of integration
% output:
%  x_next        state after one rk4 step

    k1 = ode_shuttle(x, u, p_M);
    k2 = ode_shuttle(x+h/2.*k1, u, p_M);
    k3 = ode_shuttle(x+h/2.*k2, u, p_M);
    k4 = ode_shuttle(x+h.*k3, u, p_M);
    x_next = x + h/6.*(k1+2*k2+2*k3+k4);
end


