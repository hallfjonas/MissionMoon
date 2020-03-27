function x_next = rk4_moon(x, h)
% one rk4 step
% inputs:
%  x             initial state of integration
%  p_M           position of moon
%  h             time step of integration
% output:
%  x_next        state after one rk4 step

    k1 = ode_moon(x);
    k2 = ode_moon(x+h/2.*k1);
    k3 = ode_moon(x+h/2.*k2);
    k4 = ode_moon(x+h.*k3);
    x_next = x + h/6.*(k1+2*k2+2*k3+k4);
end


