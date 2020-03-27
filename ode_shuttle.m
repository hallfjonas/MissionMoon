function [ dx ] = ode_shuttle( x, u, p_M )
% Gravitational and mass constants of earth.
AU  = 384403000;
AT  = 2360594;
G = 6.674 * 10^(-11)*AT^2/AU^3;
M_E = 5.972 * 10^(24);
M_M = 0;%7.349 * 10^(22);

p_S = x(1:2);

dx = [ x(3); ...
       x(4); ...
       -G*M_E*p_S(1)/norm(p_S)^3 - G*M_M*(p_S(1) - p_M(1))/norm(p_S - p_M)^3 + u(1); ...
       -G*M_E*p_S(2)/norm(p_S)^3 - G*M_M*(p_S(2) - p_M(2))/norm(p_S - p_M)^3 + u(2)];
end
