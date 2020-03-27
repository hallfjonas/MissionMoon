function [ dx ] = ode_moon( x )
AU  = 384403000;
AT  = 2360594;
G = 6.674 * 10^(-11)*AT^2/AU^3;
M_E = 5.972 * 10^(24);

p_M = x(1:2);

dx = [ x(3); ...
       x(4); ...
       -G*M_E*p_M(1)/norm(p_M)^3; ...
       -G*M_E*p_M(2)/norm(p_M)^3];
end

