function w_opt = create_initial_guess(duration, delta_t, p_S, v_S, p_M, v_M)
%creates an initial guess
N = floor(duration/delta_t);    % number of iterations   

n_x = 4;                % Dimensions of state, here vel & acc.
n_u = 2;                % Dimensions of control
     
AU      = 384403000;    % astronomical unit: distance ersth moon
AT      = 2551400;      % astronomical time: time of one moon circle
    
  G = 6.674 * 10^(-11)*AT^2/AU^3;     % [ m^3 / (kg s^2) ]  Gravitational Constant
M_E = 0;               % [ kg ]              Mass Earth
M_M = 0;                % [ kg ]              Mass Moon

acc_limit = 4 * 10 * (AT^2/AU);
h = 0.001;

%% Optimization problem formulation
import casadi.*


%% formulate and solve the NLP
% Dynamics (Currently using Euler Cromer)
% dynamics = @(x, u, p_M) ode_shuttle(x, u, p_S, p_M);           % ode of the systems
      F = @(x, u, p_M, h) rk4_shuttle(x, u, p_M, h);    % integrator from x_k, u_k to x_k+1

% Start with an empty NLP
w = {};             % decision variables
J = 0;              % cost
g = {};             % constraints

% elimination of initial state -> x0 is not a decision variable
x0bar = [p_S; v_S];
xk = x0bar;
lbg = [];
ubg = [];

% build decision variables, objective, and constraint
for k = 0:N-1
    % New NLP variable for the control u_k
    uk = MX.sym(['u_', num2str(k)], n_u);
    
    % collect in w
    w = {w{:}, uk};            

    % Integrate till the end of the current interval RK4
    %xnext = F(xk, uk, p_M(:, k+1), h);
    
    % Euler-Cromer Step
     pos_next = xk(1:2) + xk(3:4)*delta_t;
     vel_next = xk(3:4) + (uk)*delta_t; 
%     
     xnext = [pos_next; vel_next];
    
    % contribution of stage cost to total cost
    J = J + uk(1)^2 + uk(2)^2;
    
    % TODO include terminal state into objective 
%     if (k == N-1)
%         J = J + norm(xk(1:2) - (p_M(: , end) + [0; -(1737000+200000)/AU]));
%     end
    
    % New NLP variable for state at end of interval
    xk = MX.sym(['x_', num2str(k+1)], n_x);
    
   
    
    % collect in w
    w = {w{:}, xk};
    
    % Add dynamic constraint function,
    % constraining xk to integration result
    g = {g{:},  xk - xnext, uk, norm(xk(1:2)), norm(xk(1:2)-p_M(:,k+1))};
    
    lbg = [lbg; zeros(n_x, 1); -acc_limit; -acc_limit; (6378000 + 35786000)/AU; (1737000)/AU];  %% added constraint so shutle flies not to close to earth and moon
    ubg = [ubg; zeros(n_x, 1); acc_limit; acc_limit; inf; inf];
    
    
end

% Terminal cost

% Terminal condition
% Land "on moon"
g_terminal = [xk(1:2) - (p_M(: , end) + [0; -(1737000)/AU]); ...
              xk(3:4) - v_M(:, end)];

% No further acceleration.
g = {g{:},  g_terminal }; 
lbg = [lbg; zeros(n_x, 1)];
ubg = [ubg; zeros(n_x, 1)];


%% Create an NLP solver
%vertcat({w}) will put all elements of w into a column vector
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% either build lbg and ubg along with g, but as vector, like
% lbg = [];
% lbg = [lbg; ...];
% lbg = [lbg; ...];
% ...
% or just use g in standard form g(x) = 0 and use
% lbg = 0; ubg = 0;

% similar for w0
% for complicated initial guess build along side w, otherwise just
% put it to some value here.
w0 = 1;

% Solve the NLP
sol = solver('x0', w0, 'lbg', lbg, 'ubg', ubg);
%%
close all
% obtain solution
w_opt = full(sol.x);
end

