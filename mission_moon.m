% Main script for solivng shuttle to moon mission
clc; clear variables; close all;
import casadi.*

%% Initialization
duration = 2/3;                 % [AT] duration of space travel
delta_t = 0.005;                % [AT] discretization times
h = delta_t;                    % Time step in RK4 integrator
N = floor(duration/delta_t);    % number of iterations

% Dimensions
n_x = 4;                % Dimensions of state, here vel & acc.
n_u = 2;                % Dimensions of control
n_p = 2;                % Dimensiion of parameters
% Astronomical units and correction variables
       AU = 384403000;    % astronomical unit: distance ersth moon
       AT = 2551400;      % astronomical time: time of one moon circle
v_correct = AT/AU;        % Correction term to change between units
a_correct = AT^2/AU;      % Correction term to change between units
        G = 6.674 * 10^(-11)*AT^2/AU^3;     % Gravitational Constant
      M_E = 5.972 * 10^(24);                % Mass Earth
      M_M = 7.347 * 10^(22);                % Mass Moon

% Starting positions and velocities (of moon and shuttle)
p_M = zeros(2,N+1);
p_M(:,1) = [1; 0];                    % [ 1 unit = 1 AU ]
p_S = [-(6378000 + 35786000)/AU; 0];   %3000006378000/AU  ;             % [m] geostationary hight and velocity
%p_S = [0; 1/4];
v_M = [0; 2*pi];                      % [ AU / AT ]
v_theo = sqrt(G*M_E/norm(p_S));
v_S = [0; -v_theo ]; %7e-04 -3100 * (AT/AU)

% Some constraint values
acc_limit = 4 * 10 * (AT^2/AU);       % Acceleration constraint to shuttle
consumption_limit = inf;                % Maximal consumption limit
alpha = 0;                              % factor for final cost
epsilon_pos_start = 0.1; %0.1                  % Radius of final distance
epsilon_vel_start = 1; %1              % Tolerance to speed diff to moon
number_of_casadi_Iter = 4;
epsilon_vel_pos_start = [epsilon_pos_start; epsilon_vel_start];

%% Dynamics
% RK4-Integrator
F_moon = @(x, h) rk4_moon(x, h);
     F = @(x, u, p_M, h) rk4_shuttle(x, u, p_M, h);

%% Simulation of moon (required for terminal conditions)
for i=1:N
    % Compute changes in velocity and acceleration.
    moon_step = F_moon([p_M(:, i); v_M(:, i)], delta_t);
    p_M(:,i+1) = moon_step(1:2);
    v_M(:,i+1) = moon_step(3:4);
end


%% formulate and solve the NLP
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
    xnext = F(xk, uk, p_M(:, k+1), delta_t);

    % Euler-Cromer Step
    % pos_next = xk(1:2) + xk(3:4)*delta_t;
    % vel_next = xk(3:4) + (-G*M_E*xk(1:2)/norm(xk(1:2))^3 - G*M_M*(xk(1:2)-p_M(:,k+1))/norm(xk(1:2)-p_M(:,k+1))^3 + uk)*delta_t;
    % xnext = [pos_next; vel_next];

    % contribution of stage cost to total cost
    J = J + uk(1)^2 + uk(2)^2;

    % New NLP variable for state at end of interval
    xk = MX.sym(['x_', num2str(k+1)], n_x);

    % collect in w
    w = {w{:}, xk};

    % Add dynamic constraint function,
    % constraining xk to integration result (norm(xk(1:2)-p_M(:,k+1))
    g = {g{:},  xk - xnext, uk(1)^2 + uk(2)^2, (xk(1:2)'* xk(1:2))};

    % Constrain shuttle to not fly through earth (moon).
    lbg = [lbg; zeros(n_x, 1); -(acc_limit)^2; ((6378000 + (35786000))/AU)^2]; %35786000
    ubg = [ubg; zeros(n_x, 1); acc_limit^2;  inf];
end

% Constrain total energy consumption
g = {g{:}, J};
lbg = [lbg; 0];
ubg = [ubg; consumption_limit];

% Terminal cost

%J = J + alpha*(xk(1:2) - (p_M(: , end) + [0; -(1737000+200000)/AU]))' * (xk(1:2) - (p_M(: , end) + [0; -(1737000+200000)/AU])) ;
%J = J + alpha*((xk(3) - v_M(1 , end))+(xk(4) - v_M(2 , end)));

% Terminal state constraints:
% Close to moon
% use epsilon distance to get close to moon
epsilon = MX.sym('epsilon', n_p);
g_terminal_pos = (xk(1:2) - (p_M(:, end) + [0; -(1737000+200000)/AU]))' * (xk(1:2) - (p_M(:, end) + [0; -(1737000+200000)/AU]))  - epsilon(1)^2;
g = {g{:}, g_terminal_pos};
lbg = [lbg; -inf];
ubg = [ubg; 0];

%constraint to not get to close to moon
g_terminal_pos = (xk(1:2) - p_M(:, end))' * (xk(1:2) - p_M(:, end)) + ((1737000)/AU)^2;
g = {g{:}, g_terminal_pos};
lbg = [lbg; 0];
ubg = [ubg; inf];

% Close to speed of moon (or perhaps to moon orbit speed)
g_terminal_vel = (xk(3:4) - v_M(:, end))' * (xk(3:4) - v_M(:, end)) - epsilon(2)^2;
g = {g{:}, g_terminal_vel};
lbg = [lbg; -inf];
ubg = [ubg; 0];

%% Initial guess.
% You have two options: Create an initial guess using the given script.
% This is the default mode and you don't have to do anything about it.
w0 = create_initial_guess(duration, delta_t, p_S, v_S, p_M, v_M);

% If you would like to create your own initial guess you can load it here.
% Remember that it has to fit the dimensions.
%load init_guess2.mat

%% Create an NLP solver
%vertcat({w}) will put all elements of w into a column vector
prob = struct('f', J, 'p', epsilon, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
options = struct;
options.ipopt.max_iter = 1000;
solver = nlpsol('solver', 'ipopt', prob, options);

%% Solve the NLP

sol = w0;
sol_all = w0;
epsilon_val = epsilon_vel_pos_start;
tic
for i = 1:number_of_casadi_Iter
    sol = solver('x0', sol, 'lbg', lbg, 'ubg', ubg, 'p', epsilon_val);
    sol = full(sol.x);
    sol_all = [sol_all, sol];
    epsilon_val = epsilon_val./2;
end
time_opti = toc;
w_opt = sol_all(:,end);

for i=1:size(sol_all,2)
    sol_all_mat(:,:,i) = reshape(sol_all(:,i), [6,N]);
end

%% Create plots.
close all;
% Plot trajectory of moon
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
subplot(3,2,[1,3]);
title('Position of shuttle and moon','Interpreter','latex');
xlabel('$p_x [AU]$','Interpreter','latex')
ylabel('$p_y [AU]$','Interpreter','latex')
axis([-1.2 1.2 -1.2 1.2]);
hold on;
viscircles([0,0],6378000/AU, 'Color','b');
%alpha_factor = 0.8/10;
for i=1:N
    plot(p_M(1,i), p_M(2,i), 'ob');
    plot(sol_all_mat(3,i,1), sol_all_mat(4,i,1), '+y');
    for j=2:(size(sol_all_mat,3)-1)
        plt = plot(sol_all_mat(3,i,j), sol_all_mat(4,i,j), '+g');  % the different optimizations
        %plt.Color(4) = 1-(0.2+(alpha_factor*j));
    end
    plot(sol_all_mat(3,i,end), sol_all_mat(4,i,end), 'or')  % final solution
%     plot(w_initial_opt_mat(3,i), w_initial_opt_mat(4,i), '+g');
%     plot(w_opt_mat(3,i), w_opt_mat(4,i), 'or');
    pause(0.1)
end

subplot(3,2,2);
hold on;
title('Control of shuttle x-direction','Interpreter','latex');
xlabel('N');
ylabel('$acceleration [\frac{m}{s^2}]$','Interpreter','latex');
stairs([0:N-1], sol_all_mat(1,:,end)/a_correct)

subplot(3,2,4);

hold on;
title('Control of shuttle y-direction','Interpreter','latex');
xlabel('N');
ylabel('$acceleration [\frac{m}{s^2}]$','Interpreter','latex');
stairs([0:N-1], sol_all_mat(2,:,end)/a_correct)

norm(w_opt(end-2:end))  %final velocity

subplot(3,2,5);
hold on;
title('Velocity of shuttle x-direction','Interpreter','latex');
xlabel('N');
ylabel('$velocity [\frac{m}{s}]$','Interpreter','latex');
plot([0:N-1], sol_all_mat(5,:,end)/v_correct)


subplot(3,2,6);
hold on;
title('Velocity of shuttle y-direction','Interpreter','latex');
xlabel('N')
ylabel('$velocity [\frac{m}{s}]$','Interpreter','latex')
plot([0:N-1], sol_all_mat(6,:,end)/v_correct)
time_opti