%%% Simulating a spinning spheroidal swimmer.

%% Setup.
% G = gamma, the shear rate of the flow.
G = 1;

% B = Bretherton parameter.
B = 0.8;

% W_par is the intrinsic spin of the swimmer about its axis of symmetry.
W_par = 20;
% W_perp is the other direction of spin (perpendicular to the axis of symmetry).
W_perp = 20;

% Swimming velocity along axis of symmetry, e_hat_1.
V1 = 1;
% Swimming velocity along axis e_hat_2.
V2 = 0;
% Swimming velocity along axis e_hat_3.
V3 = 0;

% Initial condition for swimmer position.
X0 = [0;0;0];

% Initial conditions for swimmer orientation as theta, phi, and psi.
init_theta  =    pi/3;
init_phi    = -5*pi/6;
init_psi    =  2*pi/3;

% Time interval.
t_min = 0;
t_max = 20;
ts = linspace(t_min,t_max,1e5);

%-------------
addpath(genpath('./helpers'))

% ODE options.
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

%% Derived quantities used in the asymptotic analysis.
w = W_perp / W_par; % spinning ratio.
lambda = sqrt(1 + w^2); 
% Effective Bretherton parameter.
B_eff = B*(2 - w^2)/(2*(1 + w^2));
% Effective speed.
V_hat = (V1 + w*V2) / lambda;

% The ICs vector for the full simulations.
init_full = [init_theta; init_phi; init_psi; X0];

% Generate ICs for the long-time variables, alpha_bar, mu_bar, and phi_bar.
[init_alpha_bar, init_mu_bar, init_phi_bar] = ...
    compute_initial_conditions(init_theta, init_phi, init_psi, W_perp, W_par);

% The ICs vector for the reduced simulations.
init_reduced = [init_alpha_bar; init_phi_bar; init_mu_bar; X0];

params          = struct();
params.G        = G;
params.B        = B;
params.W_perp   = W_perp;
params.W_par    = W_par;
params.w        = w;
params.B_eff    = B_eff;
params.lambda   = lambda;
params.V        = [V1; V2; V3];
params.V_hat    = V_hat;

%% Solve full ODE.
tic
[~, sol_full] = ode15s(@(t,state) ode_full(t,state,params),ts,init_full,options);
toc

theta_full  = sol_full(:,1);
phi_full    = sol_full(:,2);
psi_full    = sol_full(:,3);
x_full      = sol_full(:,4);
y_full      = sol_full(:,5);
z_full      = sol_full(:,6);

%% Solve the reduced ODE.
tic
[~, sol_reduced] = ode15s(@(t,state) ode_reduced(t,state,params),ts,init_reduced,options);
toc

alpha_bar   = sol_reduced(:,1);
phi_bar     = sol_reduced(:,2);
mu_bar      = sol_reduced(:,3);
x_bar       = sol_reduced(:,4);
y_bar       = sol_reduced(:,5);
z_bar       = sol_reduced(:,6);

Ts = sqrt(W_par^2 + W_perp^2)*ts';

%% Plotting.
plot_solution;

%% Auxiliary functions 
function d_state = ode_full(t,state,params)
%% The RHS of the full ODEs for both the angular and translational
% dynamics. Parameters are passed via the structure params.
    d_state = zeros(6,1);

    % Unpack the states.
    theta   = state(1);
    phi     = state(2);
    psi     = state(3);
    x       = state(4);
    y       = state(5);
    z       = state(6);

    % Unpack the params.
    B   = params.B;
    G   = params.G;
    w1  = params.W_perp;
    w2  = params.W_par;
    V1  = params.V(1);
    V2  = params.V(2);
    V3  = params.V(3);

    % Angular dynamics.
    f1 = -B*(sin(2*theta)*sin(2*phi))/4;
    f2 = (1 - B*cos(2*phi))/2;
    f3 = (B/2)*cos(theta)*cos(2*phi);

    d_state(1) = w1*cos(psi) + G*f1;
    d_state(2) = w1*sin(psi)/sin(theta) + G*f2;
    d_state(3) = w2 - w1*sin(psi)*cot(theta) + G*f3;

    % Translational dynamics.
    d_state(4) = V1*cos(theta) + ...
                 V2*sin(theta)*sin(psi) + ...
                 V3*sin(theta)*cos(psi);

    d_state(5) = V1*sin(phi)*sin(theta) + ...
                 V2*( cos(phi)*cos(psi) - cos(theta)*sin(phi)*sin(psi)) + ...
                 V3*(-cos(phi)*sin(psi) - cos(theta)*sin(phi)*cos(psi));

    d_state(6) = -V1*cos(phi)*sin(theta) + G*y + ...
                  V2*( sin(phi)*cos(psi) + cos(theta)*cos(phi)*sin(psi)) + ...
                  V3*(-sin(phi)*sin(psi) + cos(theta)*cos(phi)*cos(psi));
end

function d_state = ode_reduced(t,state,params)
%% The RHS of the reduced ODEs for both the auxilliary angular functions and
% average translational dynamics. Parameters are passed via the structure
% params.
    d_state = zeros(6,1);

    % Unpack the states.
    alpha_bar   = state(1);
    phi_bar     = state(2);
    mu_bar      = state(3);
    x_bar       = state(4);
    y_bar       = state(5);
    z_bar       = state(6);

    % Unpack the params.
    B_eff   = params.B_eff;
    G       = params.G;
    V_hat   = params.V_hat;

    % Auxiliary angular dynamics.
    f1 = -B_eff*(sin(2*alpha_bar)*sin(2*phi_bar))/4;
    f2 = (1 - B_eff*cos(2*phi_bar))/2;
    f3 = (B_eff/2)*cos(alpha_bar)*cos(2*phi_bar);

    d_state(1) = G*f1;
    d_state(2) = G*f2;
    d_state(3) = G*f3;

    % Average translational dynamics.
    d_state(4) =  V_hat*cos(alpha_bar);
    d_state(5) =  V_hat*sin(phi_bar)*sin(alpha_bar);
    d_state(6) = -V_hat*cos(phi_bar)*sin(alpha_bar) + G*y_bar;
end