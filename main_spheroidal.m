%%% Simulating a spinning spheroidal swimmer.

%% Setup.
% G = gamma, the shear rate of the flow.
G = 1;
% B = Bretherton constant.
B = 0.8;

% W_par is the intrinsic spin of the swimmer about its axis of symmetry.
W_par = 10;
% W_perp is the other direction of spin (perpendicular to the axis of symmetry).
W_perp = 1;

% Swimming velocity along axis of symmetry, e_hat_1.
V1 = 1;
% Swimming velocity along axis e_hat_2.
V2 = 1;
% Swimming velocity along axis e_hat_3.
V3 = 0;

% Initial condition for swimmer position.
X0 = [0;0;0];

% Initial conditions for swimmer orientation as theta, phi, and psi.
init_theta  =   pi*rand;
init_phi    = 2*pi*rand;
init_psi    = 2*pi*rand;

%-------------

% ODE options.
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

%% Derived quantities used in the asymptotic analysis.
w = W_perp / W_par; % spinning ratio.
lambda = sqrt(1 + w^2); 
% Effective Bretherton constant.
B_eff = B * (2 - w^2) / (2 * (1 + w^2));
% Effective speed.
V_hat = (V1 + w*V2) / lambda;

% The ICs vector for the full simulations.
init_full = [init_theta; init_phi; init_psi; X0];

% Generate ICs for the long-time variables, alpha_bar, mu_bar, and phi_bar.
[init_alpha_bar, init_mu_bar, init_phi_bar] = ...
    initial_conditions(init_theta, init_phi, init_psi, W_perp, W_par);

% The ICs vector for the reduced simulations.
init_reduced = [init_alpha_bar; init_phi_bar; init_mu_bar; X0];

% Time interval.
t_min = 0;
t_max = 20;
ts = linspace(t_min,t_max,1e5);

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

Ts = sqrt(W_par.^2 + W_perp.^2)*ts';

%% Plotting.
% We'll plot the asymptotic solution alongside a ribbon plot of the full solution.
figure

% Asymptotic solution.
step = 100;
plot3(x_bar(1:step:end),y_bar(1:step:end),z_bar(1:step:end),'r','LineWidth',2);
grid on
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex'); 
hold on

% Ribbon plot.
plot3(x_full,y_full,z_full,'k','LineWidth',2);
vertices = {[x_full, y_full, z_full]};
psi_mod = mod(psi_full, 2*pi) - pi;
twistangle = {[psi_mod(1);diff(psi_mod)]};
str = streamribbon(vertices, twistangle ,0.2);

shading interp
lighting gouraud
material dull
camlight
axis equal
col = parula(100);
col = [col(end:-1:1,:);col];
colormap(col)
c = colorbar;
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.Label.String = '$\psi$ mod $2\pi$';

% Components of position.
figure
hold on
plot(ts,x_full,'LineWidth',2)
plot(ts,y_full,'LineWidth',2)
plot(ts,z_full,'LineWidth',2)
plot(ts,x_bar,'LineWidth',2,'LineStyle','--')
plot(ts,y_bar,'LineWidth',2,'LineStyle','--')
plot(ts,z_bar,'LineWidth',2,'LineStyle','--')
legend('$x$','$y$','$z$','$x_0$','$y_0$','$z_0$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$(x,y,z)$, $(x_0,y_0,z_0)$','Interpreter','latex')

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
    d_state(2) = w1*sin(psi) ./ sin(theta) + G*f2;
    d_state(3) = w2 - w1*sin(psi).*cot(theta) + G*f3;

    % Translational dynamics.
    d_state(4) =  V1*sin(phi)*sin(theta) + ...
                  V2*( cos(phi)*cos(psi) - cos(theta)*sin(phi)*sin(psi)) + ...
                  V3*(-cos(phi)*sin(psi) - cos(theta)*sin(phi)*cos(psi));

    d_state(5) = -V1*cos(phi)*sin(theta) + ...
                  V2*( sin(phi)*cos(psi) + cos(theta)*cos(phi)*sin(psi)) + ...
                  V3*(-sin(phi)*sin(psi) + cos(theta)*cos(phi)*cos(psi));

    d_state(6) =  V1*cos(theta) + G*y + ...
                  V2*sin(theta)*sin(psi) + ...
                  V3*sin(theta)*cos(psi);
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
    d_state(4) =  V_hat*sin(phi_bar)*sin(alpha_bar);
    d_state(5) = -V_hat*cos(phi_bar)*sin(alpha_bar);
    d_state(6) =  V_hat*cos(alpha_bar) + G*y_bar;
end

function [alp_init, mu_init, phibar_init] = initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par )
%% Compute the appropriate initial conditions for the reduced system from the initial
% conditions of the full system and the components of the spin.

omega = W_perp/W_par;
lambda = sqrt(1 + omega^2);

% Calculate IC for alpha.
lambda_cos_alp_init = omega*sin(Init_theta)*sin(Init_psi) + cos(Init_theta);
alp_init = acos(lambda_cos_alp_init./lambda);

% Calculate IC for mu.
tan_mu_init_over_lambda = cos(Init_psi)*sin(Init_theta)./(omega*cos(Init_theta) - sin(Init_theta)*sin(Init_psi));
mu_init = atan(lambda*tan_mu_init_over_lambda);

% Use additional constraints to determine if we're on the right pi-branch for mu.
if abs(sin(alp_init)*sin(mu_init) + cos(Init_psi)*sin(Init_theta)) < 1e-10 ...
                && abs(lambda*sin(alp_init)*cos(mu_init) + omega*cos(Init_theta) - sin(Init_theta)*sin(Init_psi)) < 1e-10
else
    mu_init = mu_init + pi;
end

% Calculate phibar IC.
phibar_init = Init_phi - atan2(omega*sin(mu_init),omega*cos(alp_init).*cos(mu_init) + sin(alp_init));

% Perform a consistency check to make sure we're on the right pi-branch for
% phibar.
abs_lambda_sin_theta = sqrt((omega*cos(alp_init).*cos(mu_init) + sin(alp_init)).^2 ...
    + (omega*sin(mu_init)).^2);
if abs(abs_lambda_sin_theta*cos(Init_phi) - ...
        (cos(phibar_init).*(omega*cos(alp_init).*cos(mu_init) + sin(alp_init)) ...
        - omega*sin(phibar_init).*sin(mu_init))) < 1e-10
else
    phibar_init = phibar_init + pi;
end

end