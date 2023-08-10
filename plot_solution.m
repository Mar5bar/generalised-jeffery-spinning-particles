addpath(genpath('./helpers'))

% We'll plot the asymptotic solution alongside a ribbon plot of the full solution.
figure("Position",[0,0,1400,370])
tiledlayout(1,3);

% Components of position.
nexttile()
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

% Asymptotic solution.
nexttile();
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

% Plotting angular variables on a sphere.
nexttile()
plot_sphere(theta_full,phi_full,alpha_bar,phi_bar);