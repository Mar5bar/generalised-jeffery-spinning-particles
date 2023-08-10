function plot_sphere(theta_full,phi_full,alpha_bar,phi_bar)
%% Plot angular quantities on the unit sphere.
    hold on
    % Plot the unit sphere.
    [sx,sy,sz]=sphere(30);
    s = surf(sx,sy,sz,'FaceLighting','gouraud','FaceAlpha',0.35,'FaceColor',[1 1 1],'EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',0.4);
    s.LineWidth = 1;
    material dull

    % Plot the axes and the xy unit circle
    plot3([0 1.5],[0 0],[0 0],'k','LineWidth',1)
    plot3([0 0],[0 1.5],[0 0],'k','LineWidth',1)
    plot3([0 0],[0 0],[0 1.5],'k','LineWidth',1)
    tt = 0:0.1:2*pi+0.1;
    plot3(cos(tt),sin(tt),zeros(1,length(tt)),'k','LineWidth',1)

    % Full solution.
    plot3(sin(theta_full).*cos(phi_full),sin(theta_full).*sin(phi_full),cos(theta_full),'b')

    % Average solution.
    plot3(sin(alpha_bar).*cos(phi_bar),sin(alpha_bar).*sin(phi_bar),cos(alpha_bar),'r','LineWidth',4)

    % Mark the initial conditions. 
    plot3(sin(theta_full(1)).*cos(phi_full(1)),sin(theta_full(1)).*sin(phi_full(1)),cos(theta_full(1)),'k.','MarkerSize',30)
    plot3(sin(alpha_bar(1)).*cos(phi_bar(1)),sin(alpha_bar(1)).*sin(phi_bar(1)),cos(alpha_bar(1)),'k.','MarkerSize',30)

    view(120,25)
    camlight
    axis equal off
    grid on
    box on
    
end