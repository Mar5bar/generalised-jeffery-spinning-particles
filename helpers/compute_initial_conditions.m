function [alp_init, mu_init, phibar_init] = compute_initial_conditions(Init_theta, Init_phi, Init_psi, W_perp, W_par )
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
    abs_lambda_sin_theta = sqrt((omega*cos(alp_init)*cos(mu_init) + sin(alp_init))^2 ...
        + (omega*sin(mu_init))^2);
    if abs(abs_lambda_sin_theta*cos(Init_phi) - ...
            (cos(phibar_init)*(omega*cos(alp_init)*cos(mu_init) + sin(alp_init)) ...
            - omega*sin(phibar_init)*sin(mu_init))) < 1e-10
    else
        phibar_init = phibar_init + pi;
    end

end