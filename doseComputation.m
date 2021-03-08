function [x, dose] = doseComputation
    %% Dose computation
    % General definitions

    general.r_e = 2.8179e-7;        % Classical electron radius [cm]
    general.alpha_fsc = 1/137;      % fine-structure-constant
    general.eps_b = 16e-6;          % Binding energy 16 eV
    control.general = general;

    % Grid sizes
    control.energy_length  = 50;
    control.x_length       = 200;
    control.mu_length      = 100;
    control.max_algorithm_iterations = 400;

    % Calculation parameters
    control.energy_tag_length  = control.energy_length;
    control.types          = 4;         % {photon,electron} X {0,1}
    control.BW             = 5;
    control.ce             = 200;
    control.c0             = 1e3;
    control.omega_p        = 1;
    control.x0             = 0;
    control.eps_0          = 200;        % MeV
    control.eps_min        = 1e-3; %0.7 * control.eps_0;  %
    control.eps_max        = 1.3 * control.eps_0;  % May be proportionate to control.eps_0
    control.mu             = linspace(-1+1e-15, 1-1e-15, control.mu_length)';
    control.mu(control.mu == 0) = 1e-8;
    control.Z              = 9.21;      % Atomic numbre of the medium. Can be dependent on the location.
    control.scale_x        = 5.2e-7;
    control.x_max          = 35 / control.scale_x;
    control.x              = linspace(0, control.x_max, control.x_length)';
    control.rho            = ones(size(control.x));
    control.delta_x        = control.x(2) - control.x(1);
    control.delta_mu       = control.mu(2) - control.mu(1);
    control.delta          = 1e-15;
    control.tolerance      = 1e-10;
    control.K              = 10e3;
    control.factor         = 1e-3;
    x = control.x * control.scale_x;
    control.eps = flipud(linspace(control.eps_min, control.eps_max, control.energy_length)');
    control.draw_dose = 0;

    tic
    [S, control] = stoppingPowerComputation(control);
    crossSection = crossSectionComputation(control);
    crossSection.S = S;
    psi = fluenceComputation(crossSection, control);
    [dose, flux] = doseComputation_aux(psi, crossSection, control);
    toc

end
