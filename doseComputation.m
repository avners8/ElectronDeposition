set(groot, 'defaultFigurePosition', [100 100 900 600]); % figure size
set(groot, 'defaultTextInterpreter', 'latex'); % latex
set(groot, 'defaultLegendInterpreter', 'latex'); % latex
set(groot, 'defaultAxesTickLabelInterpreter','latex'); % latex
set(0, 'DefaultLineLineWidth', 4);
lighter = - 0.071;
Green = [0.4660 0.6740 0.1880] - lighter; Red   = [0.6350 0.0780 0.1840] - lighter; Gold  = [0.9290 0.6940 0.1250] - lighter; Blue  = [0.1290 0.3400 0.8850] - lighter;
Green_Sim = [0.4660 0.6740 0.1880] + lighter; Red_Sim   = [0.6350 0.0780 0.1840] + lighter; Gold_Sim  = [0.9290 0.6940 0.1250] + lighter; Blue_Sim  = [0.1290 0.3400 0.8850] + lighter;

%% Dose computation
% General definitions

general.r_e = 2.8179e-7;        % Classical electron radius [cm]
general.alpha_fsc = 1/137;      % fine-structure-constant
general.eps_b = 16e-6;          % Binding energy 16 eV
control.general = general;

control.energy_length  = 1000;
control.x_length       = 80;
control.mu_length      = 200;
control.max_algorithm_iterations = 100;

control.energy_tag_length  = control.energy_length;
control.types          = 4;         % {photon,electron} X {0,1}
control.ce             = 200;
control.c0             = 1e3;
control.omega_p        = 1;
control.x0             = 0;
control.eps_0          = 20;        % MeV
control.eps_min        = 1e-0;
control.eps_max        = 1.2 * 40;  % May be proportionate to control.eps_0
control.mu             = linspace(-1+1e-15, 1-1e-15, control.mu_length)';
control.mu(control.mu == 0) = 1e-8;
control.Z              = 9.21;      % Atomic numbre of the medium. Can be dependent on the location.
control.x              = linspace(0, 10e7/6, control.x_length)';
control.rho            = ones(size(control.x));
control.delta_x        = control.x(2) - control.x(1);
control.delta_mu       = control.mu(2) - control.mu(1);
control.delta          = 1e-15;
control.tolerance      = 1e-3;
control.K              = 10e0;
control.factor         = 1e-1;

control.eps = flipud(linspace(control.eps_min, control.eps_max, control.energy_length)');
control.draw_dose = 0;

tic
[S, control] = stoppingPowerComputation(control);
crossSection = crossSectionComputation(control);
crossSection.S = S;
psi = fluenceComputation(crossSection, control);
[dose, flux] = doseComputation_aux(psi, crossSection, control);
toc

if ~control.draw_dose
    load(['Sim_', num2str(control.eps_0), 'MeV.mat']);
    figure();
    plot(control.x*1e-7*1.93*4/1.1, real(dose) ./ max(abs(real(dose))), 'DisplayName', 'Dose Theo'); hold on; 
    plot(control.x*1e-7*1.93*4/1.1, real(flux) ./ max(abs(real(flux))), 'DisplayName', 'Flux Theo'); hold on; 
    plot(Sim.x_sim*1e6*1e-7, Sim.dose_sim / max(Sim.dose_sim), 'DisplayName', 'Dose Sim'); hold on;
    graphParams(['Dose  -  $\epsilon_0 = $ ', num2str(control.eps_0), 'MeV'], 'x', '$D(x)/D_{max}$');
end

function graphParams(ptitle, pxlabel, pylabel)
    grid on; title(ptitle); xlabel(pxlabel); ylabel(pylabel); set(gca, 'FontSize', 14); set(gcf,'color','w'); set(gca,'linewidth',2); legend('show');
end