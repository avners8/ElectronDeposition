function [dose, flux] = doseComputation_aux(psi, crossSection, control)
    % Compute the dose
    % Q   = (Q_gg_0 + Q_eg_0, Q_ee_0 + Q_ge_0, Q_gg_1 + Q_eg_1, Q_ee_1 + Q_ge_1)
    Q = zeros(size(psi));
    delta_eps = repmat(control.delta_eps, [1 control.x_length]);
    crossSection_Compton_photon_Tot = repmat(crossSection.Compton_photon_Tot, [1 control.x_length]);
    crossSection_electron_Tot       = repmat(crossSection.Electrons_Tot, [1 control.x_length]);
    
    Q(:,:,1) = crossSection.Compton_photon_0 * squeeze(psi(:,:,1)) .* delta_eps - crossSection_Compton_photon_Tot .* psi(:,:,1);
    Q(:,:,2) = crossSection.Compton_photon_1 * squeeze(psi(:,:,2)) .* delta_eps - crossSection_Compton_photon_Tot .* psi(:,:,2);
    
    Q(:,:,3) = crossSection.Compton_electron_0 * squeeze(psi(:,:,1)) .* delta_eps;
    Q(:,:,4) = crossSection.Compton_electron_1 * squeeze(psi(:,:,2)) .* delta_eps;
    
    S = repmat(crossSection.S, [1 1 control.types]);
    S_psi_shifted = S .* psi;
    S_psi_shifted((2:(control.energy_length)),:,:) = S_psi_shifted((1:(control.energy_length-1)),:,:);
    
    Q(:,:,3) = Q(:,:,3) + (S_psi_shifted(:,:,3) - S(:,:,3) .* psi(:,:,3)) ./ delta_eps;
    Q(:,:,4) = Q(:,:,4) + (S_psi_shifted(:,:,4) - S(:,:,4) .* psi(:,:,4)) ./ delta_eps;
    
    Q(:,:,3) = Q(:,:,3) + crossSection.Electrons_0 * squeeze(psi(:,:,3)) .* delta_eps - crossSection_electron_Tot .* squeeze(psi(:,:,3));
    Q(:,:,4) = Q(:,:,4) + crossSection.Electrons_1 * squeeze(psi(:,:,4)) .* delta_eps - crossSection_electron_Tot .* squeeze(psi(:,:,4));
    
    eps = repmat(control.eps, [1 control.x_length control.types]);
    dose = sum(sum(-eps .* Q, 3) .* delta_eps, 1);
    
    flux = sum(repmat(crossSection.S(:,1), [1 control.x_length]) .* squeeze(psi(:,:,1) + psi(:,:,3)) .* delta_eps, 1);

end