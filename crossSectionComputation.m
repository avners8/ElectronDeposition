function crossSection = crossSectionComputation(control)
    % Source: "Deterministic model for dose calculation in photon
    % radiotherapy" by Hartmut Hensel (2006)
    % Output:
    % - crossSection - struct containing the cross sections for all scattering types
    %                  The output dim are [Energy X Energy_tag] for
    %                  differential cross sections and [Energy] for total
    %                  cross sections

    crossSection = struct();
    crossSection = ComptonPhotonCrossSectionComputation(crossSection, control);
    crossSection = ComptonElectronCrossSectionComputation(crossSection, control);
    crossSection = MottCrossSectionComputation(crossSection, control);
    crossSection = MollerCrossSectionComputation(crossSection, control);
        
%     crossSection.Electrons_Tot =  crossSection.T_Moller  + crossSection.T_Mott;
    crossSection.Electrons_Tot =  control.factor*(crossSection.Moller_Tot+ crossSection.Mott_Tot);
    crossSection.Electrons_0   =  control.factor*(crossSection.Moller_0  + crossSection.Mott_0  ); 
    crossSection.Electrons_1   =  control.factor*(crossSection.Moller_1  + crossSection.Mott_1  );
    
end

function crossSection = ComptonPhotonCrossSectionComputation(crossSection, control)
    
    r_e = control.general.r_e;
    eps = repmat(control.eps, [1 control.energy_tag_length control.mu_length]);
    eps_tag = permute(repmat(control.eps, [1 control.energy_length control.mu_length]), [2 1 3]);
    mu = permute(repmat(control.mu, [1 control.energy_tag_length control.energy_length]), [3 2 1]);
    
    
    crossSection.Compton_photon = r_e^2/2 * (1 + eps_tag .* (1 - mu)).^(-2) .* (1 + mu.^2 + (eps_tag .* (1 - mu)).^2 ./ (1 + eps_tag .* (1 - mu)));
    crossSection.Compton_photon = control.factor * crossSection.Compton_photon .* (eps == (eps_tag ./ (1 - eps_tag .* (1 - mu))));
    
    crossSection.Compton_photon_0 = 2*pi* sum(crossSection.Compton_photon, 3)       * (2/control.mu_length);
    crossSection.Compton_photon_1 = 2*pi* sum(mu .* crossSection.Compton_photon, 3) * (2/control.mu_length);
    
    eps = control.eps;
    crossSection.Compton_photon_Tot = control.factor * 2*pi*r_e^2 * ((1+eps)./(eps.^2).*(2.*(1+eps)./(1+2.*eps) - log(1+2.*eps)./eps) + log(1+2.*eps)./(2*eps) - (1+3.*eps)./((1+2*eps).^2));
end

function crossSection = ComptonElectronCrossSectionComputation(crossSection, control)
    
    r_e = control.general.r_e;
    mu = permute(repmat(control.mu, [1 control.energy_tag_length control.energy_length]), [3 2 1]);
    eps_tag = permute(repmat(control.eps, [1 control.energy_length control.mu_length]), [2 1 3]);
    a = (1 + eps_tag).^2 ./((mu.^2) - 1) + 1;
    eps = 2*eps_tag.^2 ./ (2*eps_tag + a);
    
    crossSection.Compton_electron = control.factor*4 * r_e^2 * (1 + eps_tag).^2 ./ (mu.^3) ./ ((a + 2).^2) .* ...
        (1 - 2./a + 2./(a.^2) + 2 * eps.^2 ./ (a.*(a + 2.*eps_tag)));
    
    crossSection.Compton_electron_0 = 2*pi* sum(crossSection.Compton_electron, 3)       * (2/control.mu_length);
    crossSection.Compton_electron_1 = 2*pi* sum(mu .* crossSection.Compton_electron, 3) * (2/control.mu_length);
end

function crossSection = MottCrossSectionComputation(crossSection, control)
    
    Z = control.Z; % Atomic number of the medium
    r_e = control.general.r_e;
    eps = repmat(control.eps, [1 control.energy_tag_length control.mu_length]);
    mu = permute(repmat(control.mu, [1 control.energy_tag_length control.energy_length]), [3 2 1]);
    eps_tag = permute(repmat(control.eps, [1 control.energy_length control.mu_length]), [2 1 3]);
    
    eta = (pi * control.general.alpha_fsc * Z^(1/3)) .^2 ./ (eps .*(eps + 2));
    beta_2 = (eps .*(eps + 2)) ./ ((eps + 1).^2);
    crossSection.Mott = (Z .* r_e .* (1 + eps)).^2 .* (1 - beta_2 .* (1 - mu)./2) ...
                    ./ (4 .* (eps .*(eps + 2)).^2 .* (1 + 2.* eta - mu).^2);
    crossSection.Mott = control.factor*crossSection.Mott .* (eps == eps_tag);
    
    crossSection.Mott_0 = 2*pi* sum(crossSection.Mott, 3)       * (2/control.mu_length);
    crossSection.Mott_1 = 2*pi* sum(mu .* crossSection.Mott, 3) * (2/control.mu_length);

    eps = control.eps;
    eta = (pi * control.general.alpha_fsc * Z^(1/3)) .^2 ./ (eps .*(eps + 2));
    crossSection.Mott_Tot = control.factor*pi * (Z * r_e).^2 ./ (eps .*(eps + 2)) .* (...
          ((eps + 1).^2)./((pi * control.general.alpha_fsc * Z^(1/3)) .^2 .* (1 + eta)) ...
        + 1 ./ (1 + eta) ...
        + log(eta) ...
        - log(1 + eta));

    eps = repmat(control.eps, [1 control.x_length]);
    rho = permute(repmat(control.rho, [1 control.energy_length]), [2 1]);
    eta = (pi * control.general.alpha_fsc * Z^(1/3)) .^2 ./ (eps .*(eps + 2));
    beta_2 = (eps .*(eps + 2)) ./ ((eps + 1).^2);
    crossSection.T_Mott = -pi*((Z*r_e*(eps+1)).^2) .* rho ./ ((eps.*(eps+2)).^2 .* (1 + eta)) .* ...
        (1 + (beta_2).*(1 + 2.*eta) + (1 + eta).*(1 + 2*(beta_2).*eta).*(log(2*eta./(2+2.*eta))));
    crossSection.T_Mott = squeeze(crossSection.T_Mott(:,1));
end

function crossSection = MollerCrossSectionComputation(crossSection, control)
    
    eps_b = control.general.eps_b; % Binding energy
    r_e = control.general.r_e;
    eps = repmat(control.eps, [1 control.energy_tag_length control.mu_length]);
    eps_tag = permute(repmat(control.eps, [1 control.energy_length control.mu_length]), [2 1 3]) + 1e-8;
    mu = sqrt((eps./eps_tag) .* ((eps_tag + 2)./(eps + 2)));
    W = eps_tag - eps;
    
    crossSection.Moller = (r_e * (eps_tag + 1) ./ (eps_tag .* (eps_tag + 2))).^2 .* ...
                        (1./(W.^2) + 1./((eps_tag - W).^2) + 1./((eps_tag + 1).^2) ...
                        -(2*eps_tag + 1)./((eps_tag + 1).^2)./(W.*(eps_tag - W)) );

    crossSection.Moller = control.factor*crossSection.Moller .* (eps > eps_b) .* (eps < (eps_tag - eps_b));
    
    crossSection.Moller_0 = 2*pi* sum(crossSection.Moller, 3)       * (2/control.mu_length);
    crossSection.Moller_1 = 2*pi* sum(mu .* crossSection.Moller, 3) * (2/control.mu_length);

    delta_eps = repmat(control.delta_eps, [1 control.energy_tag_length]);
    crossSection.T_Moller = sum(sum((1 - mu) .* crossSection.Moller .* control.delta_mu, 3) .* delta_eps, 2);
    
    eps = control.eps;
    crossSection.Moller_Tot = control.factor*2 * pi * r_e^2 .* (eps + 1).^2 ./ (eps .* (eps + 2)) .* (...
          1./eps_b ...
        - 3./(eps - eps_b) ...
        + 2./(eps + eps_b) ...
        + (eps - 3*eps_b)./(2*(eps + 1).^2) ...
        + (2*eps + 1) ./ (eps .* (eps + 1)) .* (log( (eps + eps_b)./(eps - eps_b)) - log( (eps - eps_b)./eps_b)) ); 
    
end
