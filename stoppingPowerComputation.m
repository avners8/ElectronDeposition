function [S, control] = stoppingPowerComputation(control)

    [S, control] = stoppingPowerComputation_aux(control);
    
    % Computing delta_epsilon based on the stopping power
    % Source: "A numerical approach for a system of transport
    % equations in the field of radiotherapy" by Teddy Pichard
    
    control.delta_eps = 1 * S(:,1) * control.delta_x;
    control.eps = flipud(control.eps_min + cumsum(control.delta_eps));
    [~, eps_0_ind] = min(abs(control.eps_0 - control.eps)); eps_0_ind = eps_0_ind(1);
    control.eps(eps_0_ind) = control.eps_0;
    ind = control.eps < control.eps_max;
    control.energy_length  = sum(ind);
    control.energy_tag_length  = control.energy_length;
    control.eps = control.eps(ind); control.delta_eps = control.delta_eps(ind);
    [S, control] = stoppingPowerComputation_aux(control);
    
end

function [S, control] = stoppingPowerComputation_aux(control)
    % Computing the Moller stopping power
    % Source: "Deterministic model for dose calculation in photon
    % radiotherapy" by Hartmut Hensel (2006)
    
    eps_b = control.general.eps_b; % Binding energy
    r_e = control.general.r_e;
    rho = permute(repmat(control.rho, [1 control.energy_length]), [2 1]);
    eps_e = repmat(control.eps, [1 control.x_length]);
    
    S = 20e4*2*pi*(r_e^2)* rho .* (eps_e + 1).^2 ./ (eps_e .* (eps_e + 1)) ...
        .* ( eps_e ./ (eps_e - eps_b) + 2.*log((eps_e - eps_b)./(2.*eps_b.*(eps_e - eps_b))) ...
        + 1./(2*(eps_e + 1).^2) .* (((eps_e - eps_b).^2)./4 - eps_b.^2) - (2*eps_e + 1)./((eps_e + 1).^2) .* (log(2)) );
end