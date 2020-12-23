function [S, control] = stoppingPowerComputation(control)

    [S, control] = stoppingPowerComputation_aux(control);
    
    % Computing delta_epsilon based on the stopping power
    % Source: "A numerical approach for a system of transport
    % equations in the field of radiotherapy" by Teddy Pichard
    
    control.delta_eps = abs(control.eps(2) - control.eps(1));
   
end

function [S, control] = stoppingPowerComputation_aux(control)
    % Computing the Moller stopping power
    % Source: "Deterministic model for dose calculation in photon
    % radiotherapy" by Hartmut Hensel (2006)
    
    eps_b = control.general.eps_b; % Binding energy
    r_e = control.general.r_e;
    rho = permute(repmat(control.rho, [1 control.energy_length]), [2 1]);
    eps_e = repmat(control.eps, [1 control.x_length]);
    
    S = 1e5*2*pi*(r_e^2)* rho .* (eps_e + 1).^2 ./ (eps_e .* (eps_e + 1)) ...
        .* ( eps_e ./ (eps_e - eps_b) + 2.*log((eps_e - eps_b)./(2.*eps_b.*(eps_e - eps_b))) ...
        + 1./(2*(eps_e + 1).^2) .* (((eps_e - eps_b).^2)./4 - eps_b.^2) - (2*eps_e + 1)./((eps_e + 1).^2) .* (log(2)) );
end