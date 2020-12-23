function psi = fluenceComputation(crossSection, control)
    % Source: "A numerical approach for a system of transport
    % equations in the field of radiotherapy" by Teddy Pichard 
    % This function iterativley compute the first moments of the fluence until convergence with thw M1 approximation    
    % Output: [Energy X x(location) X type(4)]
    % psi = (psi_g_0, psi_g_1, psi_e_0, psi_e_1)
    % F   = (psi_g_1, psi_g_2, psi_e_1, psi_e_2)

    % Initialization
    psi = computeBoundaryCondition(zeros([control.energy_length control.x_length control.types]), control); 
    rho = permute(repmat(control.rho, [1 control.energy_length control.types]) , [2 1 3]);
    if control.draw_dose
        load(['Sim_', num2str(control.eps_0), 'MeV.mat']);
        figure();
    end
    diff_old = 1e80;
    % Iterative algorithm
    % We limit the algorithm to max_algorithm_iterations iterations 
    for i = 1 : control.max_algorithm_iterations
        F = compute_F(psi);
        R = computeR(psi, crossSection, control);
        
        L = (psi + F) ./ (2 * control.delta_x); 
        L(:, 2:control.x_length, :) = L(:, 1:(control.x_length - 1), :);

        U = (psi - F) ./ (2 * control.delta_x); 
        U(:, 1:(end-1), :) = U(:, 2:end, :);
        
        P = L + U + rho .* R;
        [psi_new, D] = compute_psi(psi, crossSection, P, control);
        
        % Checking if the fluence converged
        diff = max(abs(D - P), [], 'all');
        if ((diff > 1e30))% | (diff_old < diff - 0.001))
            psi_new = psi;
            diff = 0;
        end
        if (diff < control.tolerance)
            fprintf('Algorithm converged after %d iterations\n', i);
            break
        end
        diff_old = diff;
        % Plotting the dose obtained at each iteration
        if control.draw_dose
            [dose, flux] = doseComputation_aux(psi_new, crossSection, control);
            plot(control.x * control.scale_x * 1e-7, real(dose) ./ max(abs(real(dose))), 'DisplayName', 'Dose Theo'); hold on; 
            plot(control.x * control.scale_x * 1e-7, real(flux) ./ max(abs(real(flux))), 'DisplayName', 'Flux Theo'); hold on; 
            plot(Sim.x_sim*1e-1*control.x_max/10 * control.scale_x * 1e-7, Sim.dose_sim / max(Sim.dose_sim), 'DisplayName', 'Dose Sim'); hold on;
            graphParams(['Dose  -  $\epsilon_0 = $ ', num2str(control.eps_0), 'MeV'], 'x', '$D(x)/D_{max}$');
            drawnow;
            hold off;
        end
        
        psi = psi_new;

    end

    psi = psi_new;

end

function psi_bc = computeBoundaryCondition(psi, control)
    % Compute the boundary condition given by the following function
    psi_bc_func = @(x,eps,omega) control.K * exp(-control.ce .* (control.eps_0 - eps).^2) .* exp(-control.c0 .* (1 - omega).^2) .* (omega>=0);
    
    psi_bc = psi;
    x = permute(repmat(control.x, [1 control.energy_length control.mu_length]), [2 1 3]);
    eps = repmat(control.eps, [1 control.x_length control.mu_length]);
    mu  = permute(repmat(control.mu, [1 control.x_length control.energy_length]), [3 2 1]);
    
    % The electron fluence boundary conditions
    psi_bc_electron_0 = sum(psi_bc_func(0,eps,mu), 3) * control.delta_mu;
    psi_bc(:, 1, 3) = psi_bc_electron_0(:,1);
    psi_bc_electron_1 = sum(mu.* psi_bc_func(x,eps,mu), 3) * control.delta_mu; 
    psi_bc(:, 1, 4) = psi_bc_electron_1(:,1);
   
    psi_bc(:, end, 3) = 2 * control.delta * ones(size(psi_bc(:, end, 3)));
    psi_bc(:, end, 4) = zeros(size(psi_bc(:, end, 4)));
    
    psi_bc(:, 1, 1)   = zeros(size(psi_bc(:, 1, 1)));
    psi_bc(:, 1, 2)   = zeros(size(psi_bc(:, 1, 2)));
    psi_bc(:, end, 1) = zeros(size(psi_bc(:, end, 1)));
    psi_bc(:, end, 1) = zeros(size(psi_bc(:, end, 2)));

end

function R = computeR(psi, crossSection, control)
    R = zeros(size(psi));
    delta_eps = repmat(control.delta_eps, [1 control.x_length]);
    
    low_tri = tril(ones(control.energy_length));
    
	R1 = ((low_tri .* crossSection.Compton_photon_0) * squeeze(psi(:,:,1)) .* delta_eps);
    R((2:control.energy_length),:,1) = R1((1:(control.energy_length-1)),:);
    
    R2 = ((low_tri .* crossSection.Compton_photon_1) * squeeze(psi(:,:,2)) .* delta_eps);
    R((2:control.energy_length),:,2) = R2((1:(control.energy_length-1)),:);
    
    R3 = ((low_tri .* crossSection.Compton_electron_1) * squeeze(psi(:,:,1)) ...
        + (low_tri .* crossSection.Electrons_0)      * squeeze(psi(:,:,3))) .* delta_eps ...
        + crossSection.S ./ delta_eps .* squeeze(psi(:,:,3));
    R((2:(control.energy_length)),:,3) = R3((1:(control.energy_length-1)),:);
    
    R4 = ((low_tri .* crossSection.Compton_electron_1) * squeeze(psi(:,:,2)) ...
        + (low_tri .* crossSection.Electrons_1)      * squeeze(psi(:,:,4))) .* delta_eps ...
        + crossSection.S ./ delta_eps .* squeeze(psi(:,:,4));
    R((2:(control.energy_length)),:,4) = R4((1:(control.energy_length-1)),:);

    R = computeBoundaryCondition(R, control);
end

function F = compute_F(psi)
    psi0_g = psi(:,:,1);    
    psi1_g = psi(:,:,2);    
    psi0_e = psi(:,:,3);    
    psi1_e = psi(:,:,4);    
    F = zeros(size(psi));
    F(:,:,1) = psi1_g;
    F(:,:,3) = psi1_e;
    F(:,:,2) = compute_psi2(psi0_g, psi1_g);
    F(:,:,4) = compute_psi2(psi0_e, psi1_e);
end

function psi2 = compute_psi2(psi0, psi1)
    x = abs(psi1) ./ (psi0 + 1e-15);
    c0 = -0.0954823981432433; c1 = 0.229069986304953; c2 = -0.0344846229504588;
    theta1 = x.^2 + (2/3) * (1 - x.^2) + (x.^2) .* (1 - x.^2) .* (c0 + c1.*(x.^2) + c2.*(x.^4));
    chi2 = x.^2 .* theta1 + (1 - theta1);
    psi2 = psi0 .* chi2; %((1-chi2)./2 + (3*chi2 - 1)./2);% .* psi1.^2);
end

function [psi_new, D] = compute_psi(psi, crossSection, P, control)
    
    psi_new = zeros([control.energy_length control.x_length control.types]);
    delta_eps = control.delta_eps;
    
    a0 = crossSection.Compton_photon_Tot - primaryDiagonal(crossSection.Compton_photon_0) .* delta_eps;
    a1 = crossSection.Compton_photon_Tot - primaryDiagonal(crossSection.Compton_photon_1) .* delta_eps;
    b0 = -primaryDiagonal(crossSection.Compton_electron_0) .* delta_eps;
    b1 = -primaryDiagonal(crossSection.Compton_electron_1) .* delta_eps;
    c0 = crossSection.Electrons_Tot - primaryDiagonal(crossSection.Electrons_0) .* delta_eps;
    c1 = crossSection.Electrons_Tot - primaryDiagonal(crossSection.Electrons_1) .* delta_eps;
    
    a0 = repmat(a0, [1 control.x_length]); a1 = repmat(a1, [1 control.x_length]);
    b0 = repmat(b0, [1 control.x_length]); b1 = repmat(b1, [1 control.x_length]);
    c0 = repmat(c0, [1 control.x_length]); c1 = repmat(c1, [1 control.x_length]);
    
    delta_eps = repmat(control.delta_eps, [1 control.x_length]);
    c0 = c0 + crossSection.S ./ delta_eps;
    c1 = c1 + crossSection.S ./ delta_eps;

    rho = permute(repmat(control.rho, [1 control.energy_length]) , [2 1]);
    
    alpha0 = 1 ./ (1/control.delta_x + rho .* a0);
    alpha1 = 1 ./ (1/control.delta_x + rho .* a1);
    gamma0 = 1 ./ (1/control.delta_x + rho .* c0);
    gamma1 = 1 ./ (1/control.delta_x + rho .* c1);
    beta0  = -b0 .* alpha0 .* gamma0;
    beta1  = -b1 .* alpha1 .* gamma1;
    
    psi_new(:,:,1) = alpha0 .* P(:,:,1);
    psi_new(:,:,2) = alpha1 .* P(:,:,2);
    psi_new(:,:,3) = beta0  .* P(:,:,1) + gamma0 .* P(:,:,3);
    psi_new(:,:,4) = beta1  .* P(:,:,2) + gamma1 .* P(:,:,4);
    
    psi_new = computeBoundaryCondition(psi_new, control);
    
    D = psi / control.delta_x;
    D(:,:,1) = D(:,:,1) + a0 .* psi(:,:,1) .* rho;
    D(:,:,2) = D(:,:,2) + a1 .* psi(:,:,2) .* rho;
    D(:,:,3) = D(:,:,3) + (b0 .* psi(:,:,1) + c0 .* psi(:,:,3)) .* rho;
    D(:,:,4) = D(:,:,4) + (b1 .* psi(:,:,2) + c1 .* psi(:,:,4)) .* rho;
    
end

function pd = primaryDiagonal(M)
    pd = M(sub2ind(size(M),1:size(M,1),1:size(M,2)))';
end

function graphParams(ptitle, pxlabel, pylabel)
    grid on; title(ptitle); xlabel(pxlabel); ylabel(pylabel); set(gca, 'FontSize', 14); set(gcf,'color','w'); set(gca,'linewidth',2); legend('show');
end
