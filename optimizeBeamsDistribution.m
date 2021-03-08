function [optimized_angles, optimized_loc, optimized_boundaries] = optimizeBeamsDistribution(beam, num_beams, mask)
    opts = optimoptions('fmincon', ...
    'OptimalityTolerance', 0, ...
    'StepTolerance', 0, ...
    'MaxFunctionEvaluations', inf,...
    'MaxIterations', 1e6);

    max_angle = pi/6;
    max_loc   = 5;%max(beam.x) / 4;
    
    angle_lb =  -max_angle * ones(num_beams, 1);
    loc_lb    = -max_loc  * ones(num_beams, 1);
    lb = [angle_lb];% loc_lb];
    
    angle_ub = max_angle * ones(num_beams, 1);
    loc_ub   = max_loc   * ones(num_beams, 1);
    ub = [angle_ub];% loc_ub];
    
    fixed_locations = reshape(repmat(linspace(-max_loc, max_loc, 11), [4 1]), 1, 44)';
    
    boundaries = repmat((1:4)', [num_beams 1]); boundaries = boundaries(1:num_beams);
    
    entropy = @(A,B) - sum(1 / sum(B,'all') .* log((A / sum(A, 'all') + 1e-10)), 'all');
    mean_nz = @(A) sum(A, 'all') / sum(A ~= 0, 'all');
%     loss_func = @(A,B) sum((A - A.*B).^2 - abs(A .* log(A + 1e-3)) .* (1 - B), 'all');
    loss_env = @(A,B) sqrt(sum((A .* (1 - B) - mean_nz(A .* (1 - B))).^2 / mean_nz(A .* (1 - B)).^2, 'all'));
    
    lambda1 = 40; lambda2 = 10;
    loss_func = @(A,B) (lambda1 * sum(-(A.*B).^2, 'all') + lambda2 * sum((A.*(1-B)).^2, 'all'))/sum(A,'all');%...
%           - 0.5*entropy(A(B==1), B) - 0.5*entropy(A(B==0), 1-B);
%     loss = @(param) loss_func(TwoD_Array(beam, 0, param((num_beams+1):(2*num_beams)), param(1:num_beams), boundaries).dose, mask);
    loss = @(param) loss_func(TwoD_Array(beam, 0, fixed_locations, param, boundaries).dose, mask);
    
    optimized_boundaries = boundaries;
    optimized_loc = fixed_locations;
    num_iterations = 1;
    min_val = 1e80;
    for i = 1:num_iterations
%         angles_sp   = max_angle * rand(num_beams, 1) - max_angle/2;
%         angles_sp   = atan((beam.focal_length - 2.5) ./ (fixed_locations));
        angles_sp   = -(mod(atan((beam.focal_length + 5*2.5) ./ fixed_locations), pi) - pi/2);
%         loc_sp      = max_loc * rand(num_beams, 1) - max_loc;

        param_sp = [angles_sp];% loc_sp];
        [optimized_param,min_val_cur,exitflag1,output_min] = fmincon(loss, param_sp, [], [], [], [], lb, ub, [], opts);
        min_val_cur
        fprintf('Iteration = %d out of %d, value: %f \n', i, num_iterations, min_val_cur);
        if(min_val_cur < min_val)
            min_val = min_val_cur;
%             optimized_angles = optimized_param(1:num_beams);
%             optimized_loc    = optimized_param((num_beams+1):(2*num_beams));
            optimized_angles = optimized_param;
        end
    end

end