classdef TwoD_PencilBeam
    properties(Access = public)
        x           % The first axis
        y           % The second axis
        dose        % The dose along those axis
        width       % Width of the beam [cm]
        location    % Initial location of the beam [x_ini,y_ini]
        angle       % Oriontation of the beam [cos(angle), sin(angle)]
        boundary    % Boundry on which the beam is located [1-4]
        oneD_beam   % The 1D beam used for the pencil beam technique
    end
    
    methods(Access = public)
        function self = TwoD_PencilBeam(oneD_beam, width, location, angle, boundary)
           self.oneD_beam = oneD_beam; self.width = width; self.location = location; self.angle = angle; self.boundary = boundary;
           self.x = oneD_beam.x; self.y = oneD_beam.x;
           self.dose = twoD_PencilBeamdoseComputation(self); 
        end
        function plot(self)
            figure();
            X = repmat(self.x,  [1 length(self.y)]);
            Y = repmat(self.y' - 17.5, [length(self.x) 1]);
            contourf(X, Y, real(self.dose) ./ max(abs(real(self.dose)), [], 'all'), 'DisplayName', 'Dose Theo'); hold on; 
            graphParams('Dose  2D', 'x', 'y', '$D(x)/D_{max}$', true);
        end
    end
    methods(Access = private)
        function dose = twoD_PencilBeamdoseComputation(obj)          
            dose = zeros(length(obj.x));
            dose(end/2,:) = obj.oneD_beam.dose;
            base_dose = dose;
            angular_spread = 100e-3/2.3; % Conversion from FWHM to STD
            effective_angles = linspace(- 4*angular_spread, 4*angular_spread, 61);
            for i = 1 : length(effective_angles)
                dose = dose + rotate(obj, base_dose, effective_angles(i))*exp(-effective_angles(i)^2/(2*(angular_spread)^2));
            end
            dose = rotate(obj, dose, obj.angle);
            dose = expand(obj, dose, obj.width);
            dose = medfilt2(dose,[5 5]);
            dose = reflect(obj, dose, obj.boundary);
            dose = translate(obj, dose, obj.location);
        end
        function reflection = reflect(obj, dose, boundary)
            switch boundary
                case 1
                    reflection = permute(dose, [2 1]);
                    reflection = fliplr(reflection);
                case 2
                    reflection = fliplr(dose);
                    reflection = flipud(reflection);
                case 3
                    reflection = permute(dose, [2 1]);
                    reflection = flipud(reflection);
                otherwise
                    reflection = dose;
            end

        end
        function dose = translate(obj, dose, location)
            if (obj.boundary == 1) | (obj.boundary == 3)
                loc_x = 0; loc_y = location;
            else
                loc_x = location; loc_y = 0;
            end
            dose = shift(obj, dose, floor(loc_x / max(obj.x) * length(obj.x)), floor(loc_y / max(obj.y) * length(obj.y)));
        end
        function expansion = expand(obj, mat, width)
            expansion = mat;
            for i = 1 : (floor(width/max(obj.x) * length(obj.x))/2)
                expansion = expansion + (shift(obj, mat, i, 0) + shift(obj, mat, -i, 0));%*exp(-i^2/100);
            end            
        end
        function rot = rotate(obj, mat, angle)
            rot = zeros(size(mat));
            vec = (1:(length(obj.x) * length(obj.y)));
            r_vec = rotate_aux(obj, vec', angle);
            rot(r_vec) = mat(vec);
        end
        
        function rot_loc = rotate_aux(obj, loc, angle)
            R_mat = [cos(angle) -sin(angle); sin(angle) cos(angle)];
            lambda_x = 1;
            lambda_y = 0;
            x_cor = mod(loc, length(obj.x)) - lambda_x*length(obj.x)/2;
            y_cor = floor(loc / length(obj.y)) - lambda_y*length(obj.y)/2;
            cor = [x_cor,y_cor]';
            rot_cor = R_mat*cor;
            rot_cor(1,:) = floor(rot_cor(1,:) + lambda_x*length(obj.x)/2);
            rot_cor(2,:) = floor(rot_cor(2,:) + lambda_y*length(obj.y)/2);
                      
            mask_x = (rot_cor(1,:) >= 1) & (rot_cor(1,:) < length(obj.x));
            mask_y = (rot_cor(2,:) >= 1) & (rot_cor(2,:) < length(obj.y));
            mask   = mask_x & mask_y;
            
            rot_loc = (rot_cor(2,:) * length(obj.y) + rot_cor(1,:)).* mask + (1-mask);
        end

        function shifted_mat = shift(obj, mat, shift_x, shift_y)
            if shift_x >= 0 
                shift_x = shift_x + 1;
            else
                shift_x = shift_x - 1;
            end
            if shift_y >= 0 
                shift_y = shift_y + 1;
            else
                shift_y = shift_y - 1;
            end

            shifted_mat = zeros(size(mat));
            if ((shift_x >= 0) && (shift_y >= 0))
                shifted_mat((shift_x:end),(shift_y:end)) = mat((1:(end - shift_x + 1)),(1:(end - shift_y + 1)));
            elseif ((shift_x >= 0) && (shift_y < 0))
                shift_y = -shift_y;
                shifted_mat((shift_x:end),(1:(end - shift_y + 1))) = mat((1:(end - shift_x + 1)),(shift_y:end));
            elseif ((shift_x < 0) && (shift_y >= 0))
                shift_x = -shift_x;
                shifted_mat((1:(end - shift_x + 1)),(shift_y:end)) = mat((shift_x:end),(1:(end - shift_y + 1)));
            else
                shift_x = -shift_x; shift_y = -shift_y;
                shifted_mat((1:(end - shift_x + 1)),(1:(end - shift_y + 1))) = mat((shift_x:end),(shift_y:end));
            end
        end

    end
end