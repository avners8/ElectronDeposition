classdef Focused_Beam
    properties(Access = public)
        x           % The first axis
        y           % The second axis
        dose        % The dose in the material
        width       % Width of the beam [cm]
        location    % Initial location of the beam [delta]
        angle       % Oriontation of the beam [cos(angle), sin(angle)]
        boundary    % Boundry on which the beams are located [1-4]
        oneD_beam   % The 1D beam used for the pencil beam
        focal_length 
    end
    
    methods(Access = public)
       function self = Focused_Beam(oneD_beam, width, location, angle, boundary, focal_length)
           self.width = width; self.location = location; self.angle = angle; self.boundary = boundary;
           self.oneD_beam = oneD_beam; self.focal_length = focal_length;
           self.x = oneD_beam.x; self.y = oneD_beam.x;
           self.dose = computeFocusedBeamDose(self);
       end
       function dose = computeFocusedBeamDose(self)
           dose = zeros([length(self.x), length(self.y)]);
           internal_width = 1;
           width_vec = linspace(-self.width / 2, self.width / 2, 21);
           for i = 1 : length(width_vec) 
               internal_angle = mod(atan(self.focal_length ./ width_vec(i)), pi) - pi/2;
               twoD_beam = TwoD_PencilBeam(self.oneD_beam, internal_width, self.location + width_vec(i), internal_angle, self.boundary); 
               dose = dose + twoD_beam.dose ./ max(abs(twoD_beam.dose), [], 'all');
           end
       end
       function dose = focusedBeamTransform(self, loc, angle, boundary)
           dose = reflect(self, self.dose, 2);
           dose = rotate(self, dose, angle);
           dose = medfilt2(dose,[5 5]);
           dose = circshift(dose,floor(loc / (max(self.x)) * length(self.x)),2);
           dose = reflect(self, dose, boundary);
       end
       function reflection = reflect(self, dose, boundary)
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
       function rot = rotate(obj, mat, angle)
            rot = zeros(size(mat));
            vec = (1:(length(obj.x) * length(obj.y)));
            r_vec = rotate_aux(obj, vec', angle);
            rot(r_vec) = mat(vec);
       end
        
       function rot_loc = rotate_aux(obj, loc, angle)
            R_mat = [cos(angle) -sin(angle); sin(angle) cos(angle)];
            x_cor = mod(loc, length(obj.x)) - 2*length(obj.x)/2;
            y_cor = floor(loc / length(obj.y)) - 1*length(obj.y)/2;
            cor = [x_cor,y_cor]';
            rot_cor = R_mat*cor;
            rot_cor(1,:) = floor(rot_cor(1,:) + 2*length(obj.x)/2);
            rot_cor(2,:) = floor(rot_cor(2,:) + 1*length(obj.y)/2);
                      
            mask_x = (rot_cor(1,:) >= 1) & (rot_cor(1,:) < length(obj.x));
            mask_y = (rot_cor(2,:) >= 1) & (rot_cor(2,:) < length(obj.y));
            mask   = mask_x & mask_y;
            
            rot_loc = (rot_cor(2,:) * length(obj.y) + rot_cor(1,:)).* mask + (1-mask);
       end
       function plot(self)
            figure();
            X = repmat(self.x,  [1 length(self.y)]);
            Y = repmat(self.y' - max(self.x)/2, [length(self.x) 1]);
            levels = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50];
            contourf(X, Y, real(self.dose) ./ max(abs(real(self.dose)), [], 'all'), 'DisplayName', 'Dose Theo'); hold on;             
%             contour(X, Y, 50*real(self.dose) ./ max(abs(real(self.dose)), [], 'all'), levels, 'DisplayName', 'Dose Theo'); hold on; 
            graphParams('Dose  2D', 'x', 'y', '$D(x)/D_{max}$', true);
%             ylim([-2 2]); xlim([0 40]);
        end
    end
end