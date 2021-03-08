classdef TwoD_Array
    properties(Access = public)
        x           % The first axis
        y           % The second axis
        num_beams   % Number of beams
        dose        % The dose in the material
        widths      % Width of the beams [cm]
        locations   % Initial location of the beams [delta]
        angles      % Oriontation of the beams [cos(angle), sin(angle)]
        boundaries  % Boundry on which the beams are located [1-4]
        focused_beam   % The 1D beam used for the pencil beam
        twoD_beams
    end
    
    methods(Access = public)
       function self = TwoD_Array(focused_beam, widths, locations, angles, boundaries)
           self.widths = widths; self.locations = locations; self.angles = angles; self.boundaries = boundaries;
           self.num_beams = length(locations);
           self.focused_beam = focused_beam;
           self.x = focused_beam.x; self.y = focused_beam.x;
           self.dose = zeros([length(self.x), length(self.y)]);
           self.twoD_beams = {};
           for i = 1 : self.num_beams
                self.dose = self.dose + focusedBeamTransform(self.focused_beam, locations(i), angles(i), boundaries(i));
           end
       end
       function plot(self)
            figure();
            X = repmat(self.x,  [1 length(self.y)]);
            Y = repmat(self.y' - 17.5, [length(self.x) 1]);
            levels = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50];
            contourf(X, Y, real(self.dose) ./ max(abs(real(self.dose)), [], 'all'), 'DisplayName', 'Dose Theo'); hold on; 
            graphParams('Dose  2D', 'x', 'y', '$D(x)/D_{max}$', true);
%             ylim([-2 2]); xlim([0 40]);
        end
    end
end