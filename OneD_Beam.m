classdef OneD_Beam
    properties
        x       % The depth axis
        dose    % The dose along this axis
    end
    
    methods
        function self = OneD_Beam()
           [self.x, self.dose] = doseComputation(); 
        end
        function plot(self)
            figure();
            plot(self.x, real(self.dose) ./ max(abs(real(self.dose))), 'DisplayName', 'Dose Theo'); hold on; 
            graphParams('Dose 1D beam', 'x', '$D(x)/D_{max}$', '', false);
        end
    end
end