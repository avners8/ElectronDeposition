%% main
set(groot, 'defaultFigurePosition', [100 100 900 600]); set(groot, 'defaultTextInterpreter', 'latex'); set(groot, 'defaultLegendInterpreter', 'latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'DefaultLineLineWidth', 4);

%% 1D Dose optimization

oneD_beam = OneD_Beam(); 
plot(oneD_beam);

%% Focused Beam

width      = 1;
location   = 0; angle = 0; boundary = 1;
focal_length = 17.5;

beam = Focused_Beam(oneD_beam, width, location, angle, boundary, focal_length);
plot(beam);

%% Optimization

% Constructing the mask
num_beams = 44;
[X,Y] = meshgrid(beam.x, beam.y);
mask_x = (X >= max(beam.x)/2 - 5/2) .* (X <= max(beam.x)/2 + 5/2);
mask_y = (Y >= max(beam.y)/2 - 5/2) .* (Y <= max(beam.y)/2 + 5/2);
mask = mask_x .* mask_y;

% Optimizing
[optimized_angles, optimized_loc, optimized_boundaries] = optimizeBeamsDistribution(beam, num_beams, mask);

% Plotting
beams = TwoD_Array(beam, width, optimized_loc, optimized_angles, optimized_boundaries);
plot(beams);
rectangle('Position',[15 -2.5 5 5], 'EdgeColor', 'w', 'LineWidth', 3);

theDose = medfilt2(beams.dose,[5 5]);
energy_ratio = sum(beams.dose .* mask, 'all') / sum(beams.dose,'all') * 100;
peak_ratio = (max(theDose(mask==1), [], 'all') - mean(theDose(mask==1), 'all')) / (mean(theDose(mask==1), 'all') - min(theDose(mask == 1), [], 'all')) * 100;
fprintf('Energy ratio: %0f, Peak to valley ratio: %0f\n', energy_ratio, peak_ratio);

%%
rep = 11;
location   = repmat(linspace(-2.1, 2.1, rep), [1 4])'; 

boundary   = repelem([1 2 3 4], rep)';
angle = zeros(size(location));
widths = 0.8 * ones(size(location));
theBeam = TwoD_Array(beam, widths, location, angle, boundary);
plot(theBeam);
rectangle('Position',[15 -2.5 5 5], 'EdgeColor', 'k', 'LineWidth', 3);

theDose = medfilt2(theBeam.dose,[5 5]);
energy_ratio = sum(theBeam.dose .* mask, 'all') / sum(theBeam.dose,'all') * 100;
peak_ratio = (max(theDose(mask==1), [], 'all') - mean(theDose(mask==1), 'all')) / (mean(theDose(mask==1), 'all') - min(theDose(mask == 1), [], 'all')) * 100;
fprintf('Energy ratio: %0f, Peak to valley ratio: %0f\n', energy_ratio, peak_ratio);
