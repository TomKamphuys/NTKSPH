%% ----- Simulate Measurements -----
%{
  Initial release: 2021-04-23

  Create some simulated measurements.
  Five different types of measurement grids are avaliable.
   'planar' - Parallel surfaces
   'full_cyl' - Concentric full cylinders
   'half_cyl' - Concentric half cylinders
   'full_sph' - Concentric full spheres
   'hald_sph' - Concentric half spheres
%}
%% Parameters: Grid, frequencies, sources, save figures or not
clear all
close all
% Save the simulated measurements plots or not?
save_figs = false;

% Valid grid types are: 'planar', 'full_cyl','full_cylTopBot', 'half_cyl', 'full_sph',
%  and 'half_sph'.
grid_type = 'full_cylTopBot';  

% List of frequencies to simulate.
% The following example gives: 20 to 20480 Hz in log scale (10 octaves),
%  2 frequencies per octave, totalling 21 frequencies.

freqs = logspace(log10(20), log10(20*(2^10)), 21);
freqs = logspace(log10(500), log10(2000), 3);
freqs = [50 1000 ];

% Save frequencies of the simulated measurments to 'frequencies.cvs'
csvwrite('frequencies.csv', freqs');

% The point sources: Each element of the cell array is one point source.
% Format is: [ Amplitude, x, y, z, phase in radians ]
pt_srcs = {
%            [  .5,  0.6  0.25,  0,   0] ...
%            [ -.5,  0.43,  0.25,  -0.0,  pi/2] ...

            [  .1,  0.2,  0.4,  -0.4,  0.0] ...
            [  .1,  0.25,  0.4,  -0.3,  0.0] ...
            [  .1,  0.3,  0.4,  -0.2,  0.0] ...
            [  .1,  0.35,  0.4,  -0.1,  0.0] ...
            [  .1,  0.4,  0.4,  0.0,  0.0] ...
            [  .1,  0.45,  0.4,  0.1,  0.0] ...
            [  .1,  0.5,  0.4,  0.2,  0.0] ...
            [  .1,  0.55,  0.4,  0.3,  0.0] ...
            [  .1,  0.6,  0.4,  0.4,  0.0] ...
     
           }; 

% Parameters for 'planar' grid (parallel planar measurement surfaces)
y_range_pln = [-0.5 0.5];
y_pnts_pln = 25;
z_range_pln = [-0.50 0.50];
z_pnts_pln = 25;
x_pln = [0.25 0.30];

% Parameters for 'full_cyl' or 'half_cyl' grid (concentric cylindrical
%  measurement surfaces, either all around or one sided)
z_range_cyl = [0 0];
z_pnts_cyl = 1;
th_pnts_cyl = 360/2;
radii_cyl = [3.4  ];
r_pnts_cap=0;

%=================================================================================================================================
% Parameters for 'full_cylTopBot'(concentric cylindrical
%  measurement surfaces, either all around with a flat round top and bottom cap at the postition of the lowest and highest z-value and a second cap at the same distance as defined at radii_cyl)
z_range_cyl = [-0.75 0.75];
z_pnts_cyl = 16;
th_pnts_cyl = 360/5;
radii_cyl = [1 1.05];
r_pnts_cap=8;
%=================================================================================================================================

% Parameters for 'full_sph' or 'half_sph' grid (concentric spherical
%  measurement surfaces, either all around or one sided)
n_approx_sph = 1250;
radii_sph = [0.45 0.50];

%% Generate the simulated measurement grid according to the grid type
%
if strcmpi(grid_type, 'planar')
    % Double parallel planar grid (parallel to YZ plane) 
    [x_meas, y_meas, z_meas] = planar(y_range_pln, y_pnts_pln, ...
                                      z_range_pln, z_pnts_pln, x_pln);
elseif strcmpi(grid_type, 'full_cyl')
    % Concentric full cylindrical (full circle) grid
    [x_meas, y_meas, z_meas] = conc_full_cyl(z_range_cyl, z_pnts_cyl, ...
                                             th_pnts_cyl, radii_cyl);                                       
elseif strcmpi(grid_type, 'full_cylTopBot')
    % Concentric full cylindrical (full circle) grid
    [x_meas, y_meas, z_meas] = conc_full_cylTopBot(z_range_cyl, z_pnts_cyl, ...
                                             th_pnts_cyl, radii_cyl,r_pnts_cap);
elseif strcmpi(grid_type, 'half_cyl')
    % Concentric half cylindrical (half circle) grid
    [x_meas, y_meas, z_meas] = conc_half_cyl(z_range_cyl, z_pnts_cyl, ...
                                             th_pnts_cyl, radii_cyl);
elseif strcmpi(grid_type, 'full_sph')
    % Concentric full spherical (all around) grid
    [x_meas, y_meas, z_meas] = conc_full_sph(n_approx_sph, radii_sph);
elseif strcmpi(grid_type, 'half_sph')
    % Concentric half spherical (one sided) grid
    [x_meas, y_meas, z_meas] = conc_half_sph(n_approx_sph, radii_sph);
else
error('%s is not a valid measurement',grid_type);
end

% Save measurement point coordinates to 'coordinates.csv'
csvwrite('coordinates.csv', [x_meas y_meas z_meas]);

%% Generate two disc shaped plotting grids on XY and XZ planes
%  These 2 grids are for visualizing the sound field radiated by
%   the point sources.
[x_d1, y_d1, z_d1] = cart_disc_grid('xz');
[x_d2, y_d2, z_d2] = cart_disc_grid('xy');

%% Generate and plot the simulated measurements

% Allocate storage to store the measurements
p_meas = zeros(length(x_meas), length(freqs));

fig1 = figure('Name', 'Simulated Measurements', ...
              'Position', [100 190 850 700]);
ax1 = subplot(1, 1, 1);

for m = 1:length(freqs)
    omega = 2 * pi * freqs(m);
    
    % Display the point sources as black 5 point stars
    for n = 1:length(pt_srcs)
        scatter3(ax1, pt_srcs{n}(2), pt_srcs{n}(3), ...
                 pt_srcs{n}(4), ...
                 60, 'k', 'p', 'MarkerFaceColor', 'k');
        if n == 1
            hold on;
        end
    end
    
    % Plot the XZ disc grid
    p_d1 = sim_meas_cart(x_d1, y_d1, z_d1, omega, pt_srcs);
    p_d1plot =20*log10(abs(p_d1)); %plot in dB values
    
    surf(ax1, x_d1, y_d1, z_d1, (p_d1plot), ...
         'FaceColor', 'interp', 'EdgeColor', 'none');
    if m == 1
        colormap('jet');
        axlims = plot_extent(x_d1, y_d1, z_d1);
        xlim(ax1, [-axlims axlims]);
        ylim(ax1, [-axlims axlims]);
        zlim(ax1, [-axlims axlims]);
    end
    axis(ax1, 'equal');
%     caxis([0 2.5]);
     caxis([-40 6]);
    colorbar(ax1, 'eastoutside');
    view(ax1, [130 20]);
    
    % Plot the XY disc grid
    p_d2 = sim_meas_cart(x_d2, y_d2, z_d2, omega, pt_srcs);
    p_d2plot =20*log10(abs(p_d2)); %plot in dB values

    surf(ax1, x_d2, y_d2, z_d2, p_d2plot, ...
         'FaceColor', 'interp', 'EdgeColor', 'none');
    
    % Plot the measurement points
    p_meas(:, m) = sim_meas_cart(x_meas, y_meas, z_meas, omega, pt_srcs);
%     scatter3(ax1, x_meas, y_meas, z_meas, ...
%              abs(p_meas(:, m)), abs(p_meas(:, m)));

    scatter3(ax1, x_meas, y_meas, z_meas,10,[0.8 0.8 0.8],'filled');
          
    title_str = sprintf('f = %.0f Hz', freqs(m));
    title(ax1, title_str);
    drawnow;
    if save_figs
        fname = sprintf('sim_f%.0f.png', freqs(m));
        saveas(fig1, fname, 'png');
    end
    pause(0.5);
    hold off;
end

% Save simulated measurements to 'measurements.csv'.
csvwrite('measurements.csv', p_meas);

main