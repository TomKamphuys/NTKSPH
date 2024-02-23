close all

%{
   An example script to illustrate the process of calculating the fitting
    coefficients, and how to use these fitting coefficients to reconstruct
    the sound pressure field at arbitrary locations.

   No soundfield separation is used for the reconstruction in this example.

   To simplify matters, this script perform the fitting and reconstruction
    at 1 frequency only. Also, it is assumed that the optimal expansion
    order is known already.

   For a 2D plane only trying to replicate what I think is Gedlee's method
%}
%%
%{
   Step 0 - Generate the simulated measurements
    We generate the measurement grid first.
    Then we simulate the measurements at the grid points.
%}

% Temperature is 20 deg C.
temp = 273.15 + 20;

% In this example I'm using Kessito's 'full_cylTopBot' grid.
%z_range_cyl = [-0.75 0.75];
%z_pnts_cyl = 16;
%th_pnts_cyl = 360/5;
%radii_cyl = [1 1.05];
%r_pnts_cap=8;


number = 201;
angle = linspace(-pi, pi, number)';
r = 1.0;% + 0.5.*randn(size(angle)); %rem((1:number)', 2); %randn(size(angle)); % m
Nmax = 6;


x_meas = r.*cos(angle);
y_meas = zeros(size(x_meas));
z_meas = r.*sin(angle);

%[x_meas, y_meas, z_meas] = conc_full_cylTopBot(z_range_cyl, z_pnts_cyl, ...
%                                        th_pnts_cyl, radii_cyl,r_pnts_cap);

% Simulate a dipole source using 2 interior point sources.
%
pt_srcs = {[ 0.5,  0.05  0.0, 0.0, 0] ...
           [ 0.5, -0.05, 0.0, 0.0, 0]};

%pt_srcs = {[ 1.0, 0.0, 0.0, 0.0, 0.0]};



% To make things simpler, simulate 1 frequency only
freqs = 2000;
omega = 2*pi*freqs;

% Simulate the measurements
p_meas = sim_meas_cart(x_meas, y_meas, z_meas, omega(1), pt_srcs);

% Plot the measurements and the locations of the points sources
plotting = true;
if plotting
    fig1 = figure('Name', 'Simulated Measurements and Reconstructions without SFS', ...
                  'Position', [100 190 1600 700]);
    ax1 = subplot(1, 6, [1 4]);
    for n = 1:length(pt_srcs)
        scatter3(ax1, pt_srcs{n}(2), pt_srcs{n}(3), ...
                 pt_srcs{n}(4), ...
                 60, 'k', 'p', 'MarkerFaceColor', 'k');
        if n == 1
            hold on;
        end
    end
    scatter3(ax1, x_meas, y_meas, z_meas, 10, dB_SPL(p_meas), 'filled');
    axis(ax1, 'equal');
    view(ax1, [130 20]);
    colormap('jet');
    caxis([85, 95]);
    colorbar(ax1, 'eastoutside');
    title(ax1, 'Simulated Measurements and Reconstructions without SFS');
end
%%
%{
   Step 1 - Generate the measurement PSI matrix
    For simplicity, we assume that we have already determined that the
    optimal expansion order N is 12.

   Note: Because typically the measurement point coordinates are read from
         CSV files, they are read in as 1 by m row vectors. The functions
         to generate the PSI matrix expects the coordinates to be in
         1 by m row vectors.
         However, the functions to generate the measurement grids return
         m by 1 column vectors. Therefore, we need to transpose the
         measurement coordinates 'x_meas', 'y_meas' and 'z_meas' column
         vectors into row vectors!
%}
x_meas = x_meas';
y_meas = y_meas';
z_meas = z_meas';

% Assume we already know the optimal expansion order is 14.
%
%Nmax = 14;

% Convert the measurement points coordinates from Cartesian to polar
[phi_meas, theta_meas, r_meas] = cart2sph (x_meas, y_meas, z_meas);
theta_meas = pi/2 - theta_meas;

% Generate the measurement PSI matrix
[PSI_mat, Nmax] = sph_PSI_mix_2d(r_meas, theta_meas, phi_meas, omega, Nmax, temp);

%%
%{
   Step 2 - Compute the fitting coefficients
    The column vector 'CD_vec' are the fitting coefficients
%}
[CD_vec, res] = lstsq_solve(PSI_mat, p_meas);

%%
%{
   Step 3 - Reconstruction grid and the reconstruction PSI matrix
    We will use points on a vertical and a horizontal line, forming a
     cross pattern that is parallel to the YZ plane.
%}
y_hor = linspace(-0.75, 0.75, 21);
z_hor = zeros(1, 21);

y_ver = zeros(1, 21);
z_ver = linspace(-0.75, 0.75, 21);

% Remember: 'x_recon', 'y_recon', 'z_recon' must be 1 by m row vectors!
x_recon = repmat(1.1, 1, 21);
y_recon = [z_hor];
z_recon = [y_hor];



% Convert the reconstruction points coordinates from Cartesian to polar
[phi_recon, theta_recon, r_recon] = cart2sph (x_recon, y_recon, z_recon);
theta_recon = pi/2 - theta_recon;

PSI_recon = sph_PSI_mix_2d(r_recon, theta_recon, phi_recon, omega, Nmax, temp);

%%
%{
   Step 4 - Reconstruction without soundfield separation


%}
p_recon = PSI_recon * CD_vec;


% Plot the reconstruction
if plotting
    scatter3(ax1, x_recon, y_recon, z_recon, 60, dB_SPL(p_recon), '+', ...
             'linewidth', 2);
    hold off;
end

%%
%{
   Verify the accuracy of the SFS reconstruction.
%}

% Calculate the reference sound pressures at the reconstruction points
%  for the interior sources only.
%int_srcs = {[ 0.5,  0.05  0.0, 0.0, 0] ...
%            [ 0.5, -0.05, 0.0, 0.0, pi]};
int_srcs = pt_srcs;
p_ref = sim_meas_cart(x_recon, y_recon, z_recon, omega(1), int_srcs);

p_err = p_ref - reshape(p_recon, 1, []);
ax2 = subplot(1, 6, [5 6]);
line(ax2, 1:length(p_err), dB_SPL(p_recon), 'color', 'b');
hold on;
line(ax2, 1:length(p_err), dB_SPL(p_ref), 'color', 'g');
line(ax2, 1:length(p_err), dB_SPL(p_err), 'color', 'r');
pbaspect(ax2, [1 1 1]);
ax2_title = sprintf('Reconstruction vs Reference, Nmax = %d', Nmax);
title(ax2, ax2_title);
ylabel(ax2, 'dB SPL');
legend(ax2, 'Reconstruction', 'Reference', 'Error', ...
       'Location', 'southoutside');
hold off;

figure

[X, Y, outRef, outRecon] = plane_reconstruct_2d(2, 2, 200, pt_srcs, CD_vec, omega, Nmax, temp);
pcolor(X, Y, outRef)
shading flat
title('Reference')
hold on
plot(x_meas, z_meas, 'k*')
axis equal
colorbar

minimum = min(outRef(:));
maximum = max(outRef(:));

figure
pcolor(X, Y, outRecon)
shading flat
title('Reconstructed' )
axis equal
caxis([minimum, maximum])
colorbar
