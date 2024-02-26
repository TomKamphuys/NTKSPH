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

% size(vertical, horizontal)
% measurements size(m, 1) 1 is 'one'
% PSI_mat      size(m, c)
% coefs        size(c, 1) 1 is 'one


number = 31;
angle = linspace(-pi, pi, number)';
%angle = deg2rad([0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 180])'; % gedlee
r = 1.0.*ones(size(angle));%  + 0.1 .*randn(size(angle)); %rem((1:number)', 2); %randn(size(angle)); % m
Nmax = 15;


[x_meas, y_meas] = pol2cart(angle, r);


% Simulate a dipole source using 2 interior point sources.
%
pt_srcs = {[ 0.5,  0.4,  0.000, 0] ...
           [ 0.5, -0.4, -0.000, pi]};
%pt_srcs = {[ 0.5,  0.0,  0.0, 0]};


% To make things simpler, simulate 1 frequency only
freqs = 2000;
omega = 2*pi*freqs;

% Simulate the measurements
p_meas = sim_meas_cart_2d(x_meas, y_meas, omega(1), pt_srcs);

% Plot the measurements and the locations of the points sources
plotting = true;
if plotting
    fig1 = figure('Name', 'Simulated Measurements and Reconstructions without SFS', ...
                  'Position', [100 190 1600 700]);
    ax1 = subplot(1, 6, [1 4]);
    for n = 1:length(pt_srcs)
        scatter(ax1, pt_srcs{n}(2), pt_srcs{n}(3), ...
                 60, 'k', 'p', 'MarkerFaceColor', 'k');
        if n == 1
            hold on;
        end
    end
    scatter(ax1, x_meas, y_meas, 10, dB_SPL(p_meas), 'filled');
    axis(ax1, 'equal');
    colormap('jet');
%    caxis([85, 95]);
    colorbar(ax1, 'eastoutside');
    title(ax1, 'Simulated Measurements and Reconstructions without SFS');
end


% Generate the measurement PSI matrix
[PSI_mat, Nmax] = sph_PSI_mix_2d(r, angle, omega, Nmax, temp);
disp(cond(PSI_mat));


%%
%{
   Step 2 - Compute the fitting coefficients
    The column vector 'CD_vec' are the fitting coefficients
%}
[CD_vec, res] = lstsq_solve(PSI_mat, p_meas);
%CD_vec = -.13*[0 i]';
%disp(CD_vec);

%%
%{
   Step 3 - Reconstruction grid and the reconstruction PSI matrix
    We will use points on a vertical and a horizontal line, forming a
     cross pattern that is parallel to the YZ plane.
%}

amount = 201;
r_recon = 1.*ones(1,amount)';
a_recon = linspace(-pi, pi, amount)';

[x_recon, y_recon] = pol2cart(a_recon, r_recon);

PSI_recon = sph_PSI_mix_2d(r_recon, a_recon, omega, Nmax, temp);


%   Step 4 - Reconstruction without soundfield separation
p_recon = PSI_recon * CD_vec;


int_srcs = pt_srcs;
p_ref = sim_meas_cart_2d(x_recon, y_recon, omega(1), int_srcs);

p_err = p_ref - p_recon;
ax2 = subplot(1, 6, [5 6]);
plot(ax2, 1:length(p_err), dB_SPL(p_recon), 'color', 'b');
hold on;
plot(ax2, 1:length(p_err), dB_SPL(p_ref), 'color', 'g');
plot(ax2, 1:length(p_err), dB_SPL(p_err), 'color', 'r');
pbaspect(ax2, [1 1 1]);
ax2_title = sprintf('Reconstruction vs Reference, Nmax = %d', Nmax);
title(ax2, ax2_title);
ylabel(ax2, 'dB SPL');
legend(ax2, 'Reconstruction', 'Reference', 'Error', ...
       'Location', 'southoutside');
hold off;

figure
subplot(2, 2, 1)
[X, Y, outRef, outRecon] = plane_reconstruct_2d(2, 2, 100, pt_srcs, CD_vec, omega, Nmax, temp);
pcolor(X, Y, outRef)
shading flat
title('Reference')
hold on
plot(x_meas, y_meas, 'k.')
axis equal
colorbar

minimum = min(outRef(:));
maximum = max(outRef(:));

subplot(2, 2, 2)
pcolor(X, Y, outRecon)
shading flat
title('Reconstructed' )
axis equal
caxis([minimum, maximum])
colorbar

subplot(2, 2, 3)
bar( abs(CD_vec));

subplot(2, 2, 4)
pcolor(X, Y, outRecon - outRef)
shading flat
title('Difference')
axis equal
caxis([0, maximum/5])
colorbar
