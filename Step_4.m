%% ----- Step 4 -----
%{
  Initial release: 2021-04-18
  Updated:         2021-04-22

  Use the measurement PSI matrix and measured pressure to compute the
   spherical wave expansion function coefficients.
  For each maximum order (i.e. n = 0, 1, 2, ..., Nmax), we need to
   compute the expansion coefficients separately.
  The residuals (mean square error) returned by the least squares
   calculations is just one indicator of the fitting accuracy.
  Convergence of reconstructions at the CTA-2034 grid is also used as
   a separate (and more reliable) indicator of the fitting accuracy.
%}
% Start with reconstruction with order n = 0; pre-allocate storage for
%  arrays

% xoffset=x-offsetposition(1,1);
% yoffset=y-offsetposition(1,2);
% zoffset=z-offsetposition(1,3);

n = 0;
residuals = zeros(1, Nmax+1); % Array to store fitting residuals
conv = zeros(1, Nmax); % Array to store convergence metric
p2034 = zeros(p2034_len, Nmax + 1); % Store CTA2034 reconstructions

% Solve for the coefficients using least sqaures: coeff = PSI \ press
tic;
[CD_vec, res] = lstsq_solve(PSI_mat(:, 1:2*(n+1)^2), p_meas(:, f_bin));
fprintf('N = %d, ', n)
fprintf('f = %.0f Hz, compute time = %.3f; \n', freqs(f_bin), toc)
residuals(1) = res;

% Reconstruct at the 2 plotting grids p1 & p2. The reconstruction
%  results are reshaped back to the same matrix dimensions as the
%  plotting grid coordinates for surface plotting.
p1 = PSI_p1(:, 1:2*(n+1)^2) * CD_vec;
p1 = reshape(p1, p1_shape);
p1plot =20*log10(abs(p_d1)); %plot in dB values

p2 = PSI_p2(:, 1:2*(n+1)^2) * CD_vec;
p2 = reshape(p2, p2_shape);
p2plot =20*log10(abs(p_d1)); %plot in dB values

% Reconstruct at the CTA2034 spinorama grid
p2034(:, 1) = PSI_2034(:, 1:2*(n+1)^2) * CD_vec;

% Setup the figure and axes, and plot the n = 0 reconstruction results
fig2 = figure('Name', 'Reconstructions', ...
              'Position', [400 140 1400 700]);

ax1 = subplot(2, 5, [1, 8]);

hold on;
x_ext = 2 * [-1 1];
y_ext = 2 * [-1 1];
z_ext = 2 * [-1 1];
line(ax1, x_ext, [0,0], [0,0], 'LineWidth', 1, 'Color', 'r');
line(ax1, [0,0], y_ext, [0,0], 'LineWidth', 1, 'Color', 'g');
line(ax1, [0,0], [0,0], z_ext, 'LineWidth', 1, 'Color', 'b');

%% ==shift postion offset back for plotting
x_p1=x_p1+offsetposition(1,1);
y_p1=y_p1+offsetposition(1,2);
z_p1=z_p1+offsetposition(1,3);

x_p2=x_p2+offsetposition(1,1);
y_p2=y_p2+offsetposition(1,2);
z_p2=z_p2+offsetposition(1,3);
%%  =============================

surf1 = surf(ax1, x_p1, y_p1, z_p1, p1plot, 'FaceColor', 'interp');
surf1.EdgeColor = 'none';
surf2 = surf(ax1, x_p2, y_p2, z_p2, p2plot, 'FaceColor', 'interp');
surf2.EdgeColor = 'none';

    % Plot the measurement points
    p_meas(:, m) = sim_meas_cart(x_meas, y_meas, z_meas, omega, pt_srcs);
%     scatter3(ax1, x_meas, y_meas, z_meas, ...
%              abs(p_meas(:, m)), abs(p_meas(:, m)));

    scatter3(ax1, x_meas, y_meas, z_meas,10,[0.8 0.8 0.8],'filled');

% If 'caxlims' is defined in 'params.mat', set colormap limits to
%  values of 'caxlims' in 'params.mat'.
% If 'cmap' is defined in 'params.mat', set colormap type to
%  value of 'cmap' in 'params.mat'.
if exist('params.mat', 'file')
    var_info = who('-file', 'params.mat');
    if ismember('caxlims', var_info)
        load('params.mat', 'caxlims');
        caxis(ax1, caxlims);
    end
    if ismember('cmap', var_info)
        load('params.mat', 'cmap');
        colormap(ax1, cmap);
    end
end
cbar = colorbar(ax1, 'eastoutside');

xlabel(ax1, 'X');        % label x axis
ylabel(ax1, 'Y');        % label y axis
zlabel(ax1, 'Z');        % label z axis
ax1.Title.String = sprintf('N = %d;  f = %.0f Hz\n', n, freqs(f_bin));
set(ax1, 'DataAspectRatio', [1 1 1]);
grid(ax1, 'on');
view(ax1, [130 20]);

ax2 = subplot(2, 5, [4, 10]);



hold off;

% Repeat reconstruction from n = 1 to Nmax and update plots
for n = 1:Nmax
    pause(0.5);
    tic;
    [CD_vec, res] = lstsq_solve(PSI_mat(:, 1:2*(n+1)^2), p_meas(:, f_bin));
    fprintf('N = %d, ', n)
    fprintf('f = %.0f Hz, compute time = %.3f; \n', freqs(f_bin), toc)
    residuals(n+1) = res;
    
    p1 = PSI_p1(:, 1:2*(n+1)^2) * CD_vec;
    p1 = reshape(p1, p1_shape);
    p1plot =20*log10(abs(p1)); %plot in dB values

    p2 = PSI_p2(:, 1:2*(n+1)^2) * CD_vec;
    p2 = reshape(p2, p2_shape);
    p2plot =20*log10(abs(p2)); %plot in dB values
    
    % Reconstruct at CTA2034 grid and compute convergence metrics
    p2034(:, n+1) = PSI_2034(:, 1:2*(n+1)^2) * CD_vec;
    conv(n) = norm(p2034(:, n) - p2034(:, n+1)) / p2034_len;
    
    % Update surface plots
    ax1.Title.String = sprintf('N = %d;  f = %.0f Hz\n', n, freqs(f_bin));
    surf1.CData = p1plot;
    surf2.CData = p2plot;
    
    % Plot least squares residuals and convergence metrics
    cla(ax2);
    semilogy(ax2, 1:n+1, residuals(1:n+1), 'bs-');
    hold(ax2, 'on');
    semilogy(ax2, 1:n, conv(1:n), 'rs-');
    ax2.Title.String = 'Fitting Residuals and CTA Convergence Metrics';
    legend(ax2, 'Residuals', 'CTA Convergence');
    xlabel(ax2, 'Expansion Order N');
    ylabel(ax2, 'Residuals/Convergence');
    grid(ax2, 'on');
    xlim(ax2, [0 Nmax]);
    ylim(ax2, [1e-10 1e2]);
    pbaspect(ax2, [1, 1, 1]);
    drawnow;
end

%% Error plot======================================================================================
% Allocate storage to store the measurements
p_meas = zeros(length(x_meas), length(freqs));

fig1 = figure('Name', 'Error', ...
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
%     p_d1 = sim_meas_cart(x_d1, y_d1, z_d1, omega, pt_srcs);
%     p_d1plot =20*log10(abs(p_d1)); %plot in dB values
    
    surf(ax1, x_d1, y_d1, z_d1, (p_d1plot-p1plot), ...
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
     caxis([-6 6]);
    colorbar(ax1, 'eastoutside');
    view(ax1, [130 20]);
    
    % Plot the XY disc grid
%     p_d2 = sim_meas_cart(x_d2, y_d2, z_d2, omega, pt_srcs);
%     p_d2plot =20*log10(abs(p_d2)); %plot in dB values

    surf(ax1, x_d2, y_d2, z_d2, p_d2plot-p2plot, ...
         'FaceColor', 'interp', 'EdgeColor', 'none');
    
    % Plot the measurement points
    p_meas(:, m) = sim_meas_cart(x_meas, y_meas, z_meas, omega, pt_srcs);
%     scatter3(ax1, x_meas, y_meas, z_meas, ...
%              abs(p_meas(:, m)), abs(p_meas(:, m)));

    scatter3(ax1, x_meas, y_meas, z_meas,10,'white','filled');
          
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
%======================================================================================