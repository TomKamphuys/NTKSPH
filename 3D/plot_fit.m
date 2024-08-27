function fit_error = plot_fit(r, theta, z, p, freq, Nmax)
  % plot_fit calculated the fit from the measurments (input arguments),
  % calculates the error of the fit and plots and saves the results in an image

r = r + 0.03;

[x, y, z] = cyl2cart(r, theta, z);

% Temperature is 20 deg C.
temp = 273.15 + 20;

freqs = freq;
omega = 2*pi*freqs;

% Simulate the measurements
p_meas = p(:,round(freqs/10));
x_meas = x';
y_meas = y';
z_meas = z';

figure(1)
subplot(2, 3, 1)
volume = dB_SPL(p_meas');
scatter(theta, z, [], volume - max(volume(:)), 'filled')
title('Raw measurments')
colormap jet
caxis([-25 0])


x_meas = x_meas' - 0.043;
y_meas = y_meas';
z_meas = z_meas';

[phi_meas, theta_meas, r_meas] = cart2sph_phys(x_meas, y_meas, z_meas);

sph_harm = calc_angular_part(phi_meas', theta_meas', Nmax);

kr = calc_kr(r, freq, temp);

outgoing = calc_radial_part(kr, Nmax) .* sph_harm;
incoming = calc_radial_part_in(kr, Nmax) .* sph_harm;

total = [outgoing, incoming];
[CD_vec, res] = lstsq_solve(total, p_meas);

fit_error = calc_error(p_meas, total, CD_vec);

[X, Y, outRef, outRecon] = cylwall_reconstruct(3, 100, CD_vec, omega, Nmax, temp);

subplot(2, 3, 2)
pcolor(X, Y, outRef - max(outRef(:)))
shading flat
title('Fitted Measurement' )
colormap jet
caxis([-25 0])

subplot(2, 3, 3)
pcolor(X, Y, outRecon - max(outRecon(:)))
shading flat
title('Outgoing only' )
colormap jet
caxis([-25 0])
%colorbar

subplot(2, 3, 4)
%text(0, 0, sprintf('Frequency: %d Hz; fit error: %4.2f', freq, fit_error));
title(sprintf('Frequency: %d Hz; fit error: %d dB', freq, fit_error))
axis off

subplot(2, 3, 5)
bar(abs(get_outgoing_coefficients(CD_vec)))
title('Outgoing coefficients')

subplot(2, 3, 6)
n = size(CD_vec, 1);
bar(abs(CD_vec(n/2+1:end)))
title('Incoming coefficients')

print(sprintf('%d.png', round(freq)), '-dpng');


endfunction
