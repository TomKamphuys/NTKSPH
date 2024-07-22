function fit_error = plot_fit(r, theta, z, p, freq, Nmax)


[x, y, z] = cyl2cart(r, theta, z);

% Temperature is 20 deg C.
temp = 273.15 + 20;

% To make things simpler, simulate 1 frequency only
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


x_meas = x_meas';
y_meas = y_meas';
z_meas = z_meas';

[phi_meas, theta_meas, r_meas] = cart2sph(x_meas, y_meas, z_meas);
theta_meas = pi/2 - theta_meas;

[PSI_mat, Nmax] = sph_PSI_mix(r_meas, theta_meas, phi_meas, omega, Nmax, temp);

[CD_vec, res] = lstsq_solve(PSI_mat, p_meas);

fit_error = calc_error(r_meas, theta_meas, phi_meas, p_meas, CD_vec, omega, Nmax, temp);


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
bar(abs(CD_vec(1:2:end)))
title('Outgoing coefficients')

subplot(2, 3, 6)
bar(abs(CD_vec(2:2:end)))
title('Incoming coefficients')

print(sprintf('%d.png', round(freq)), '-dpng');


endfunction
