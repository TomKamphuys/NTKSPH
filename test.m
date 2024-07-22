%[p, r, theta, z] = read_nfs_measurements();

freqs = 200;
Nmax = 5;

load 2ddata.mat;

z = z-0.3;
r = r + 0.03;


[x, y, z] = cyl2cart(r, theta, z);

% Temperature is 20 deg C.
temp = 273.15 + 20;



omega = 2*pi*freqs;

% Simulate the measurements
p_meas = p(:,round(freqs/10));
x_meas = x';
y_meas = y';
z_meas = z';

figure
%subplot(2, 3, 1)
%scatter(theta, z, [], dB_SPL(p_meas'), 'filled')
%title('Raw measurments')
%colormap jet

%figure
%scatter(x, y, [], dB_SPL(p_meas'), 'filled')
%colormap jet

% Plot the measurements and the locations of the points sources
plotting = false;
if plotting
    fig1 = figure('Name', 'Simulated Measurements and Reconstructions without SFS', ...
                  'Position', [100 190 1600 700]);

    scatter3(x_meas, y_meas, z_meas, 10, dB_SPL(p_meas), 'filled');
    axis('equal');
    view([130 20]);
%    colormap('hot');
%    caxis([85, 95]);
    colorbar('eastoutside');
    colormap jet
    title('Simulated Measurements and Reconstructions without SFS');
end

x_meas = x_meas';
y_meas = y_meas';
z_meas = z_meas';

[phi_meas, theta_meas, r_meas] = cart2sph(x_meas, y_meas, z_meas);
theta_meas = pi/2 - theta_meas;

for ind = 1:Nmax
  disp(ind)
  [PSI_mat, Nmax] = sph_PSI_mix(r_meas, theta_meas, phi_meas, omega, ind, temp);

  [CD_vec, ~] = lstsq_solve(PSI_mat, p_meas);

  [error_db(ind), error_percentage(ind)] = calc_error(r_meas, theta_meas, phi_meas, p_meas, CD_vec, omega, Nmax, temp);

  if error_percentage < 1
    break;
  end
endfor
disp(error_db)
disp(error_percentage)
[~, index] = min(error_percentage);
[PSI_mat, Nmax] = sph_PSI_mix(r_meas, theta_meas, phi_meas, omega, ind, temp);
[CD_vec, ~] = lstsq_solve(PSI_mat, p_meas);



[X, Y, outRef, outRecon] = cylwall_reconstruct(.3, 50, CD_vec, omega, Nmax, temp);

%figure
subplot(2, 3, 2)
pcolor(X, Y, outRef)
shading flat
title('Fitted Measurement' )
%axis equal
%caxis([-50 -20]);
colormap jet

colorbar

%figure
subplot(2, 3, 3)

pcolor(X, Y, outRecon)
shading flat
title('Outgoing only' )
%axis equal
%caxis([-50 -20]);
colormap jet
colorbar

%figure
subplot(2, 3, 5)
bar(abs(CD_vec(1:2:end)))
title('Outgoing coefficients')

%figure
subplot(2, 3, 6)
bar(abs(CD_vec(2:2:end)))
title('Incoming coefficients')

