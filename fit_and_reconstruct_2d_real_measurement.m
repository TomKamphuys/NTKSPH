close all
%clear all

clear p_recon;
clear bla;


% Temperature is 20 deg C.
temp = 273.15 + 20;

% size(vertical, horizontal)
% measurements size(m, 1) 1 is 'one'
% PSI_mat      size(m, c)
% coefs        size(c, 1) 1 is 'one

[f, p, angles] = read_speaker_measurement();
r = ones(size(angles));
Nmax = 15;

amount = 360;


freqs = f(1:1:end);

r_recon = 1.*ones(1,amount)';
a_recon = linspace(-pi, pi, amount)';
[x_recon, y_recon] = pol2cart(a_recon, r_recon);

for ind = 1:length(freqs)
  omega = 2*pi*freqs(ind);
  p_meas = p(:,ind);

  PSI_mat = sph_PSI_mix_2d(r, angles, omega, Nmax, temp);

  CD_vec = lstsq_solve(PSI_mat, p_meas);

  bla(ind,:) = abs(CD_vec);

  PSI_recon = sph_PSI_mix_2d(r_recon, a_recon, omega, Nmax, temp);

  p_recon(ind,:) = PSI_recon * CD_vec;

endfor

subplot(2, 2, 1)
pcolor(f, rad2deg(angles), dB_SPL(p))
shading flat
colormap(jet)
%set(gca,'xscale','log');
xlabel('frequency [Hz]')
ylabel('Angle [deg]')
colorbar
title('Measurement')
%ylim([-45, 45])


subplot(2, 2, 2)
pcolor(freqs, rad2deg(a_recon), dB_SPL(p_recon'))
shading flat
colormap(jet)
%set(gca,'xscale','log');
xlabel('frequency [Hz]')
ylabel('Angle [deg]')
colorbar
title('Reconstructed')
%ylim([-45, 45])

subplot(2, 2, 3)
pcolor(freqs, 1:(Nmax+1), (bla./repmat(max(bla, [], 2), 1, Nmax+1))');
shading flat
xlabel('frequency [Hz]')
ylabel('order [-]')
colorbar
title('Normalized Coefficients')

%subplot(2, 2, 4)
%pcolor(freqs, rad2deg(a_recon), dB_SPL(p) - dB_SPL(p_recon'))
%shading flat
%colormap(jet)
%%set(gca,'xscale','log');
%xlabel('frequency [Hz]')
%ylabel('Angle [deg]')
%colorbar
%title('Difference')
%ylim([-45, 45])



