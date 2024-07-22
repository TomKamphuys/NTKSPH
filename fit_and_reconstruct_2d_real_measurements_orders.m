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
Nmax = 25;

amount = 360;

f_index = 500;
freqs = f(f_index);

r_recon = 1.*ones(1,amount)';
a_recon = linspace(-pi, pi, amount)';
[x_recon, y_recon] = pol2cart(a_recon, r_recon);

bla = zeros(Nmax, Nmax);

for n = 1:Nmax
  omega = 2*pi*freqs;
  p_meas = p(:,f_index);

  PSI_mat = sph_PSI_mix_2d(r, angles, omega, n, temp);

  CD_vec = lstsq_solve(PSI_mat, p_meas);

  bla(n,1:n+1) = abs(CD_vec);

  PSI_recon = sph_PSI_mix_2d(r_recon, a_recon, omega, n, temp);

  p_recon(n,:) = PSI_recon * CD_vec;

endfor





