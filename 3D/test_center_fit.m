function error_percentage = test_center_fit(corr, ind)

  x_corr = corr(1);
  z_corr = corr(2);
  r_corr = corr(3);

  load('22082024_measurement.mat');

  r = r + 0.03;

  beam_offset = 0.025; % meter
  arm_offset = 0.0023; % meter
  arm_angle = 1.46; % degrees
  theta = correct_for_setup(theta, r, beam_offset, arm_offset, arm_angle); % setup is not perfect

  N = 12;
  temp = 273.15 + 20;

  z = z;
  r = r - r_corr;

  [x, y, z] = cyl2cart(r, theta, z);

  x_meas = x - x_corr; % acoustic center correction
  y_meas = y;
  z_meas = z - z_corr; % acoustic center correction

  [phi, theta, r] = cart2sph_phys(x_meas', y_meas', z_meas');

  sph_harm = calc_angular_part(phi, theta, N);

%  ind = 1000;
%  f = f(ind);

  kr = calc_kr(r, f(ind), temp);

%  for ind = 1:length(f)
    p_meas = p(:, ind);

    outgoing = calc_radial_part(kr, N) .* sph_harm;
    incoming = calc_radial_part_in(kr, N) .* sph_harm;

    total = [outgoing, incoming];

    CD_vec = total \ p_meas;

%    CD(:, ind) = get_outgoing_coefficients(CD_vec);
%  endfor

  omega = 2*pi*f(ind);

  [error_dB, error_percentage] = calc_error(p_meas, total, CD_vec);

endfunction

