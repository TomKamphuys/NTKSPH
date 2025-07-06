function [CD, fit_error, CD_tot] = convert_to_coefficients(p, f, r, theta, z, N)

  temp = 273.15 + 20;

  beam_offset = 0.025 + 0.014; % meter
  arm_length = 0.79; % meter
  arm_offset = 0.013 + 0.01; % meter
  arm_angle = 1.0; % degrees
  theta = correct_for_setup(theta, r, arm_length, beam_offset, arm_offset, arm_angle); % setup is not perfect

  [x, y, z] = cyl2cart(r, theta, z);

  x_meas = x + 0.06; % - 0.043; % - 0.13; % acoustic center correction
  y_meas = y;
  z_meas = z; % - 0.09; % acoustic center correction

  [phi, theta, r] = cart2sph_phys(x_meas', y_meas', z_meas');

  sph_harm = calc_angular_part(phi, theta, N);

  h = waitbar(0);

  for ind = 1:length(f)

    waitbar(ind/length(f));

    kr = calc_kr(r, f(ind), temp);

    p_meas = p(:, ind);

    outgoing = calc_radial_part(kr, N) .* sph_harm;
    incoming = calc_radial_part_in(kr, N) .* sph_harm;

    total = [outgoing, incoming];

    CD_vec = total \ p_meas;

    fit_error(ind) = calc_error(p_meas, total, CD_vec);

    CD(:, ind) = get_outgoing_coefficients(CD_vec);
    CD_tot(:,ind) = CD_vec;
  endfor

  close(h)

endfunction



