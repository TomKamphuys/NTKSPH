function sph_viewer(CD)

  amount = 200;
  y_size = 0.3;
  Nmax = 12; % TODO from CD (size)
  temp = 273.15 + 20;


  theta = linspace(-pi, pi, amount);
  z = linspace(-y_size, y_size, amount);

  [THETA, Z] = meshgrid(theta, z);

  [X, Y, Z] = cyl2cart(0.3, THETA, Z);

  x_recon = X(:)';
  y_recon = Y(:)';
  z_recon = Z(:)';

    % Convert the reconstruction points coordinates from Cartesian to polar
  [phi_recon, theta_recon, r_recon] = cart2sph(x_recon, y_recon, z_recon);
  theta_recon = pi/2 - theta_recon;

  index = 1;
  CD_vec = CD(:,index);

  while (true)

    freqs = index*10;
    omega = 2*pi*freqs;

    figure(1)
    PSI_recon = sph_PSI_mix(r_recon, theta_recon, phi_recon, omega, Nmax, temp);
    p_recon = PSI_recon(:, 1:2:end) * CD_vec;

    out_recon = dB_SPL(reshape(p_recon, size(X)));
    pcolor(X, Y, out_recon)
    shading flat
    title('Outgoing only' )
    %axis equal
    %caxis([-50 -20]);
    colormap jet
    colorbar

    a = input('input');
    if a == 1
      CD_vec = CD(:,index+1);
    else
      CD_vec = CD(:,index-1);
    end
  endwhile

endfunction
