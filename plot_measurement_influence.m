function plot_measurement_influence(r, theta, z, PSI_mat)
  [x, y, z] = cyl2cart(r, theta, z);


  A = PSI_mat;
  B = A'*A;
  C = B^-1;
  D = C*A';

  orders = size(D, 1);
  rows = ceil(sqrt(orders/2));

  minimum = 999999999999999999;
  maximum = -999999999999999999;

  figure
  for ind = 2:2:orders
%    figure(ind/2)
    subplot(rows, rows, ind/2)
%    scatter(theta, z, [], abs(D(ind,:)), 3, 'filled')
    influence = abs(D(ind,:));
    minimum = min(minimum, min(influence));
    maximum = max(maximum, max(influence));

    scatter3(x, y, z, 3, influence, 'filled');
    axis off;
    axis equal;
    colormap jet
    colorbar
  end

  for ind = 2:2:orders
%    figure(ind/2)
    subplot(rows, rows, ind/2)
    caxis([minimum, maximum]);
  end

endfunction
