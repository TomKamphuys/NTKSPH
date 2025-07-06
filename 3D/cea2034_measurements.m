function [phi, theta, r] = cea2034_measurements()

  r = 1 * ones(74, 1); % m
  phi = [-180:10:180 180*ones(1,10) zeros(1,18) 180*ones(1,9)]'/180*pi;
  theta = [zeros(1, 37) 90:-10:0 10:10:180 170:-10:90]'/180*pi;

endfunction
