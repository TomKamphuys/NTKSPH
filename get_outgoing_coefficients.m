function out = get_outgoing_coefficients(CD)

  n = size(CD, 1);

  out = CD(1:n/2,:);

endfunction
