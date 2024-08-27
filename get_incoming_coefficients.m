function out = get_incoming_coefficients(CD)

  n = size(CD, 1);

  out = CD(n/2 + 1:end,:);

endfunction
