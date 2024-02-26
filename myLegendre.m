function out = myLegendre(n, x)

  out = 0;
  for k = 0:floor(n/2)
    out = out + (-1)^k .* factorial(2.*n-2.*k) ./ (2.^n .* factorial(k) .* factorial(n-k) .* factorial(n-2*k)) .* x.^(n-2.*k);
  end


endfunction
