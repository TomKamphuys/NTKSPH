function [hor, ver] = take_cea2034_measurements(CD, freqs)

  N = calcNmax(size(CD, 1));

  [phi, theta, bla] = cea2034_measurements();
  r = 1;

  [spl, phase] = take_virtual_measurement(CD, phi, theta, r, freqs);

  hor.angle = -180:10:180;
  ver.angle = -180:10:180;

  hor.f = repmat(freqs', [37, 1]);
  ver.f = hor.f;
  hor.spl = spl(1:37, :);
  ver.spl = spl(38:end, :);
  hor.phase = phase(1:37,:);
  ver.phase = phase(38:end,:);

endfunction
