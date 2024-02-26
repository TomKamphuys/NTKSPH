function [f, p, angles] = read_speaker_measurement()

  angles = (0:10:180)'; % should depend on data in directory
  decimation_amount = 100;

  for ind = 1:length(angles)

%    waitbar(ind/length(angles));

    fileContent = importdata(sprintf('./Data/hor %d.txt', angles(ind)), ' ', 14);
    data = fileContent.data;

    % * Freq(Hz) SPL(dB) Phase(degrees)

    f = data(:,1);
    pt = to_pressure(data(:,2)) .* exp(i.*deg2rad(data(:,3)));
    f = f(1:decimation_amount:end);
    pt = pt(1:decimation_amount:end);

    p(ind,:) = pt';

  endfor

  p = [flipud(p); p(2:end,:)];
  angles = deg2rad([flipud(-angles); angles(2:end)]);


endfunction



