load data.mat

filename = 'soundpressure-high.gif';
DelayTime = 0.1;

z(r < 0.26) = [];
theta(r < 0.26) = [];

t = theta;
t(t>pi) = t(t>pi) - 2*pi;

f = figure(1)

for ind = 1:50:2000
  p_meas = p(:,ind);
  p_meas(r < 0.26) = [];

  scatter(t, z, [], dB_SPL(p_meas), 'filled');
  title(['index: ' num2str(ind)])
  caxis([-40, -5])
  drawnow

  % Image Processing
  % Assign plot to a frame
  frame = getframe(f);
  % Convert frame to RGB image (3 dimensional)
  im = frame2im(frame);
  % Transform RGB samples to 1 dimension with a color map "cm".
  [imind,cm] = rgb2ind(im);
  if i == 1;
      % Create GIF file
      imwrite(imind,cm,filename,'gif','DelayTime', DelayTime , 'Compression' , 'lzw');
  else
      % Add each new plot to GIF
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', DelayTime , 'Compression' , 'lzw');
  end


endfor

