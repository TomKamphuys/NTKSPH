
clear data;

for ind = 1:3
  Nmax=2*ind;
  fit_and_reconstruct_2d_real_measurement;

  print(sprintf('frame%d.png', ind), 'png');

%  data(:,:,ind) = dB_SPL(p_recon');
end


im = imread ("animation.pdf", "Index", "all");
imwrite (im, "animation.gif", "DelayTime", .5)

%maxi = max(max(max(data)));
%mini = min(min(min(data)));
%data = uint8(((data - mini) / (maxi-mini))*2^8);
%
%imwrite(data(:,:,1),'animGif.gif','gif','writemode','overwrite',...
%        'LoopCount',inf,'DelayTime',0);
%
%%Loop through and write the rest of the frames
%for ii=2:size(data,3)
%     imwrite(data(:,:,ii),'animGif.gif','gif','writemode','append','DelayTime',0)
%end


