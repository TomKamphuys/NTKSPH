%freqs = 20:10:800;
freqs = 810:100:20000;

for ind = 1:length(freqs)
  fit_error(ind) = plot_fit(r, theta, z, p, freqs(ind), 8);
endfor

figure
plot(freqs, fit_error);
xlabel('Frequency [Hz]')

ylabel('Fit error [dB]')
