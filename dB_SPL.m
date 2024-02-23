function dB_out = dB_SPL(x)
% Convert sound pressure 'x', which may be complex, into dB SPL. 
dB_out = 20 * log10(abs(x) / 20e-6);
end

