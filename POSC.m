function corr = POSC(recSig, hrtfSig)
[recR] = gccphat(recSig(:, 1), recSig(:, 2));
[hrtfR] = gccphat(hrtfSig(:, 1), hrtfSig(:, 2));
corr = sum(real(recR) .* real(hrtfR));
end


function r12 = gccphat(x, xref)
Ncorr = 2 * size(x, 1) - 1;
NFFT = 2 ^ nextpow2(Ncorr);
R12 = fft(x, NFFT) .* conj(fft(xref, NFFT)) ./ abs(fft(x, NFFT) .* conj(fft(xref, NFFT)));
r12_temp = fftshift(ifft(R12), 1);
r12 = r12_temp(NFFT / 2 + 1 - (Ncorr - 1) / 2 ...
    : NFFT / 2 + 1 + (Ncorr - 1) / 2, :);
end
