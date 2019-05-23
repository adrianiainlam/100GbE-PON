function x = CDCompensation(xCD, D, lambda, z, Tsamp)
  %% Chromatic dispersion compensation.
  %% Params:
  %%  - xCD: received waveform with CD
  %%  - D: dispersion coefficient (ps / (nm km))
  %%  - lambda: wavelength (nm)
  %%  - z: length of fibre (km)
  %%  - Tsamp: sampling time (s)
  %% Output:
  %%  - x: xCD after performing CD compensation

  %% Convert everything to SI base units
  c = 299792458; % m/s
  D = D * 1e-6; % s/m^2
  lambda = lambda * 1e-9; % m
  z = z * 1e3; % m

  %% Discrete FIR filter:
  %% N: filter length;   k: filter index;   h: filter coefficients.
  N = 2 * floor(abs(D) * lambda^2 * z / (2 * c * Tsamp^2)) + 1;
  k = -floor(N / 2) : floor(N / 2);
  h = exp(-1j * pi * c * Tsamp^2 * k .^ 2 / (D * lambda^2 * z));

  %% Perform filtering in frequency domain
  len_fft = max(length(xCD), length(h));
  H = fft(h, len_fft);
  XCD = fft(xCD, len_fft);
  x = ifft(H.' .* XCD);

  %% Re-order due to circular convolution
  l = (N - 1) / 2;
  if l > 0
    x = [x(l:end); x(1:l-1)];
  else
    x = [x(end); x(1:end-1)];
  end
end
