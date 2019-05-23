function [x, y] = ssfs(xin, yin, D, lambda, z, dz, Tsamp, gamma, alpha)
  %% Split-step Fourier solver, simulates chromatic dispersion,
  %% attenuation, Kerr effect, and power splitting / amplification.
  %% Params:
  %%  - xin, yin: input waveform (x and y polarizations)
  %%  - D: dispersion coefficient (ps / (nm km))
  %%  - lambda: wavelength (nm)
  %%  - z: length of fibre (km)
  %%  - dz: step size (km)
  %%  - Tsamp: sampling time (s)
  %%  - gamma: Non-linear coefficient (W^-1 / km)
  %%  - alpha: attenuation (dB / km)
  %% Output:
  %%  - x, y: output waveform (both polarizations)

  %% Convert everything to SI base units
  c = 299792458; % m/s
  D = D * 1e-6; % s/m^2
  lambda = lambda * 1e-9; % m
  z = z * 1e3; % m
  gamma = gamma * 1e-3; % watt^-1 / m
  dz = dz * 1e3; % m
  alpha = alpha / 10 * log(10); % Np/km
  alpha = alpha * 1e-3; % Np/m

  stepnum = z / dz;

  %% Frequency response of CD
  n = length(xin);
  fs = 1 / Tsamp;
  omega = (2*pi * fs / n * [(0 : floor((n-1)/2)), (-ceil((n-1)/2) : -1)]).';
  dispDFT = exp(-1j * omega.^2 * D * lambda^2 * dz / (4 * pi * c));

  %% Convenient variables to reduce typing
  hhz = 1j * (8/9) * gamma * dz; % Kerr phase shift per power
  attn = -alpha * dz / 2; % attenuation

  %% Initial Kerr half step
  P = abs(xin) .^ 2 + abs(yin) .^ 2;
  x = xin .* exp(P .* hhz / 2 + attn / 2);
  y = yin .* exp(P .* hhz / 2 + attn / 2);

  for i = 1 : stepnum
    %% CD in frequency domain
    xDFT = fft(x);
    yDFT = fft(y);
    x = ifft(xDFT .* dispDFT);
    y = ifft(yDFT .* dispDFT);

    %% Kerr effect in time domain
    P = abs(x) .^ 2 + abs(y) .^ 2;
    x = x .* exp(P .* hhz + attn);
    y = y .* exp(P .* hhz + attn);

    %% Power split after 40km
    if i * dz == 40e3
      splitnum = 1024; % energy factor, amplitude factor is sqrt of this
      x = x ./ sqrt(splitnum);
      y = y ./ sqrt(splitnum);
      %% Splitter loss - 1:4 coupler * 5 levels, 0.3 dB per level
      %% so total loss is 1.5 dB
      x = x ./ sqrt(10 ^ 0.15);
      y = y ./ sqrt(10 ^ 0.15);
    end
  end
  %% Final Kerr effect has overshot by half step, so cancel this
  P = abs(x) .^ 2 + abs(y) .^ 2;
  x = x .* exp(-P .* hhz / 2 - attn / 2);
  y = y .* exp(-P .* hhz / 2 - attn / 2);
end
