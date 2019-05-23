function x = txFilter(modData, rolloff, span, sps)
  %% Transmitter pulse-shaping (root raised cosine) filter.
  %% Inputs:
  %%  - modData: modulated data
  %%  - rolloff: rolloff factor in root raised cosine filter.
  %%  - span: filter span (number of symbols)
  %%  - sps: samples per symbol
  %% Output:
  %%  - x: pulse-shaped waveform

  %% Construct filter object
  txfilter = comm.RaisedCosineTransmitFilter...
               ('Shape', 'Square root', ...
                'RolloffFactor', rolloff, ...
                'FilterSpanInSymbols', span, ...
                'OutputSamplesPerSymbol', sps, ...
                'Gain', sqrt(sps)); % so that output has energy 1

  %% Extract filter coefficients
  coef = coeffs(txfilter);
  %% Upsample data and perform filtering in frequency domain
  filter_fft = fft(coef.Numerator, length(modData) * sps);
  modData_fft = fft(upsample(modData, sps));
  x = ifft(modData_fft .* filter_fft.');
  %% Re-order signal due to circular convolution
  l = (length(coef.Numerator) - 1) / 2;
  x = [x(l:end); x(1:l-1)];
end
