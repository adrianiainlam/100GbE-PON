function r = rxFilter(y, rolloff, span, sps)
  %% Receiver matched (root raised cosine) filter,
  %% downsampling to 2 samples/sym.
  %% Inputs:
  %%  - y: received waveform
  %%  - rolloff: rolloff factor in root raised cosine filter
  %%  - span: filter span (number of symbols)
  %%  - sps: Input samples per symbol
  %% Output:
  %%  - r: filtered signal

  %% Construct filter object
  rxfilter = comm.RaisedCosineReceiveFilter...
                 ('Shape', 'Square root', ...
                  'RolloffFactor', rolloff, ...
                  'FilterSpanInSymbols', span, ...
                  'InputSamplesPerSymbol', sps, ...
                  'Gain', 1 / sqrt(sps));

  %% Perform filtering in frequency domain
  coef = coeffs(rxfilter);
  filter_fft = fft(coef.Numerator, length(y));
  y_fft = fft(y);
  rs = ifft(y_fft .* filter_fft.');

  %% Re-order signal due to circular convolution
  l = (length(coef.Numerator) - 1) / 2;
  rr = [rs(l:end); rs(1:l-1)];

  %% Downsample
  r = downsample(rr, sps/2, 2);
end
