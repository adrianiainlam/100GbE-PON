function demodData = deqpskdemod(modData)
  %% Input: Rx samples from a DE-QPSK constellation.
  %% Output: Data as integers between 0 to 3 inclusive.
  %% Parameter formats are the same as MATLAB's "pskdemod" function.

  %% MATLAB has no easy way to perform just the decision stage,
  % so here we demodulate and then remodulate to get the closest symbols.
  demod1 = pskdemod(modData, 4, 0, 'gray');
  remod = pskmod(demod1, 4, 0, 'gray');

  clear demod1; % save some memory

  [~, wavelength_channels] = size(modData);
  %% Instead of looping through the symbols, it will be faster
  % to vectorize the operation. So we create an array that is
  % delayed by 1 sample.
  delayed = [ones(1, wavelength_channels); remod(1:end-1, :)];
  demodData = uint8(pskdemod(remod .* conj(delayed), 4, 0, 'gray'));
end
