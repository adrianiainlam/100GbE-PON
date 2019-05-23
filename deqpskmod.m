function modData = deqpskmod(data)
  %% Input: Unmodulated data as integers between 0 to 3 inclusive.
  %% Output: Differentially encoded QPSK symbols.
  %% Parameter formats are the same as MATLAB's "pskmod" function.
  modData = pskmod(data, 4, 0, 'gray');
  numSymbs = size(data, 1);
  for i = 2 : numSymbs
    modData(i, :) = modData(i, :) .* modData(i - 1, :);
  end
end
