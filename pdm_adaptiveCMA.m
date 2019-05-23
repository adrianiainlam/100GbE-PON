function [x, y] = pdm_adaptiveCMA(rx, ry)
  %% Perform adaptive equalization using CMA.
  %% Input: rx, ry: Both polarizations of received signal
  %% Output: x, y: Equalizaed signal

  taps = 19; % Number of taps. Should be odd.
  mu = 1e-3; % Convergence parameter for gradient descent.

  hxx = zeros(taps, 1);
  hxy = zeros(taps, 1);
  hyx = zeros(taps, 1);
  hyy = zeros(taps, 1);
  %% hxx: real indices  -K, ..., 0, ..., K.   K = floor(taps/2)
  %%    MATLAB indices   1      1+K     taps

  %% Initialize hxx, hxx[0] = 1, hxx[k] = hxx[-k] = 0
  hxx(ceil(taps/2)) = 1;
  hxy(ceil(taps/2)) = 1;
  hyx(ceil(taps/2)) = 1;
  hyy(ceil(taps/2)) = 1;

  numSymbs = length(rx);

  %% Normalize to unit power.
  rx = rx / sqrt(mean(abs(rx) .^ 2));
  ry = ry / sqrt(mean(abs(ry) .^ 2));

  x = zeros(numSymbs, 1);
  y = zeros(numSymbs, 1);

  %% Run CMA twice so that the first symbols were also equalized
  for loops = 1:2
    %% Loop through each symbol
    for it = 1:numSymbs
      %% Construct block of length equal to filter length (taps)
      if it <= (taps - 1) / 2;
        %% If near the start, prepend zeros
        xp = [zeros((taps - 1) / 2 - it + 1, 1); rx(1:it + (taps - 1) / 2)];
        yp = [zeros((taps - 1) / 2 - it + 1, 1); ry(1:it + (taps - 1) / 2)];
      elseif it + (taps - 1) / 2 > numSymbs
        %% If near the end, append zeros
        xp = [rx(it - (taps - 1) / 2 : end); zeros(it + (taps - 1) / 2 - numSymbs, 1)];
        yp = [ry(it - (taps - 1) / 2 : end); zeros(it + (taps - 1) / 2 - numSymbs, 1)];
      else
        %% Just slice the signal
        xp = rx(it - (taps - 1) / 2 : it + (taps - 1) / 2);
        yp = ry(it - (taps - 1) / 2 : it + (taps - 1) / 2);
      end

      %% Filtering
      xout = sum(hxx .* xp) + sum(hxy .* yp);
      yout = sum(hyx .* xp) + sum(hyy .* yp);
      x(it) = xout;
      y(it) = yout;

      %% Caculate error
      ex = 1 - abs(xout) ^ 2;
      ey = 1 - abs(yout) ^ 2;

      %% Update filter by gradient descent
      hxx = hxx + mu * ex * xout * conj(xp);
      hxy = hxy + mu * ex * xout * conj(yp);
      hyx = hyx + mu * ey * yout * conj(xp);
      hyy = hyy + mu * ey * yout * conj(yp);

      %% If both filters converge to the same polarization,
      % re-initialize the filters.
      if sum(abs(hxx - hyx)) < 0.01 && sum(abs(hxy - hyy)) < 0.01
        hxx = 0.5 * (hxx + flipud(conj(hyy)));
        hyy = conj(flipud(hxx));
        hxy = 0.5 * (hxy - conj(flipud(hyx)));
        hyx = -conj(flipud(hxy));
      end
    end
  end
end
