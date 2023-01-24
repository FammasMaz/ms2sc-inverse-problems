function X = rms(x)
X = sqrt(mean(x(isfinite(x)).^2));
