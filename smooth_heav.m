function h = smooth_heav()
h = @(x, x0, eps) (1+tanh((x-x0)/eps))/2;
% h = @(x) (x>x0).*(1-exp((x-x0)/eps));
end