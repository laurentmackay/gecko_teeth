function c = coupling_strength_unbounded(i,j,kappa)
c=exp(-abs(i-j)/kappa)/(2*kappa);
end