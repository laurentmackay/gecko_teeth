function c = coupling_strength(i,j,L,kappa,N)
s=cosh((L*(-1+i+j-N))/(N*kappa));
d=cosh((L-abs(L*(-i+j)/N))/kappa);
c=csch(L/kappa)*(s+d)/(2*kappa);
end