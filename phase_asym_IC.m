function theta = phase_asym_IC(N,phi_l,delta)

theta=zeros(N,1);



for i=3:2:N

    theta(i) = theta(i-2)+delta*2*pi;

end

for i=2:2:N

    theta(i) = theta(i-1)+phi_l*2*pi;

end

theta = mod(theta, 2*pi);
end