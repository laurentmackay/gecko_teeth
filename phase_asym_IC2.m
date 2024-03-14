function theta = phase_asym_IC(N,phi_l,chi)

theta=zeros(N,1);



for i=3:2:N

    theta(i) = theta(i-2)+chi*2*pi;
    disp(i)
end

for i=2:2:N-1

    theta(i) = mod((theta(i-1)+theta(i+1))/2,2*pi)+phi_l*2*pi;

end
% if mod(N,2)==1
%  theta(i) =theta(N-1)+phi_l*2*pi;
% end

theta = mod(theta, 2*pi);
end