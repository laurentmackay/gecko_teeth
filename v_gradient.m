function grad = v_gradient(gamma, N)

x=(0:(N-1))'/(N-1);
grad=1+gamma*2*abs(x-0.5);

end