function grad = v_gradient(gamma, N)

x=(0:(N-1))'/(N-1);
grad=1+gamma*(x-0.5);

end