for i = 1:8
 [x0,n] = FHN_sensing_vectors(i, f0=0); 
 disp(['X0(' num2str(i) ',:)=' sprintf('(/%0.5fd0, %0.5fd0/)',x0)]); 
 disp(['NORMAL(' num2str(i) ',:)=' sprintf('(/%0.5fd0, %0.5fd0/)',n) ])
disp(" ")
end