function h = cartesian_phase_inhibitor_ode_nondim( N, gamma, eps, kappa, opts)
arguments
    N (1,1) {mustBeNumeric}
    gamma (1,1) {mustBeNumeric}
    eps (1,1) {mustBeNumeric}
    kappa (1,1) {mustBeNumeric}
    opts.T0 (1,1) {mustBeNumeric} = 1;
    opts.no_self (1,1) {mustBeNumericOrLogical} = false
end 
i=1:N;
locations = 2*((i-1)/(N-1) - 0.5);
T=opts.T0*(1+gamma*abs(locations'));
omega_0 = 1./T;
mu=1.0;
cij = coupling_strength_unbounded(i,i',kappa);
if opts.no_self
    cij(sub2ind([N,N],i,i))=0;
end
cij=(cij'./sum(cij))';
% cij=N*(cij/sum(cij(:)));
cij
% eps=eps/max(cij(:));

    function dydt = inner(t,sol)

        x=sol(1:N);
        y=sol(N+1:end);
        r=sqrt(x.^2+y.^2);
        theta=atan2(y,x);
        % disp(theta);

        theta_r=theta(3:end);
        theta_l=theta(1:end-2);
        theta_dot =  omega_0 - eps*(cij*secrete(theta)).*sense(theta);
        root = (mu-r.^2);
        % disp(root)
        % disp([secrete(theta(2)); (secrete(theta_r)+secrete(theta_l))/2; secrete(theta(end-1))].*sense(theta))
        % disp(theta_dot)
        dydt = [ root.*x - y.*theta_dot; root.*y + x.*theta_dot];
        % dydt = [  - y.*theta_dot;  x.*theta_dot];
    end


h=@inner;
end