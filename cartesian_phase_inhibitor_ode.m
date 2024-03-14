function h = cartesian_phase_inhibitor_ode( mu, T, eps, m)


omega_0 = (2*pi)./T;

    function dydt = inner(t,sol)
        N=length(sol)/2;
        x=sol(1:N);
        y=sol(N+1:end);
        r=sqrt(x.^2+y.^2);
        theta=atan2(y,x);
        % disp(theta);

        theta_r=theta(3:end);
        theta_l=theta(1:end-2);
        % disp(sense(theta))
        theta_dot =  omega_0 - eps*[secrete(theta(2)); (secrete(theta_r)+secrete(theta_l))/2; secrete(theta(end-1))].*sense(theta);
        % theta_dot =  omega_0 - eps*[secrete(theta(2)); (secrete(theta_r)+secrete(theta_l)); secrete(theta(end-1))].*sense(theta);
        root = (mu-r.^2);
        % disp([secrete(theta(2)); (secrete(theta_r)+secrete(theta_l))/2; secrete(theta(end-1))].*sense(theta))
        % disp(theta_dot)
        dydt = [ root.*x - y.*theta_dot; root.*y + x.*theta_dot];
    end


h=@inner;
end