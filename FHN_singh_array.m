function h = FHN_array(n,D,b)
a = 0.139;
k = 0.6;
e = 1e-3;
d2v=zeros(n,1);
    function duvdt = inner(~,uv)

        u = uv(1:n);
        v = uv((n+1):end);

        
        d2v(1) = (2.0*v(1)-v(2)-v(end));
        for i = 2:n-1
            d2v(i) = (2.0*v(i)-v(i+1)-v(i-1));
        end
        d2v(n) = (2.0*v(n)-v(1)-v(n-1));

        
        duvdt = [u.*(1-u).*(u-a)-v;
            e*(k*u-v-b) - D*d2v/2.0];
    end

h=@inner;
end