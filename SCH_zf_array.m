function h = FHN_array(n,D,b)
a = 0.1;
e = 1e-2;
d2v=zeros(n,1);
    function duvdt = inner(~,uv)

        u = uv(1:n);
        v = uv((n+1):end);

        
        d2v(1) = (v(1)-v(2));
        for i = 2:n-1
            d2v(i) = (2.0*v(i)-v(i+1)-v(i-1));
        end
        d2v(n) = (v(n)-v(n-1));

        rx = u.*u.*v;
        duvdt = [e*(a-u+rx);
            e*(b-rx) - D*d2v];
    end

h=@inner;
end