function h = FHN_array(n,D,mask)
         a = 0.7; 
         b = 0.8;
         g = 1/12.5;
         i = 0.5;
    function dvwdt = inner(~,vw)

         v = vw(1:n);
         w = vw((n+1):end);
         dv = [ -(v(end)-v(1)); diff(v) ];
         dw = [ -(w(end)-w(1)); diff(w) ];
         dvwdt = [v-v.^3/3 - w + i;
                  g*(v+a-b*w)]+D*[dv;dw];     
    end

    h=@inner;
end