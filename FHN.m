function dvwdt = FHN(~,vw)
        
         a = 0.7; 
         b = 0.8;
         g = 1/12.5;
         i = 0.5;
         v = vw(1);
         w = vw(2);
         dvwdt = [v-v^3/3 - w + i;
                  g*(v+a-b*w)];          
end