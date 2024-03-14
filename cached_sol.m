function [outputArg1,outputArg2] = cached_sol(ode_rhs_and_gamma, D, p ,L, tfinal)

ws=functions(ode_rhs_and_gamma).workspace{1}
c = struct2cell(ws)

str="";
for i=c
if isnumeric(i)
    str=[str strjoin(string(i))]
end
end

end