function str = u_str(u)
data =  [ (1:length(u)); u ];
str = sprintf( [ '{' repmat('%d : %2.9E, ', 1, length(u)) '}'], data(:));
end