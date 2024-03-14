function H = ClosureHash(fn)
% nm = functions(fn)
Engine = java.security.MessageDigest.getInstance('MD5');
H = CoreHash(fn, Engine, fn);
H = sprintf('%.2x', H);   % To hex string
function H = CoreHash(Data, Engine, func)
% Consider the type of empty arrays:
S = [class(Data), sprintf('%d ', size(Data))];
Engine.update(typecast(uint16(S(:)), 'uint8'));
H = double(typecast(Engine.digest, 'uint8'));
if isa(Data, 'struct')
   n = numel(Data);
   if n == 1  % Scalar struct:
      F = sort(fieldnames(Data));  % ignore order of fields
      for iField = 1:length(F)
          field = Data.(F{iField});
          if ~isa(field, 'function_handle') || ~isequal(field, func) 
            H = bitxor(H, CoreHash(Data.(F{iField}), Engine, func));
          end
      end
   else  % Struct array:
      for iS = 1:n
         H = bitxor(H, CoreHash(Data(iS), Engine, func));
      end
   end
elseif isempty(Data)
   % No further actions needed
elseif isnumeric(Data)
   Engine.update(typecast(Data(:), 'uint8'));
   H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
elseif ischar(Data) || isstring(Data) % Silly TYPECAST cannot handle CHAR
    if isstring(Data)
    Data = char(Data);
    end
   Engine.update(typecast(uint16(Data(:)), 'uint8'));
   H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
elseif iscell(Data)
   for iS = 1:numel(Data)
      H = bitxor(H, CoreHash(Data{iS}, Engine, func));
   end
elseif islogical(Data)
   Engine.update(typecast(uint8(Data(:)), 'uint8'));
   H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
elseif isa(Data, 'function_handle')

    H = bitxor(H, CoreHash(functions(Data), Engine, Data));
else
   warning(['Type of variable not considered: ', class(Data)]);
end