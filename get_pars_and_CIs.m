function varargout = get_pars_and_CIs(fits)
nms=coeffnames(fits{1});

for i =1:length(nms)
    nm=nms{i};
    varargout{2*i-1} = cellfun(@(cfit) cfit.(nm),fits);
    varargout{2*i} = cell2mat(cellfun(@(cfit) get_col(confint(cfit),i),fits ,'UniformOutput',false));
end
end