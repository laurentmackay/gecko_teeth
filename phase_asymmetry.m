function [chi, chis] =phase_asymmetry(phi_l, phi_r, phi_ls, phi_rs)

% chi_func = @(r,l) ((sin(pi*r).*sin(pi*l)).^4).*sin(pi*(r-l))/pi;
chi_func = @(r,l) sin(pi*(r-l))/pi;
chi = chi_func(phi_r,phi_l);
chis = cellfun(@(r,l) chi_func(cell2mat(r),cell2mat(l)), phi_rs, phi_ls, 'UniformOutput', false);

end