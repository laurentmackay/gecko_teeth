function [psi, psis] =phase_checker(phi_l, phi_r, phi_ls, phi_rs)

% psi_func = @(r,l) (1-abs(cos(pi*(r+l)/2))).*(1-abs(sin(pi*(r-l)/2)));
% psi_func = @(r,l) 4*abs(sin(pi*(r-l)).*sqrt(sin(r).*sin(l)))/3;
psi_func = @(r,l) (sin(pi*r).*sin(pi*l)).^3;
psi = psi_func(phi_r,phi_l);
psis = cellfun(@(r,l) psi_func(cell2mat(r),cell2mat(l)), phi_rs, phi_ls, 'UniformOutput', false);

end