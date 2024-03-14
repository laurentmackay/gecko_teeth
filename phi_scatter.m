function h=phi_scatter(phi_r, phi_l, varargin)

h=scatter(phi_r, phi_l, varargin{:}); 

hold('on'); 
plot([0,1],[0,1],'-k'); 
hold('off');
xlabel('\phi_r','FontSize',14);
ylabel('\phi_l','FontSize',14,'Rotation',0);
end