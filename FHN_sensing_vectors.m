function [x0,n] = FHN_sensing_vectors(i, opts)
arguments 
    i (1,1) {mustBeInteger}
    opts.f0 (1,1) {mustBeNumeric} = 0

end
%all this stuff is highly parameter dependent, but it works well with
%k=0.6, a=0.169, and b in [0.1, 0.2] with e <= 1e-2
% 1->1
% 2->3
% 3->5
% 4->7
% 5->8
% 6->2
% 7->4
% 8->6

switch i
    case 1
        n=-[1,1.5]';
        x0=[0.75, 0];
    case 3
        n = [0.01,-1]';
        x0 = [0.0,0.11];
    case 5
        n = [1,2]';
        x0 = [0.0,0.04];

    case 7
        n = -[1,-30]';
        x0 = [0.0,0.005];

    case 8
     n = -[1,-4]';
     x0 = [0.8,0.04];
    case 2
     n = -[1,8]';
     x0 = [0.8,0.08];
    case 4
     n = [1,-6]';
     x0 = [0.0,0.08];
    case 6
    n = [1,9]';
    x0 = [0.0,0.02];


end
x0(2)=x0(2)+opts.f0;
n=n/norm(n);