function [oind, tpeak,phi_l, phi_r, phi_ls, phi_rs] =get_phases(peaks)

N=length(peaks);



phi_rs = cell(N,1);
phi_rs{1} = cell([]);
phi_rs{end} = cell([]);

phi_ls = cell(N,1);
phi_ls{1} = cell([]);
phi_ls{end} = cell([]);


tpeak = {};
oind={};

phi_l={};
phi_r={};
for i=1:N

    if i>1
        left = peaks{i-1};
    else
        left=[];
    end

    if i<N
        right =  peaks{i+1};
    else
        right=[];
    end
    
    for j = 1:length(peaks{i})

        curr = peaks{i}(j);

        i_prev_right = find(right<=curr,1,'last');
        i_prev_left = find(left<=curr,1,'last');

        i_next_right = find(right>curr,1);
        i_next_left = find(left>curr,1);

        if ~isempty(i_prev_right) && ~isempty(i_next_right)
            c=curr - right(i_prev_right);
            d = right(i_next_right)-right(i_prev_right);
            r = c/d;

        else
            r=nan;
        end

        if ~isempty(i_prev_left)  && ~isempty(i_next_left)
            a= curr - left(i_prev_left);
            b = left(i_next_left)-left(i_prev_left);
            l = a/b;
        else
            l=nan;
        end

        phi_r{end+1} = r;
        phi_rs{i}{end+1} = r;

        phi_l{end+1} = l;
        phi_ls{i}{end+1} = l;



        tpeak{end+1}=curr;
        oind{end+1} = i;

    end
end

oind= cell2mat(oind);
tpeak=cell2mat(tpeak);
phi_r = cell2mat(phi_r);
phi_l = cell2mat(phi_l);


end