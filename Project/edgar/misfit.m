function val_misfit = misfit(dsigma,t_range,epsilon_data,varargin)
    for i = 1:3
        if length(varargin{i}) == 1
            index_r0 = i;
        end
    end
    index_range = setdiff(1:3,index_r0);
    range_1 = varargin{index_range(1)};
    range_2 = varargin{index_range(2)};

    parameters = zeros(1,3);
    val_misfit = zeros(length(range_1),length(range_2));
    
    for r1 = range_1
        i = find(range_1 == r1);
        for r2 = range_2
            j = find(range_2 == r2);
            parameters(index_range) = [r1 r2];
            parameters(index_r0) = varargin{index_r0};
            val_misfit(i,j) = integral(@(t) discrepancy(t,dsigma,t_range,epsilon_data,parameters),t_range(1),t_range(end),"Waypoints",t_range(2:end-1));
        end
    end
end