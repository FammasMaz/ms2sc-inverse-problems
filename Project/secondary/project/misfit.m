function val_misfit = misfit(x,dsigma,t_range,epsilon_data)
    val_misfit = integral(@(t) discrepancy(t,x,dsigma,t_range,epsilon_data),t_range(1),t_range(end),"Waypoints",t_range(2:end-1));
end