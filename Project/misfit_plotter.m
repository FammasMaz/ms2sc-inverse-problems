function misfit_plotter(a_range, b_range, c_range, epsilon_exp, t, beam, dsigma, z_axis)
    if z_axis == 'E'
        ran = a_range;
        value = zeros(size(ran));
        for i=1:length(ran)   
            x = [ran(i); beam.K; beam.n]; 
            value(i)=misfit_sig(x, epsilon_exp, t, beam, dsigma);
        end
    elseif z_axis == 'K'
        ran = b_range;
        value = zeros(size(ran));
        for i=1:length(ran)   
            x = [beam.E; ran(i); beam.n]; 
            value(i)=misfit_sig(x, epsilon_exp, t, beam, dsigma);
        end
    elseif z_axis == 'n'
        ran = c_range;
        value = zeros(size(ran));
        for i=1:length(ran)   
            x = [beam.E; beam.K; ran(i)]; 
            value(i)=misfit_sig(x, epsilon_exp, t, beam, dsigma);
        end
    end
    
    figure
    plot(ran, value);
    xlabel(z_axis);
    ylabel('Misfit');
    title(['Misfit Values for Different Values of ', z_axis]);
    saveas(gcf, ['assets/misfit_', z_axis, '.png']);
end