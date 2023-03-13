function time_effect(x0, beam, dsigma, lambda, arg)
    t0 = 0;
    t1 = 80;
    step = [0.1, 0.2, 0.5, 1, 2];
    misfit = zeros(1, length(step));
    for i=1:length(step)
        tspan = t0:step(i):t1;
        epsilon_exp = forana(beam, tspan', dsigma);
        misfit(i) =  misfit_sig(x0, epsilon_exp, tspan', beam, dsigma, lambda, arg);

    end
    figure
    plot(step, misfit, 'o-')
    xlabel('step size')
    ylabel('misfit')
    title('misfit vs step size')
    saveas(gcf, ['assets/misfit_vs_step_size', arg, '.png'])
end
