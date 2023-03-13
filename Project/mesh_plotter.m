function mesh_plotter(misfit_values, a_range, b_range, c_range, z_axis, beam, arg)

    if z_axis == 'K'
        idx = find(b_range == beam.K);
        misfit_slice = squeeze(misfit_values(:,idx,:))';
        xlab = 'E';
        ylab = 'n';
        f = a_range;
        g = c_range;
    elseif z_axis == 'n'
        idx = find(c_range == beam.n);
        misfit_slice = squeeze(misfit_values(:,:,idx));
        xlab = 'K';
        ylab = 'e';
        f = b_range;
        g = a_range;
    elseif z_axis == 'E'
        idx = find(a_range == beam.E);
        misfit_slice = squeeze(misfit_values(idx,:,:))';
        xlab = 'K';
        ylab = 'n';
        f = b_range;
        g = c_range;
    else
        error('Invalid ranges');
    end    
    [A_mesh, b_mesh] = meshgrid(f, g);
    figure
    mesh(A_mesh, b_mesh, misfit_slice); % plot for n = b_range(1)
    xlabel(xlab);
    ylabel(ylab);
    zlabel('Misfit');
    title([ 'Misfit for ', xlab, ' and ', ylab]);
    saveas(gcf, strcat('assets/misfit_', xlab, '_', ylab, arg, '.png'));
end