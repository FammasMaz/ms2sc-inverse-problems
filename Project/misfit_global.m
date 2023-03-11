function [misfit_values, E_opt, K_opt, n_opt] = misfit_global(eps_exp, tn, beam, dsigma, E_range, K_range, n_range)
options = optimset('Display','off');
misfit_min = Inf;
misfit_values = zeros(length(E_range),length(K_range),length(n_range));
for i = 1:length(E_range)
    for j = 1:length(K_range)
        for k = 1:length(n_range)
            beamtemp = beam;
            beamtemp.E = E_range(i);
            beamtemp.K = K_range(j);
            beamtemp.n = n_range(k);
            [t, epsilon] = ode45(@forward_sigma, tn, 0, [], beamtemp, dsigma);
            %epsilon = forana(beamtemp, tn, dsigma);
            misfit_temp = quad(@discrepancy, tn(1), tn(end), [], 0, eps_exp, epsilon, n_range(k), tn);

            % theta = epsilon(:,1); % model prediction
            % thetaexp = eps_exp(:,1); % observation
            % RMSE = sqrt(mean((theta - thetaexp).^2));

            % misfit_temp = RMSE;

            misfit_values(i,j,k) = misfit_temp; % store misfit value
            if misfit_temp < misfit_min
                misfit_min = misfit_temp;
                E_opt = E_range(i);
                K_opt = K_range(j);
                n_opt = n_range(k);
            end
        end
    end
end

