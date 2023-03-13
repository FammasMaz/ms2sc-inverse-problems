function eps_exp = low_pass_filter(tn, eps_exp, lambda)
        u = (1/lambda)*ones(lambda,1);
        eps_fil = conv(eps_exp, u, 'same');
        eps_exp = eps_fil;
end