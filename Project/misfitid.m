function misfitid = misfitid(beam, tn, epsiexp, dsigma)
    epsilon = forana(beam, tn, dsigma);
    misfitid = sum((epsilon - epsiexp).^2);
end