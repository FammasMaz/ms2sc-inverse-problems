function epsilon=forana(beam, tn, dsigma)
    epsilon = ((dsigma.*tn)/beam.E) + ((dsigma.*tn).^beam.n.*tn)/(beam.K.^beam.n.*(beam.n+1));
end
