function [vals]=prior_valsf_beta(v_f,prior)
%this function calculates the negative log of the betapdf numerator
%denominator of the beta pdf cancels in calculation of acceptance ratio so we have omitted it here to streamline calculation

vals=sum(((1-prior.fid_alpha)*log(v_f)) + ((1-prior.fid_beta)*log(1-v_f)));
end