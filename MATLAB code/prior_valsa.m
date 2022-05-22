function [val]=prior_valsa(v_theta,prior)
%this function calculates the -log(pdf) at each of the values within theta 

val=0.5*(sum((v_theta([4,5,6]).^2)./(prior.sigs^2)) + sum((v_theta([10,11,12]).^2)./(prior.sigb^2)));
end

