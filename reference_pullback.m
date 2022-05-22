function [log_rp]=reference_pullback(tparam_prime,tparam)

%This function calculates the reference pullback term for bounded variables. 
%The correction is the same for [-pi,pi] and [0,1] bounded variables and so this function can be used for both cases. 
%tparam_prime = proposed parameter value
%tparam = current parameter value
%Here we calculate the log of the reference pullback (log_rp) as to avoid numerical errors. = log(ref_pullback) = log( (e^x/e^y) * ((e^y +1)/(e^x +1))^2)

%this calculation has already been simplified assuming param_prime and param are real values

log_rp= sum( tparam_prime - tparam + (2*log( (exp(tparam) +1)./(exp(tparam_prime) +1))) );

end
