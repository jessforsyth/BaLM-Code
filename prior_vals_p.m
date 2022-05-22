function [val]=prior_vals_p(v_p,data,prior)
%calculate the negative log of the prior density for momenta=v_p 
C=prior.mom_C^-1;

val=0.5*((v_p(:,1)'*C*v_p(:,1)) + (v_p(:,2)'*C*v_p(:,2)) + (v_p(:,3)'*C*v_p(:,3)));

end

