function [val]=neg_log_PP(Y1,n1,fY2,F,Psi,nu,P)

X=Y1-fY2(:,P(1:n1)); %difference between Y1 and current transformed Y2 (re-ordered according to permutation vector P) 
val=log(det(Psi+(X*F*X')))*((nu+n1)/2); %this is the negative log of the posterior predictive, F= fidelity matrix. 

end

