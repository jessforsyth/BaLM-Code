function [val]=rel_posterior_nodef(x,general,data,prior,samples,samp_fid,MLM)

[A,b]=matrix_gen(x(1:12));

if samp_fid==1
    F=diag(x(13:end).^2);
else
    F=eye(data.n1);
end

P=MLM;%samples.min_P;

Y2d=data.Y2;
fY2=A*Y2d +b;

X=data.Y1-fY2(:,P(1:data.n1));
nlogpp=log(det(general.Psi+(X*F*X')))*((general.nu+data.n1)/2); %this is the negative log of the posterior predictive

priors_affine=0.5*(sum((x([4,5,6]).^2)./(prior.sigs^2)) + sum((x([10,11,12]).^2)./(prior.sigb^2)));
priors_fid= sum(((1-prior.fid_alpha)*log(x(13:end))) + ((1-prior.fid_beta)*log(1-x(13:end))));
if samp_fid==1
    fid_normconst=-general.d * sum(log(x(13:end)));
else
    fid_normconst=0;
end

prior_terms=priors_affine + priors_fid;
val=nlogpp + fid_normconst + prior_terms;

end