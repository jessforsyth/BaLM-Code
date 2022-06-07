function [val]=rel_posterior_withdef(x,general,data,prior,samples,samp_fid,MLM)

[A,b]=matrix_gen(x(1:12));
p=x(13:12+(3*data.n2));

if samp_fid==1
    F=diag(x(12+(3*data.n2)+1:end).^2);
else
    F=eye(data.n1);
end

P=MLM;

y0=[single(reshape(data.Y2',3*data.n2,1));single(p')];
[~,yout]=ode45(@(t,y) p_deform_mex(t,y,single(prior.kernel),single(data.n)),single([0,1]),y0);  %tspan=[0,1];
yfinal=yout(end,:)';
Y2d=double(reshape(yfinal(1:(3*data.n2)),data.n2,3)');

fY2=A*Y2d +b;

X=data.Y1-fY2(:,P(1:data.n1));
nlogpp=log(det(general.Psi+(X*F*X')))*((general.nu+data.n1)/2); %this is the negative log of the posterior predictive

p=reshape(p',[data.n2,3]);

C=diag(ones(1,data.n2)*prior.sigmom^2);
priors_momenta= 0.5*((p(:,1)'*C*p(:,1)) + (p(:,2)'*C*p(:,2)) + (p(:,3)'*C*p(:,3)));
priors_affine=0.5*(sum((x([4,5,6]).^2)./(prior.sigs^2)) + sum((x([10,11,12]).^2)./(prior.sigb^2)));
priors_fid= sum(((1-prior.fid_alpha)*log(x(12+(3*data.n2)+1:end))) + ((1-prior.fid_beta)*log(1-x(12+(3*data.n2)+1:end))));
if samp_fid==1
    fid_normconst=-general.d * sum(log(x(12+(3*data.n2)+1:end)));
else
    fid_normconst=0;
end

prior_terms=priors_affine + priors_fid + priors_momenta;
val=nlogpp + fid_normconst + prior_terms;

end