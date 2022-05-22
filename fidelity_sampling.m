function [data, general, beta, prior, samples, arrays,tempering]=fidelity_sampling(it,data, general, beta, prior, samples, arrays,saving,tempering,samp_def)
%Sampling on the fidelity terms that make up the matrix D, which is n1xn1
%these terms describe how much we 'pay attention to cells' 
%terms close to 0 indicate that the cells are not weighted highly, i.e. for
%an untrue match

%First check the acceptance rate is within tolerance limits. 
if mod(it,general.fc)==0
    general.alpha_f(general.fc)=samples.alpha_f;
    general.alpha_bar_f=mean(general.alpha_f);
    if general.alpha_bar_f > (general.acc_opt + general.acc_tol)
        beta.f=beta.f*(1+beta.adj); 
    elseif general.alpha_bar_f < (general.acc_opt - general.acc_tol)
        beta.f=beta.f*(1-beta.adj);
    end
else
    general.alpha_f(mod(it,general.fc))=samples.alpha_f;
end

%Propose random walk on ftau (transformed params)
vftau=samples.ftau + (beta.f*normrnd(0,prior.fid_propstd,1,data.n1));
vf=1./(exp(vftau)+1); %transform back to fidelity params 
vf(vf==1)=1-1e-10;
vf(vf==0)=1e-10;
F=diag(vf).^2; %generate fidelity matrix

%we will assume for now that this module will ALWAYS be run after transformation sampling
n_log_PP=neg_log_PP(data.Y1,data.n1,samples.fY2,F,general.Psi,general.nu,samples.P); %calculate the nlog of the posterior predictive
newprior=prior_valsf_beta(vf,prior);            %nlog of prior density on fidelity parameters 
prior_terms=  samples.prior_termf - newprior;   %negative logs of priors
log_rp= reference_pullback(vftau,samples.ftau); %reference pullback to account for us sampling on the transformed parameter not the fidelities directly
nlog_normconst=fid_norm_const(general,vf);      %negative log of normalisation constant factor

%calculate acceptance ratio 
alpha=min(1,exp( ((1/tempering.T)*(samples.n_log_PP - n_log_PP + samples.nlog_normconst - nlog_normconst)) + prior_terms + log_rp) );

if alpha>rand()
    %accept
    samples.n_log_PP=n_log_PP;
    samples.f=F;
    samples.ftau=vftau;
    samples.prior_termf=newprior;
    samples.nlog_normconst=nlog_normconst;
    
end
samples.alpha_f=alpha;

%append values to arrays and save data
if mod(it,general.write_f)==0
    write_your_data('f',saving, general, beta, tempering,samp_def,arrays,samples,data)
elseif  mod(it,general.thin)==0 && mod(it,general.write_f)~=0
        itap=(it/general.thin)-(floor(it/general.write_f)*(general.write_f/general.thin));
        arrays.f(itap,:)=single(sqrt(diag(samples.f)));
        arrays.accr_f(itap)=single(general.alpha_bar_f);
        arrays.beta_f(itap)=single(beta.f);
end

end
   
