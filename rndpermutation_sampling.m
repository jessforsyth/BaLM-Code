function [data, general, beta, prior, samples, arrays,tempering]=rndpermutation_sampling(it,data, general, beta, prior, samples, arrays,saving,tempering,samp_def)
%this funcion woks to sample a single swap of two cells per iteration in order to sample from the permutation state space. 

%first check the acceptance rate (but we don't do anything about it if it
%is outside of preferred limits!)
if mod(it,general.fc)==0
    general.alpha_P(general.fc)=samples.alpha_P;
    general.alpha_bar_P=mean(general.alpha_P);
else
    general.alpha_P(mod(it,general.fc))=samples.alpha_P;
end

%propose two random points within P to swap
k=randi(data.n2-data.n_refpoints);
l=randi(data.n2-data.n_refpoints);
while k==l %don't allow the swapping of the same point
    l=randi(data.n2-data.n_refpoints); 
end
%perfom swap of terms k and l to generate permutation vector proposal Q
Q=samples.P;
Q(k)=samples.P(l);
Q(l)=samples.P(k);

n_log_PP=neg_log_PP(data.Y1,data.n1,samples.fY2,samples.f,general.Psi,general.nu,Q); %calculate the nlog of the posterior predictive
%calculate the acceptance ratio
alpha=min(1,exp((1/tempering.T)*(samples.n_log_PP - n_log_PP)));

%now figure out if we should accept or reject Q (the permutation vector proposal)
if alpha>rand()
    %accept
    samples.P=Q;
    samples.n_log_PP=n_log_PP;
end

%now save out any data
samples.alpha_P=alpha;

%append values to arrays
if mod(it,general.write_f)==0
    write_your_data('p',saving, general, beta, tempering,samp_def,arrays,samples,data)
elseif mod(it,general.thin)==0 && mod(it,general.write_f)~=0
    itap=(it/general.thin)-(floor(it/general.write_f)*(general.write_f/general.thin));
    arrays.accr_p(itap)=single(general.alpha_bar_P);
end

%now count the number of matches
if tempering.on==0 %only count if sampling at T==1
    for k=1:data.n1 
        arrays.matches(k,samples.P(k))=arrays.matches(k,samples.P(k))+1;
    end
end

end  