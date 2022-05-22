function [samples, data]=min_check_nlog_posterior(it,samples, samp_def, samp_fid, data,prior,general)
%this function keeps a track of the avg nlog of the posterior and the
%minimum value found over the course of sampling (i.e. the deepest point
%within the state space visited by the chain) 

if it==1 
    %i.e. if performing the first comparison
    if samp_def==0
        samples.prior_termd=0;
    end
    if samp_fid==0
        samples.prior_termf=0;
    end
end

nlog_prior_terms=samples.prior_terma + samples.prior_termd + samples.prior_termf; %calculate the nlog of the prior densities of current pos
samples.nlog_posterior=samples.n_log_PP + samples.nlog_normconst + nlog_prior_terms; %calculate the curren nlog of the posterior (prior and posterior pred)

if it==1
    samples.min_nlog_posterior=samples.nlog_posterior;
end

if samples.nlog_posterior<=samples.min_nlog_posterior 
    %comparison to identify minimum visited nlog posterior position
    
    samples.min_nlog_posterior=samples.nlog_posterior; %store corresponding vals
    samples.min_theta=samples.theta;
    samples.min_A=samples.A;
    samples.min_b=samples.theta(10:12);
    samples.min_P=samples.P;
    
    if samp_def==1
        samples.min_p=samples.p;
    else
        samples.min_p=0;
    end
    
    if samp_fid==1
        samples.min_fid=sqrt(diag(samples.f));
    else
        samples.min_fid=zeros(data.n1);
    end
    
end

%now calculate the average nlog posterior
if it>general.N1/2   %i.e. only caculate average of final half of N1 T=1 samples, i.e. region with no burn in! 
    it_half=it-general.N1/2;
    samples.avg_nlog_posterior= samples.avg_nlog_posterior + ((samples.nlog_posterior - samples.avg_nlog_posterior)/it_half);

    samples.avg_A=samples.avg_A + ((samples.A - samples.avg_A)./it_half);
    samples.avg_b=samples.avg_b + ((samples.theta(10:12)-samples.avg_b)./it_half);


    if samp_def==1
        samples.avg_p=samples.avg_p + ((samples.p-samples.avg_p)./it_half);
    else
        samples.avg_p=0;
    end

    if samp_fid==1
        samples.avg_fid=samples.avg_fid + ((sqrt(diag(samples.f))-samples.avg_fid)./it_half);
        samples.avg_nlog_normconst= fid_norm_const(general,samples.avg_fid); 
    else
        samples.avg_fid=zeros(data.n1);
    end
end


end
