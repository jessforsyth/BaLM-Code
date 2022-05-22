function [data, general, beta, prior, samples, arrays,saving,tempering]=run_sampling_T1(filename,samp_def,samp_trans,samp_perm,samp_fid,data, general, beta, prior, samples, arrays,saving,tempering)

tempering.on=0;
tempering.T=1; %fix temperature at T==1
general.writing_count=0; 
general.C_transformation_temp=single(zeros(general.num_trans)); %re-start covariance storing for adaptive MCMC (easier just to wipe the vals here)
samples.avg_transformation=zeros(general.num_trans,1); %re-start average parameter value store
general.check_adj = 1.1;         %amount to scale the number of iterations between each adaptation of covariance
general.check_C=general.write_f; %when to first adapt proposal covariance, this value must be greater than fc 

%Sample at T=1
for it=1:general.N1
    
    if mod(it, general.N1/10)==0 
        %progress update in Chainx_log.txt
        f=fopen(filename,'a');
        fwrite(f,sprintf('Progress : %3.0f%% at %s \n',(it/general.N1)*100,datetime('now')));
        fclose(f);
    end

    if samp_trans==1
        %transformation sampling
        [data, general, beta, prior, samples, arrays,tempering]=transformation_sampling(it,samp_def,data, general, beta, prior, samples, arrays,saving,tempering);       
    end
    if samp_fid==1
        %fidelity sampling (data selection)
        [data, general, beta, prior, samples, arrays,tempering]=fidelity_sampling(it,data, general, beta, prior, samples, arrays,saving,tempering,samp_def);
    end
    if samp_perm==1
        %permutation sampling
        [data, general, beta, prior, samples, arrays,tempering]=rndpermutation_sampling(it,data, general, beta, prior, samples, arrays,saving,tempering,samp_def);
    end
    
    %whilst sampling at T==1 track progress of chain by calculating the average nlog of the posterior and keeping a store of the minimum nlog of the posterior density reached during sampling. 
    [samples, data]=min_check_nlog_posterior(it,samples, samp_def, samp_fid, data,prior,general);
    %calculate the cell to match distance for this current position within state space. 
    [samples,arrays]=calc_cell_to_match_dist(it,data,samples,saving,general,beta,tempering,samp_def,arrays);

end
%save matches matrix as text file once you have finished the sampling at T==1
f=fopen(strcat(saving.directory,char('matches_matrix.txt')),'w');
fprintf(f,saving.matches,arrays.matches);    
fclose(f);

end