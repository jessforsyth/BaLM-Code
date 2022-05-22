function [saving]=file_generation(saving)

    if ~exist(saving.directory,'dir')
       mkdir(saving.directory)
    end
    
    saving.files={'theta_array.bin','av_acc_t.bin','beta_t.bin','d_array.bin', ...
    'av_acc_f.bin','beta_f.bin','f_prior.bin','permutation_array.bin',...
    'av_acc_P.bin','A_matrix.bin','temp_av_acc_t.bin', 'temp_beta_t.bin',...
    'temp_av_acc_P.bin','temp_av_acc_f.bin','temp_beta_f.bin','temperature.bin', 'temp_fidelity_terms.bin', ...
    'px.bin','py.bin','pz.bin',...
    'cov.bin','covariance_log.txt','cell2matchdist.bin'}; %should tidy this up, as you get empty files if you aren't sampling on all parameters but these can be deleted at some point!

    for file=1:length(saving.files)
        filename=char(strcat(saving.directory,saving.files(file)));
        f=fopen(filename,'w');
        fclose(f);
    end

end