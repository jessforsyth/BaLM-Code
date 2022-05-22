function [data, general, beta, prior, samples, arrays,saving,tempering]=run_tempering(it,filename,samp_def,samp_trans,samp_perm,samp_fid,data, general, beta, prior, samples, arrays,saving,tempering) 
%Here we run each of the modules during tempered sampling. This continues
%until T<=1, then we exit this function and move onto sampling at T==1. 

if tempering.on==1  
    while tempering.T>1 && it<(general.N*5) %the second conditionis so that if you get stuck at a temperature for a long time, the code will exit and allow you to debug. e.g. really bad trapping!
        it=it+1;
        if mod(it, general.N/10)==0  
            %progress update written to Chainx_log.txt file
            f=fopen(filename,'a');
            fwrite(f,sprintf('Tempering Progress : %3.0f at %s \n',tempering.T,datetime('now')));
            fclose(f);
        end
        
        if samp_trans==1  
            %sample on affine and if samp_def==1 also non-linear deformation
            [data, general, beta, prior, samples, arrays,tempering]=transformation_sampling(it,samp_def,data, general, beta, prior, samples, arrays,saving,tempering);       
        end
        if samp_fid==1    
            %sample on fidelity term
            [data, general, beta, prior, samples, arrays,tempering]=fidelity_sampling(it,data, general, beta, prior, samples, arrays,saving,tempering,samp_def);
        end
        if samp_perm==1   
            %sample on permutation vector
            [data, general, beta, prior, samples, arrays,tempering]=rndpermutation_sampling(it,data, general, beta, prior, samples, arrays,saving,tempering,samp_def);
        end
    end

    %if T==1 or the previous while loop has some condition met, make sure
    %we finish sampling up till we perform a final write of the data
    %(otherwise we would loose the most recent samples!)
    if mod(it,general.write_f)~=0        
        for it=it+1:((floor(it/general.write_f)+1)*general.write_f)
            if mod(it, general.N/10)==0 
                f=fopen(filename,'a');
                fwrite(f,sprintf('Tempering Progress : %3.0f at %s \n',tempering.T,datetime('now')));
                fclose(f);
            end

            if samp_trans==1  %sample on affine and if samp_def==1 also deformation
                [data, general, beta, prior, samples, arrays,tempering]=transformation_sampling(it,samp_def,data, general, beta, prior, samples, arrays,saving,tempering);       
            end
            if samp_fid==1    %sample on fidelity term
                [data, general, beta, prior, samples, arrays,tempering]=fidelity_sampling(it,data, general, beta, prior, samples, arrays,saving,tempering,samp_def);
            end
            if samp_perm==1   %sample on permutation vector
                [data, general, beta, prior, samples, arrays,tempering]=rndpermutation_sampling(it,data, general, beta, prior, samples, arrays,saving,tempering,samp_def);
            end
        end
    end  
    general.tempN=it; %this tells us how many samples we actually performed during tempering. 
    
    %Plot figures
    f=fopen(filename,'a');
    fwrite(f,sprintf('Tempering Finished at %s \n',datetime('now')));
    fclose(f); 
    general.writing_count_temp=general.writing_count;%so we know how many times we wrote to the bin files
    figure_generation(tempering,data, general, samples, saving, samp_trans, samp_def, samp_fid, samp_perm, prior,arrays, filename) 
end  
f=fopen(filename,'a');
fwrite(f,sprintf('Sampling at T=1 at %s \n',datetime('now')));
fclose(f); 

end