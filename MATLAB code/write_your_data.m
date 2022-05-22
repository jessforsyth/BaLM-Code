function write_your_data(cur_samp,saving, general, beta, tempering,samp_def,arrays,samples,data)
%Writing out data for affine, deformation, fidelity and permutation sampling.
%'a'=transformation, 'f'=fidelity, 'p'=permutation, 'd'=cell to match distances
%This function stops us having to save ALL samples in array during
%sampling, and stops us having to write out every single iteration. 

finalentry=(general.write_f/general.thin);
if cur_samp=='a' && tempering.on==1
    arrays.accr_t(finalentry)=single(general.alpha_bar_t);
    arrays.beta_t(finalentry)=single(beta.t);
    arrays.temp(finalentry)=single(tempering.T);

    fileIDA=fopen(strcat(saving.directory, 'temp_av_acc_t.bin'),'a');
    fileIDB=fopen(strcat(saving.directory, 'temp_beta_t.bin'),'a');
    fileIDC=fopen(strcat(saving.directory,'temperature.bin'),'a');

    fwrite(fileIDA,arrays.accr_t','single');
    fwrite(fileIDB,arrays.beta_t','single');
    fwrite(fileIDC,arrays.temp','single');

    fclose(fileIDA);
    fclose(fileIDB);
    fclose(fileIDC);
elseif cur_samp=='a' && tempering.on==0
    arrays.accr_t(finalentry)=single(general.alpha_bar_t);
    arrays.beta_t(finalentry)=single(beta.t);
    arrays.temp(finalentry)=single(tempering.T);
    arrays.theta(finalentry,:)=single(samples.theta);
    arrays.A(finalentry,:)=single(samples.A);

    fileID=fopen(strcat(saving.directory, 'theta_array.bin'),'a');
    fileID2=fopen(strcat(saving.directory, 'av_acc_t.bin'),'a');
    fileID3=fopen(strcat(saving.directory, 'beta_t.bin'),'a');
    fileID4=fopen(strcat(saving.directory,'A_matrix.bin'),'a');


    fwrite(fileID,arrays.theta','single');
    fwrite(fileID2,arrays.accr_t','single');
    fwrite(fileID3,arrays.beta_t','single');
    fwrite(fileID4,arrays.A','single');

    fclose(fileID);
    fclose(fileID2);
    fclose(fileID3);
    fclose(fileID4);

    if samp_def==1
        arrays.px(finalentry,:)=single(samples.p(1,:));
        arrays.py(finalentry,:)=single(samples.p(2,:));
        arrays.pz(finalentry,:)=single(samples.p(3,:));

        fileID5=fopen(strcat(saving.directory,'px.bin'),'a');
        fileID6=fopen(strcat(saving.directory,'py.bin'),'a');
        fileID7=fopen(strcat(saving.directory,'pz.bin'),'a');
        fwrite(fileID5,arrays.px','single');
        fwrite(fileID6,arrays.py','single');
        fwrite(fileID7,arrays.pz','single');
        fclose(fileID5);
        fclose(fileID6);
        fclose(fileID7);
    end
end
    
if cur_samp=='f' && tempering.on==1
    arrays.accr_f(finalentry)=single(general.alpha_bar_f);
    arrays.beta_f(finalentry)=single(beta.f);

    fileIDA=fopen(strcat(saving.directory,'temp_av_acc_f.bin'),'a');
    fileIDB=fopen(strcat(saving.directory,'temp_beta_f.bin'),'a');

    fwrite(fileIDA,arrays.accr_f','single');
    fwrite(fileIDB,arrays.beta_f','single');

    fclose(fileIDA);
    fclose(fileIDB);
elseif cur_samp=='f' && tempering.on==0
    arrays.accr_f(finalentry)=single(general.alpha_bar_f);
    arrays.beta_f(finalentry)=single(beta.f);
    arrays.f(general.write_f/general.thin,:)=single(sqrt(diag(samples.f)));

    fileID=fopen(strcat(saving.directory,'d_array.bin'),'a');
    fileID2=fopen(strcat(saving.directory,'av_acc_f.bin'),'a');
    fileID3=fopen(strcat(saving.directory,'beta_f.bin'),'a');

    fwrite(fileID,arrays.f','single');
    fwrite(fileID2,arrays.accr_f','single');
    fwrite(fileID3,arrays.beta_f','single');

    fclose(fileID);
    fclose(fileID2);
    fclose(fileID3);
end

if cur_samp=='p' && tempering.on==1
    arrays.accr_p(finalentry)=single(general.alpha_bar_P);

    fileID=fopen(strcat(saving.directory,'temp_av_acc_P.bin'),'a');
    fwrite(fileID,arrays.accr_p','single');
    fclose(fileID);
elseif cur_samp=='p' && tempering.on==0
    arrays.accr_p(finalentry)=single(general.alpha_bar_P);

    fileID=fopen(strcat(saving.directory,'av_acc_P.bin'),'a');
    fwrite(fileID,arrays.accr_p','single');
    fclose(fileID);
end

if cur_samp=='d' && tempering.on==0
   arrays.cell2matchdist(general.write_f/general.thin,:)=single(samples.cell2matchdist);

   fileID=fopen(strcat(saving.directory,'cell2matchdist.bin'),'a');
   fwrite(fileID,arrays.cell2matchdist','single');
   fclose(fileID);

end
       

end