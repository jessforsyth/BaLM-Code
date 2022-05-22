function run_matching(chain)
%This function now runs alt he aspects of the pipeline, e.g. load in data,
%calculate start temperature, call the function to run tempered sampling,
%sampling at T==1, then figure generation and data save out. 

rng(chain) %specify seed for random number generation as chain number, this can be removed if you don't need reproducible results

tic %start timer

%What sampling do you want to perform 1=yes, 0=no
samp_trans=1; %affine transformation ESSENTIAL must==1
samp_def=0;   %non-linear transformation
samp_fid=1;   %data selection
samp_perm=1;  %permutation sampling ESSENTIAL must==1

%import specific files for each chain (define in import_files)
[data,general, beta, prior, samples, arrays,saving, tempering]=import_files(chain);

%initiate tempering variables
tempering.on=1; %1=include tempering, 0=no tempering
if tempering.on==1
    %calculate start temperature and cooling rate
    [tempering]=tempering_data(general.Nt,tempering, data, prior, general,samp_trans,samp_def,samp_perm,samp_fid);
else
    tempering.T=1;
end

%opening files to overwrite previous data
[saving]=file_generation(saving);
filename=strcat(saving.directory,sprintf('Chain%d_log.txt',chain)); %the Chainx_log.txt file serves as your log file! 
f=fopen(filename,'w');
fwrite(f,sprintf('Chain %d : Sampling on - def : %d, affine : %d, fidelity : %d, permutation : %d  \n', chain,samp_def, samp_trans, samp_fid, samp_perm));
fclose(f);
f=fopen(filename,'a');
fwrite(f,sprintf('Chain %d : temperature=%.3f \n',chain,tempering.T));
fclose(f);

%initialise values of parameters etc. 
[data, general, beta, prior, samples, arrays,saving]=initialise_values(data, general, beta, prior, samples, arrays,saving,samp_trans,samp_def,samp_fid,samp_perm,tempering);

it=0; %initiate iteration counter at 0

%Start tempering
[data, general, beta, prior, samples, arrays,saving,tempering]=run_tempering(it,filename,samp_def,samp_trans,samp_perm,samp_fid,data, general, beta, prior, samples, arrays,saving,tempering) ;

%Perform sampling at T==1
[data, general, beta, prior, samples, arrays,saving,tempering]=run_sampling_T1(filename,samp_def,samp_trans,samp_perm,samp_fid,data, general, beta, prior, samples, arrays,saving,tempering);

time=toc; %time taken for sampling

%save out parameter.m file for future reference 
save(strcat(saving.directory,'parameters.mat'),'data','general','beta','prior','samples','arrays','saving','tempering')  %this saves all of the parameters in a .mat file
        %to access parameters .....
        %parameters=load('parameters.mat')
        %e.g. parameters.data.Y1 

%Plot figures
figure_generation(tempering,data, general, samples, saving, samp_trans, samp_def, samp_fid, samp_perm, prior,arrays, filename) 
        
%final update to Chainx_log.txt file
f=fopen(filename,'a');
fwrite(f,sprintf('\nChain %d : finished in %.3f seconds\n',chain,time));
fclose(f); 

end