%post-processing estimation of maximum likely match (MLM)
%then find MAPs on transformation and fidelity parameters conditioned on
%the MLM

clc
clear

close all 

samp_fid=0;
samp_def=0;
def_test=0;

startingFolder = pwd;
root = uigetdir(startingFolder);
allSubfolders = genpath(root);
subFolders = regexp(allSubfolders, ';', 'split'); %first folder is the original selected folder

diary(strcat(root,'\MAPest_logs.txt'));

for k = 2 :  length(subFolders)-1
    % Get this subfolder.
	thisSubFolder = subFolders{k}
    %now load parameters file
    load(fullfile(thisSubFolder,'\parameters.mat'))
    saving.directory=thisSubFolder;
    
    %load matches file and find most likely matches (MLM)
    fileID=fopen(fullfile(thisSubFolder,'matches_matrix.txt'),'r');
    matches=fscanf(fileID,'%d', [data.n1,data.n2]);
    fclose(fileID);
    
    costM=1-(matches./general.N1);
    MLM=matchpairs(costM',10); %identify most likely match matrix using matchpairs algorithm
    %for each column in A find rows, so we need to use transpose of our matches matrix
    MAPs.P=MLM(:,1)';
    
    %MAP estimation
    if samp_def==0
        
        [X,samples]=map_estimation_no_def(general,data,prior,samples,samp_fid,saving,MLM);
        MAPs.params=X;
        [MAPs.A,MAPs.b]=matrix_gen(X(1:12));
        MAPs.fid=X(13:end);
        %MAPs.P=samples.min_P;

        MAPs.min_nlog_posterior=samples.opt_min_nlog_posterior;
        MAPs.Y2=(MAPs.A*data.Y2)+MAPs.b;
        save(strcat(saving.directory,'\MAPS.mat'),'MAPs') 
        
    else
        [X,samples]=map_estimation_with_def(general,data,prior,samples,samp_fid,saving,MLM);
        MAPs.params=X;
        [MAPs.A,MAPs.b]=matrix_gen(X(1:12));
        MAPs.p=X(13:12+(3*data.n2));
        %MAPs.P=samples.min_P;
        MAPs.fid=X(12+(3*data.n2)+1:end);

        MAPs.min_nlog_posterior=samples.opt_min_nlog_posterior;
        %apply def
        y0=[single(reshape(data.Y2',3*data.n2,1));single(MAPs.p')];
        [~,yout]=ode45(@(t,y) p_deform_mex(t,y,single(prior.kernel),single(data.n)),single([0,1]),y0);  %tspan=[0,1];
        yfinal=yout(end,:)';
        Y2d=double(reshape(yfinal(1:(3*data.n2)),data.n2,3)');
        %perform affine transformation

        MAPs.Y2=(MAPs.A*Y2d)+MAPs.b;
        save(strcat(saving.directory,'\MAPS.mat'),'MAPs') 
    end    
    
end

subFolders'

diary off