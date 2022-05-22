function [general]=check_covariance(general,saving,it)
%check if covariance is positive semi-definite and symmetric (helps accounts for
%rounding errors)
    symtf=issymmetric(general.C_transformation);
    if symtf==0 
        %if C is not symmetric
        f=fopen(strcat(saving.directory,'covariance_log.txt'),'a'); %write to C log file
        fwrite(f,sprintf('Iteration : %d \n Error: The covariance matrix is not symmetric. Performed 0.5(C + CT) to adjust \n',it));
        fclose(f);
        general.C_transformation=0.5.*(general.C_transformation + general.C_transformation'); %make C symmetric
        general.C_transformation_temp=general.C_transformation; %and assign
    end
    
    smeig=min(eig(general.C_transformation));%calculate the smallest eigenvalue of C
    if smeig<=0 %if we have negative eigenvalues then 
        f=fopen(strcat(saving.directory,'covariance_log.txt'),'a'); %write to C log file
        fwrite(f,sprintf('Iteration : %d \n Error: The covariance matrix has negative eigenvalues and is not positive semi-definite. abs(min eigenvalue %f x 1.05) added to diagonals. \n',it,abs(smeig)));
        fclose(f);
        general.C_transformation = general.C_transformation + ((1.05 * abs(smeig)).*eye(general.num_trans)); %add on small value (1.05*smeig) to make C positive semi-def
    end

    general.C_transformation_temp=general.C_transformation; %re-assign

end