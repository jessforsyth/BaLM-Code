function [data, general, beta, prior, samples, arrays,tempering]=transformation_sampling(it,samp_def,data, general, beta, prior, samples, arrays,saving,tempering)
%this function works to sample on the 12 parameters that describe a 3D
%affine trasnformation (and deformation momenta if samp_def==1). We use adaptive Metropolis Hastings
%MCMC, and adapt beta (step size) during sampling. 

%First check the acceptance rate is within tolerance limits. 
if mod(it,general.fc)==0 %only check if at some mlutiple of the check frequency
    general.alpha_t(general.fc)=samples.alpha_t; %final alpha value to array
    general.alpha_bar_t=mean(general.alpha_t);   %calculate the average value of the acceptance ratios for previous fc iterations
    
    if general.alpha_bar_t > (general.acc_opt + general.acc_tol)
        beta.t=beta.t*(1+beta.adj); %accrate too high, increase step size
    elseif general.alpha_bar_t < (general.acc_opt - general.acc_tol)
        beta.t=beta.t*(1-beta.adj); %accrate too low, decrease step size
    elseif tempering.on==1 && general.alpha_bar_t > (general.acc_opt - general.acc_tol) && general.alpha_bar_t < (general.acc_opt + general.acc_tol)
        %if in the correct accrate regime, and tempering , reduce temp
        tempering.T=max(1,tempering.T*(tempering.tc));
        beta.t=beta.t*tempering.betac; %adjust beta in response to temp decrease
        beta.f=beta.f*tempering.betac; %adjust fidelity beta as well!
    end
else
    general.alpha_t(mod(it,general.fc))=samples.alpha_t; %if not at check frequency, just append to array
end

%perform random walk proposal (on unbounded params i.e. transformed params)
v = samples.transformation' + (beta.t .* mvnrnd(zeros(1,general.num_trans),general.C_transformation)) ;
v_theta=v(1:12);  %proposed affine transformation parameters
if samp_def==1
    v_p=v(13:end);%proposed momenta
end

%convert angular components of v from -inf -->inf variables back to bounded angles .... 
to_check=[1,2,3,7,8,9];
v_angles=((2*pi)./(exp(v_theta(to_check))+1)) - pi;
v_theta_ang=v_theta;
v_theta_ang(to_check)=v_angles ;

%generate matrices and vector for affine transformation
[theta_A,theta_b]=matrix_gen(v_theta_ang);

%if performing deformation sampling - do this here using the function below...... 
if samp_def==1
    %apply deformation
    y0=[single(reshape(data.Y2',3*data.n2,1));v_p'];
    [~,yout]=ode45(@(t,y) p_deform_mex(t,y,single(prior.kernel),single(data.n)),single([0,1]),y0);  %tspan=[0,1];
    yfinal=yout(end,:)';
    Y2d=double(reshape(yfinal(1:(3*data.n2)),data.n2,3)');
    %perform affine transformation
    fY2d=(theta_A*Y2d)+theta_b;
    %calculate components of acceptance ratio 
    priortermd=prior_vals_p(reshape(v_p,[data.n2,3]),data,prior); %prior density on momenta
    n_log_PP=neg_log_PP(data.Y1,data.n1,fY2d,samples.f,general.Psi,general.nu,samples.P); %i.e. calculate the nlogPP with deformed coords
    priorterma=prior_valsa(v_theta_ang,prior); %prior density on affine parameters (remember to use angles not transformed params
    all_priors=samples.prior_terma - priorterma + samples.prior_termd - priortermd; %compute entire prior density comparing proposed and current densities 
else
    %perform affine transformation
    fY2=(theta_A*data.Y2)+theta_b;
    %calculate components of alpha
    n_log_PP=neg_log_PP(data.Y1,data.n1,fY2,samples.f,general.Psi,general.nu,samples.P); %i.e. calculate the nlogPP with non-deformed coords
    priorterma=prior_valsa(v_theta_ang,prior);
    all_priors=samples.prior_terma - priorterma; %compare prior densities for current and previous affine parameter sets
end

log_rp= reference_pullback(v(to_check),samples.transformation(to_check)');  %to account for the sampling on the transformed params rather than the bounded parameters 

%now calculate acceptance ratio (with tempered posterior pred.)
alpha=min(1, exp(((1/tempering.T)*(samples.n_log_PP - n_log_PP)) + all_priors + log_rp));

if alpha>rand()
    %accept parameters
    samples.A=reshape(theta_A',1,[]);
    samples.theta=v_theta_ang;
    samples.n_log_PP=n_log_PP;
    samples.prior_terma=priorterma;
    samples.transformation=v';
    
    if samp_def==1
        samples.p=reshape(v_p,[data.n2,3])';
        samples.fY2=fY2d;
        samples.prior_termd=priortermd;
    else
        samples.fY2=fY2;
    end
end

samples.alpha_t=alpha;

%append values to arrays and save data
if mod(it,general.write_f)==0  
    %we only write data every 10e3 iterations, but you can change this in the parameter.m file 
    general.writing_count=general.writing_count+1;
    write_your_data('a',saving, general, beta, tempering,samp_def,arrays,samples,data);
elseif  mod(it,general.thin)==0 && mod(it,general.write_f)~=0 
    %if you aren't writing data, but you are storing this current iteration
    %(i.e. you're not on a thinned iteration) store vals in arrays
    itap=(it/general.thin)-(floor(it/general.write_f)*(general.write_f/general.thin));
    arrays.theta(itap,:)=single(samples.theta);
    arrays.A(itap,:)=single(samples.A);
    arrays.temp(itap)=single(tempering.T);
    arrays.accr_t(itap)=single(general.alpha_bar_t);
    arrays.beta_t(itap)=single(beta.t);
    if samp_def==1
        arrays.px(itap,:)=single(samples.p(1,:));
        arrays.py(itap,:)=single(samples.p(2,:));
        arrays.pz(itap,:)=single(samples.p(3,:));
    end       
end

%adaptive MCMC ONLY WHEN T==1
if tempering.on==0 && it>2 
    [general]=trans_t1covariance_gen(samples, general, it);
    if it==general.check_C %if we are checking the covariance
        if general.adaptive_on==0 
            %if this is the first time we are updating the covariance to an adaptively learned covariance -re-initialise beta.t one first adaptation
            beta.t=(2.38^2)/general.num_trans;   %
        end
        general.adaptive_on=1; %switch adaptive MCMC on
        general.C_transformation=general.C_transformation_temp; %update the proposal covariance to the learned covariance        
        [general]=check_covariance(general,saving,it); %check that the learned covariance is symmetric & pos semi-def
        
        %write to file
        fileIDD=fopen(strcat(saving.directory,'cov.bin'),'a');
        fwrite(fileIDD,diag(general.C_transformation),'single');
        fclose(fileIDD);
        arrays.cov_update=[arrays.cov_update,it];
                
        general.check_C=round(general.check_C * general.check_adj);
    end   
end

%perform rolling avg of parameter values here so we always have avg_n-1 stored for covraince estimates
samples.avg_transformation=samples.avg_transformation + (samples.transformation - samples.avg_transformation)./it; %only update avg here as we use avg(it-1) for calc of covariance


end


