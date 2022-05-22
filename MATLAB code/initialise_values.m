function [data, general, beta, prior, samples, arrays,saving]=initialise_values(data, general, beta, prior, samples, arrays,saving,samp_trans,samp_def,samp_fid,samp_perm,tempering)
%this function serves to initiate the values of required parameters, this
%could have been done inside each of the functions but seems cleaner to do
%it here. 

if samp_trans==1
    %sample from priors on affine transformation (initialise chain randomly)
    phis1=(rand(1,3)-0.5)*(2*pi());  %affine angles
    phis2=(rand(1,3)-0.5)*(2*pi());  %affine angles 2
    s=normrnd(0,prior.sigs,1,3);     %scaling parameters 
    b=normrnd(0,prior.sigb,1,3);     %translation parameters 
    samples.theta=[phis1,s,phis2,b]; %save to samples structure 
    tphis1= log( ((2*pi)./(phis1+pi)) -1 ); %transform the angle parameters
    tphis2= log( ((2*pi)./(phis2+pi)) -1 );
    [theta_A,theta_b]=matrix_gen(samples.theta); %generate affine matrix and vector
    samples.prior_terma=prior_valsa(samples.theta,prior); %calculate prior density for current values of affine parameters 
    samples.transformation=[tphis1,s,tphis2,b]';  %transformed transformation parameters 
    samples.A=reshape(theta_A',1,[]);  %store the current affine matrix 
    if samp_def==1          
        samples.p =mvnrnd(zeros(3,data.n2),prior.mom_C); %sample random momenta start position 
        %apply deformation
        y0=[single(reshape(data.Y2',3*data.n2,1));reshape(samples.p',3*data.n2,1)];
        [~,yout]=ode45(@(t,y) p_deform_mex(t,y,single(prior.kernel),single(data.n)),single([0,1]),y0);  %tspan=[0,1];
        yfinal=yout(end,:)';
        Y2d=double(reshape(yfinal(1:(3*data.n2)),data.n2,3)'); %deformed Y2 coordinates
        
        samples.prior_termd=prior_vals_p(samples.p',data,prior); %calculate prior term values

        general.num_trans=12 + (3*data.n2) ; %number of parameters in the transformation module    
        samples.transformation=[samples.transformation;samples.p(1,:)';samples.p(2,:)';samples.p(3,:)'];
        %generate the proposal covariance now
        C=zeros(general.num_trans);
        C(1:12,1:12)=diag([ones(1,3).*prior.ang_propstd.^2,ones(1,3).*prior.sigs^2,ones(1,3).*prior.ang_propstd.^2,ones(1,3).*prior.sigb^2]);
        for dim=1:3
            C(13+((dim-1)*data.n2):12+((dim)*data.n2),13+((dim-1)*data.n2):12+((dim)*data.n2))=prior.mom_C;
        end
        %apply deformation and then affine transformation
        samples.fY2=(theta_A*Y2d)+theta_b;        
    else
        samples.fY2=(theta_A*data.Y2)+theta_b; %apply affine transformation (if no non-linear deformation)
        general.num_trans= 12;  %number of parameters in the transformation module  
        %generate the proposal covariance now
        C=diag([ones(1,3).*prior.ang_propstd.^2,ones(1,3).*prior.sigs^2,ones(1,3).*prior.ang_propstd.^2,ones(1,3).*prior.sigb^2]); %calculate C to improve efficiency of sampling on Theta
    end  
    beta.t=(2.38^2)/general.num_trans; %initialise step-size for RW on transformation parameters
    general.C_transformation= C; %covariance matrix for tranformation sampling
    samples.avg_transformation=samples.transformation;
    general.sd=(2.38^2)/general.num_trans;  %adjustment to proposal covariance
end

if samp_fid==1
    f=betarnd(prior.fid_alpha,prior.fid_beta,1,data.n1); %start with random fidelity parameters drawn from prior
    f(f==1)=1-1e-10;
    f(f==0)=1e-10;
    samples.f=diag(f).^2;                             %generate fidelity matrix
    samples.nlog_normconst=fid_norm_const(general,f); %normalisation constant dependent on fidelity params
    samples.prior_termf= prior_valsf_beta(f,prior);   %prior density for current fidelity params
    samples.ftau=log((1./f)-1);                       %position in transformed space
    beta.f=(2.38^2)/(data.n1);                        %initial step size for RW
else 
    samples.f=diag(ones(1,data.n1));                  %if not sampling on fidelity, just set all params equal to 1
end

if samp_perm==1
    samples.P=[randperm(data.n2-data.n_refpoints),linspace(data.n2-data.n_refpoints,data.n2,data.n_refpoints)]; %start at random permutation vector
else
    samples.P=linspace(1,data.n2,data.n2);
end   
 
samples.n_log_PP=neg_log_PP(data.Y1,data.n1,samples.fY2,samples.f,general.Psi,general.nu,samples.P); %calculate the nlog of the posterior predictive at this point for given parameters 
end