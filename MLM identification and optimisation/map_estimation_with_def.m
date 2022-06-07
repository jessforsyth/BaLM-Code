function [X,samples]=map_estimation_with_def(general,data,prior,samples,samp_fid,saving,MLM)

thetafinal=samples.min_theta;
if samp_fid==1
    fidfinal=samples.min_fid';
else
    fidfinal=[];
end
pfinal=reshape(samples.min_p',[1,3*data.n2]); 
X0=double([thetafinal,pfinal,fidfinal]);
    
options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',1e6,'MaxIter',1e6);

if samp_fid==1
    lb=[-pi*ones(1,3),-inf*ones(1,3),-pi*ones(1,3),-inf*ones(1,3),-inf*ones(1,3*data.n2),zeros(1,data.n1)];
    ub=[pi*ones(1,3),inf*ones(1,3),pi*ones(1,3),inf*ones(1,3),inf*ones(1,3*data.n2),ones(1,data.n1)];
else
    lb=[-pi*ones(1,3),-inf*ones(1,3),-pi*ones(1,3),-inf*ones(1,3),-inf*ones(1,3*data.n2)];
    ub=[pi*ones(1,3),inf*ones(1,3),pi*ones(1,3),inf*ones(1,3),inf*ones(1,3*data.n2)];
end

tryagain='y';
while tryagain=='y'
    X = fmincon(@(x)rel_posterior_withdef(x,general,data,prior,samples,samp_fid,MLM),X0,[],[],[],[],lb,ub,[],options);
    X0=X;
    tryagain='n';
    %tryagain=input('Try again from new param mins? Press enter to move onto the next dataset. ','s');
end

samples.opt_min_nlog_posterior=rel_posterior_withdef(X,general,data,prior,samples,samp_fid,MLM);

end