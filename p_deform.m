function [dydt]=p_deform(t,y,sigkernel, n)
%Function to perform the non-linear deformation using Y2 and original
%momentum for each point. 

%unpack positions and momenta
q=y(1:3*n,1);       %position
q=reshape(q,n,3)';

p=y((3*n)+1:end,1); %momenta
p=reshape(p,n,3)';

%create empty derivative vectors
dq=single(zeros(3,n));
dp=single(zeros(3,n));

for cellj=1:n
    ut_qtj=single(zeros(3,1));      %empty array for storing vals
    grad_ut_qtj=single(zeros(3,1)); %empty array for storing vals
    for celli=1:n
        
        diff=q(:,celli)-q(:,cellj);                 %refer to 2.1 and S1.5 in the manuscript for detailed derivation of approach 
        K=exp(-((norm(diff)^2)/(2*sigkernel^2)));
        Kp=(K*p(:,celli));
        
        ut_qtj= ut_qtj + Kp ;
        grad_ut_qtj= grad_ut_qtj + ((diff/(sigkernel^2)).*Kp);
    end
    
    dp(:,cellj)=dot(-grad_ut_qtj', p(:,cellj) ); 
    dq(:,cellj)=ut_qtj;
end

%repack
dydt=[reshape(dq',3*n,1);reshape(dp',3*n,1)];

end
        
