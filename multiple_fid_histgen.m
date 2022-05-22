function multiple_fid_histgen(maxi,mini, data, burn_in, fidelity_array,prior,general)
%help to make plotting the mutliple histograms easier!
x=0:.01:1; %points for prior pdf
yprior = betapdf(x,prior.fid_alpha,prior.fid_beta);%prior.fid_C*(prior.fid_a + exp(prior.fid_b*x));
ymaxprior= betapdf(x,prior.fid_alpha+general.d,prior.fid_beta);
for i=mini:min(data.n1,maxi)
    subplot(4,4,i-mini+1)
    hf=histogram(fidelity_array(burn_in:end,i),101,'Normalization','pdf', 'DisplayStyle', 'stairs');
    hold on
    plot(x,yprior)
    hold on
    plot(x,ymaxprior)
    ylabel(num2str(i))
end
sgtitle('Fidelity Terms')
legend('Samples','Prior','Maximum Prior')

end