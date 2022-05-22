function multiple_mom_histgen(maxi, mini, data, burn_in,px_array,py_array,pz_array)
%help to make plotting the mutliple histograms easier!
for i=mini:min(maxi,data.n1)
    subplot(4,4,i-mini+1)
    h=histogram(px_array(burn_in:end,i),'Normalization','pdf','DisplayStyle','stairs');
    hold on
    h=histogram(py_array(burn_in:end,i),'Normalization','pdf','DisplayStyle','stairs');
    hold on
    h=histogram(pz_array(burn_in:end,i),'Normalization','pdf','DisplayStyle','stairs');
    ylabel(num2str(i));
end
legend('p_x','p_y','p_z')

end
          