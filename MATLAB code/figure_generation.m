function figure_generation(tempering,data, general, samples, saving, samp_trans, samp_def, samp_fid, samp_perm, prior,arrays, filename) 
%save out all of the figures for the current run. saved into saving.directory

if tempering.on==1
    %open files    
    num_entries_temp=general.writing_count_temp * (general.write_f / general.thin);
    
    fileID=fopen(strcat(saving.directory, 'temperature.bin'),'r');
    temp_temps=fread(fileID,[num_entries_temp,1],'single');
    fclose(fileID);
    
    if samp_trans==1
        fileID=fopen(strcat(saving.directory, 'temp_av_acc_t.bin'),'r');
        temp_av_acc_t=fread(fileID,[num_entries_temp,1],'single');
        fclose(fileID);
        fileID=fopen(strcat(saving.directory, 'temp_beta_t.bin'),'r');
        temp_beta_t=fread(fileID,[num_entries_temp,1],'single');
        fclose(fileID);

    else
        temp_av_acc_t=zeros(num_entries_temp,1);
        temp_beta_t=zeros(num_entries_temp,1);
    end
    if samp_perm==1
        fileID=fopen(strcat(saving.directory,'temp_av_acc_P.bin'),'r');
        temp_av_acc_P=fread(fileID,[num_entries_temp,1],'single');
        fclose(fileID);
    else
        temp_av_acc_P=zeros(num_entries_temp,1);
    end
    if samp_fid==1
        fileID=fopen(strcat(saving.directory,'temp_av_acc_f.bin'),'r');
        temp_av_acc_f=fread(fileID,[num_entries_temp,1],'single');
        fclose(fileID);
        fileID=fopen(strcat(saving.directory,'temp_beta_f.bin'),'r');
        temp_beta_f=fread(fileID,[num_entries_temp,1],'single');
        fclose(fileID);
    else
        temp_av_acc_f=zeros(num_entries_temp,1);
        temp_beta_f=zeros(num_entries_temp,1);
    end
    
    iterations=linspace(1,numel(temp_temps),numel(temp_temps))*general.thin;
    %now plot data
    figure('Name', 'Acceptance Rates + Beta Adaptation During Tempering')
    subplot(4,1,1)
    plot(temp_temps)
    set(gca, 'YScale', 'log')
    title('Temperature')
    subplot(4,1,2)
    yyaxis left 
    plot(iterations,temp_av_acc_t);
    ylabel('Avg.Acc. Ratio')
    yyaxis right
    plot(iterations,temp_beta_t)
    set(gca,'YScale','log')
    ylabel('\beta_t')
    title('Transformation Sampling')    
    subplot(4,1,3)
    plot(iterations,temp_av_acc_P)
    set(gca,'YScale','log')
    ylabel('Avg. Acc. Ratio')
    title('Permutation Sampling')    
    subplot(4,1,4)
    yyaxis left
    plot(iterations,temp_av_acc_f)    
    hold on 
    yline(general.acc_opt-general.acc_tol,'--')
    hold on
    yline(general.acc_opt+general.acc_tol,'--')
    ylabel('Avg. Acc. Ratio')
    yyaxis right
    plot(iterations,temp_beta_f)
    set(gca,'YScale','log')
    ylabel('\beta*')
    title('Fidelity Sampling')    
    savefig(gcf,strcat(saving.directory,'acc_rates_beta_temp.fig'),'compact');
       
else 
    %number of samples to write = 
    numsamps=general.writing_count * (general.write_f / general.thin) ; 
    burn_in=1;%round(0.1 * numsamps); %10% burn-in period 
    
    fileID=fopen(strcat(saving.directory,'cell2matchdist.bin'),'r');
    celldist_array=fread(fileID,[data.n1,numsamps],'single');
    celldist_array=celldist_array';
    fclose(fileID);
    
    %open files
    if samp_trans==1
        fileID=fopen(strcat(saving.directory, 'theta_array.bin'),'r');
        theta_array=fread(fileID,[12,numsamps],'single');
        theta_array=theta_array';
        fclose(fileID);
        fileID=fopen(strcat(saving.directory, 'av_acc_t.bin'),'r');
        av_acc_t=fread(fileID,[numsamps,1],'single');
        fclose(fileID);
        fileID=fopen(strcat(saving.directory, 'beta_t.bin'),'r');
        beta_t=fread(fileID,[numsamps,1],'single');
        fclose(fileID);
        fileID=fopen(strcat(saving.directory,'A_matrix.bin'),'r');
        A_matrix=fread(fileID,[9,numsamps],'single');
        A_matrix=A_matrix';
        fclose(fileID);
        fileID=fopen(strcat(saving.directory, 'cov.bin'),'r');
        cov_a=fread(fileID,[general.num_trans,numsamps],'single');
        cov_a=cov_a';
        fclose(fileID);
    else
        theta_array=zeros(numsamps,12);
        av_acc_t=zeros(numsamps,1);
        beta_t=zeros(numsamps,1);
        A_matrix=zeros(numsamps,9);
        cov_a=zeros(numsamps,general.num_trans);
    end
    if samp_def==1
        fileID=fopen(strcat(saving.directory,'px.bin'),'r');
        fileIDa=fopen(strcat(saving.directory,'py.bin'),'r');
        fileIDb=fopen(strcat(saving.directory,'pz.bin'),'r');
        px_array=fread(fileID,[data.n2,numsamps],'single');
        px_array=px_array';
        py_array=fread(fileIDa,[data.n2,numsamps],'single');
        py_array=py_array';
        pz_array=fread(fileIDb,[data.n2,numsamps],'single');
        pz_array=pz_array';
        fclose(fileID);
        fclose(fileIDa);
        fclose(fileIDb);
    else
        px_array=zeros(numsamps,data.n2);
        py_array=zeros(numsamps,data.n2);
        pz_array=zeros(numsamps,data.n2);
    end
    
    if samp_perm==1
        fileID=fopen(strcat(saving.directory,'av_acc_P.bin'),'r');
        av_acc_P=fread(fileID,[numsamps,1],'single');
        fclose(fileID);
        fileID=fopen(strcat(saving.directory,'matches_matrix.txt'),'r');
        matches=fscanf(fileID,'%d', [data.n1,data.n2]);
        fclose(fileID);
    else
        av_acc_P=zeros(numsamps,1);
        matches=zeros(data.n1,data.n2);
    end
    
    if samp_fid==1
        fileID=fopen(strcat(saving.directory,'d_array.bin'),'r');
        fidelity_array=fread(fileID,[data.n1,numsamps],'single');
        fidelity_array=fidelity_array';
        fclose(fileID);
        fileID=fopen(strcat(saving.directory,'av_acc_f.bin'),'r');
        av_acc_f=fread(fileID,[numsamps,1],'single');
        fclose(fileID);
        fileID=fopen(strcat(saving.directory,'beta_f.bin'),'r');           
        beta_f=fread(fileID,[numsamps,1],'single');   
        fclose(fileID);

    else
        fidelity_array=zeros(numsamps,data.n1);
        av_acc_f=zeros(numsamps,1);
        beta_f=zeros(numsamps,1);
    end
    
    iterations=linspace(1,numel(av_acc_t),numel(av_acc_t))*general.thin;

    %now plot data
    figure('Name', 'Acceptance Rates + Beta Adaptation at T=1')
    subplot(4,1,1)
    plot(iterations,ones(numel(av_acc_t),1))
    set(gca, 'YScale', 'log')
    title('Temperature')    
    subplot(4,1,2)
    yyaxis left 
    plot(iterations,av_acc_t);
    ylabel('Avg.Acc. Ratio')
    yyaxis right
    plot(iterations,beta_t)
    set(gca,'YScale','log')
    ylabel('\beta_t')
    title('Transformation Sampling')    
    subplot(4,1,3)
    plot(iterations,av_acc_P)
    set(gca,'YScale','log')
    ylabel('Avg. Acc. Ratio')
    title('Permutation Sampling')
    subplot(4,1,4)
    yyaxis left
    plot(iterations,av_acc_f)
    hold on 
    yline(general.acc_opt-general.acc_tol,'--')
    hold on
    yline(general.acc_opt+general.acc_tol,'--')
    ylabel('Avg. Acc. Ratio')
    yyaxis right
    plot(iterations,beta_f)
    set(gca,'YScale','log')
    ylabel('\beta*')
    title('Fidelity Sampling')    
    savefig(gcf,strcat(saving.directory,'acc_rates_beta_t=1.fig'),'compact');
    
    figure('Name','Affine Transformation Parameters')
    for i=1:9
        subplot(4,3,i)
        hA=histogram(A_matrix(burn_in:end,i),100,'Normalization','pdf');
        ylabel(['A(',num2str(i),')']);
    end
    for j=10:12
        subplot(4,3,j)
        hB=histogram(theta_array(burn_in:end,j),100,'Normalization','pdf');
        ylabel(['b(',num2str(j-9),')']);
    end
    sgtitle('Affine Transformation Parameters')
    savefig(gcf,strcat(saving.directory,'affine_trans_params.fig'),'compact') 
    
    %trace plots of affine transformation params
    figure('Name','Affine Transformation Parameters Trace Plots')
    for i=1:9
        subplot(4,3,i)
        plot(A_matrix(burn_in:end,i));
        ylabel(['A(',num2str(i),')']);
        xlabel('burn_in + it')
    end
    for j=10:12
        subplot(4,3,j)
        plot(theta_array(burn_in:end,j));
        ylabel(['b(',num2str(j-9),')']);
        xlabel('burn_in + it')
    end
    sgtitle('Affine Transformation Parameters')
    savefig(gcf,strcat(saving.directory,'affine_trans_params_trace.fig'),'compact') 
    
    figure('Name', 'Trace of variances at T==1')
    for i=1:general.num_trans
        plot([arrays.cov_update,general.N1],[cov_a(:,i);cov_a(end,i)],'o-')
        hold on
    end
    xlabel('Iteration')
    ylabel('Variance')
    savefig(gcf,strcat(saving.directory,'cov.fig'),'compact');
    
    if samp_fid==1
        %%%%EXTRA FIGURE FOR FIDELITY HISTOGRAMS
 
        figure('Name','Fidelity pdf I')
        multiple_fid_histgen(16,1, data, burn_in, fidelity_array,prior,general);
        savefig(gcf,strcat(saving.directory,'fidelity_terms.fig'),'compact')

        if data.n1>16 
            figure('Name','Fidelity pdf II')
            multiple_fid_histgen(32,17, data, burn_in, fidelity_array,prior,general);
            savefig(gcf,strcat(saving.directory,'fidelity_terms2.fig'),'compact')
        end
        if data.n1>32 
            figure('Name','Fidelity pdf III')
            multiple_fid_histgen(48,33, data, burn_in, fidelity_array,prior,general);
            savefig(gcf,strcat(saving.directory,'fidelity_terms3.fig'),'compact')
        end
        if data.n1>48
            figure('Name','Fidelity pdf IV')
            multiple_fid_histgen(64,49, data, burn_in, fidelity_array,prior,general);
            savefig(gcf,strcat(saving.directory,'fidelity_terms4.fig'),'compact')
        end
    end

    if samp_def==1
        
        figure('Name','Momenta Sampling Histogram 1')
        multiple_mom_histgen(16, 1, data, burn_in,px_array,py_array,pz_array);       
        savefig(gcf,strcat(saving.directory,'momenta.fig'),'compact');
        if data.n1>16 
            figure('Name','Momenta Sampling Histogram 2')
            multiple_mom_histgen(32, 17, data, burn_in,px_array,py_array,pz_array);
            savefig(gcf,strcat(saving.directory,'momenta2.fig'),'compact');
        end
        if data.n1>32 
            figure('Name','Momenta Sampling Histogram 3')
            multiple_mom_histgen(48, 33, data, burn_in,px_array,py_array,pz_array);
            savefig(gcf,strcat(saving.directory,'momenta3.fig'),'compact');
        end
        if data.n1>48 
            figure('Name','Momenta Sampling Histogram 4')
            multiple_mom_histgen(64, 49, data, burn_in,px_array,py_array,pz_array);
            savefig(gcf,strcat(saving.directory,'momenta4.fig'),'compact');
        end
    end
    
    %plot matches and fid heatmap with minimum -log(posterior values)
    
    figure('Name','Point Matching with min -log(posterior) vals')
    Afinal=reshape(samples.min_A,[3,3])';
    bfinal=samples.min_b';
    
    if samp_def==1
        pfinal=reshape(samples.min_p',[3*data.n2,1]);%reshape(MAP_p',[3*data.n2,1]);
        y0=[single(reshape(data.Y2',3*data.n2,1));single(pfinal)];
        [~,yout]=ode45(@(t,y) p_deform_mex(t,y,single(prior.kernel),single(data.n)),single([0,1]),y0);  %tspan=[0,1];
        yfinal=yout(end,:)';
        Y2d=double(reshape(yfinal(1:(3*data.n2)),data.n2,3)');
        %now apply affine
        transY2=(Afinal*Y2d)+bfinal;
    else
        transY2=(Afinal*data.Y2)+bfinal;
    end  
    
    subplot(1,2,1)
    scatter3(data.Y1(1,:),data.Y1(2,:),data.Y1(3,:),60,'MarkerFaceColor',[0.4940 0.1840 0.5560])
    hold on
    scatter3(transY2(1,:),transY2(2,:),transY2(3,:),140,'MarkerEdgeColor',[0.9290 0.6940 0.1250])
    hold on
    scatter3(data.Y2(1,:),data.Y2(2,:),data.Y2(3,:),60,'MarkerFaceColor',[0.9290 0.6940 0.1250])
    hold on
    legend('Y1','Transformed Y2','Y2','Location','southoutside')
    title('Original and Matched Points')
    axis equal    
    subplot(1,2,2)
    colours=linspace(0,1,data.n2);
    %[~,finalP]=max(matches,[],2);
    finalP=samples.min_P;
    for j=1:data.n1
        scatter3(data.Y1(1,j),data.Y1(2,j),data.Y1(3,j),140,colours(j)','filled')
        hold on
        scatter3(transY2(1,finalP(j)),transY2(2,finalP(j)),transY2(3,finalP(j)),140,colours(j));
        hold on
    end
    %now plot remaining 'extra' cells if any present
    lastcells=logical(1-ismember(linspace(1,data.n2,data.n2),finalP)); %extra cells in the case when n1!=n2 ones=cells to plot
    scatter3(transY2(1,lastcells),transY2(2,lastcells),transY2(3,lastcells),60,'k','filled');
    legend('Y1','Transformed Y2','Location','southoutside')
    hold off    
    axis equal
    title(sprintf('Matched Points : min nlogPP=%.4f',samples.min_nlog_posterior))
    savefig(gcf,strcat(saving.directory,'matched_points_using_nlogp.fig'),'compact');
    
    figure('Name','Permutation Sampling Heatmap and Fidelity Parameters ( min nlogP )')
    subplot(1,11,[1,2,3,4,5,6,7,8,9])
    h=heatmap(round(matches(1:data.n1,1:data.n2)./general.N1,2));
    h.Colormap=parula;
    h.XLabel='Cell in Y2';
    h.YLabel='Cell in Y1';
    h.Title='Probability of Matching';
    caxis([0,1])
    subplot(1,11,11)
    h2=heatmap(round(samples.min_fid,2));%heatmap(round(MAP_fidelity',2));
    h2.Colormap=gray;
    h2.Title='Fidelity';
    caxis([0,1])
    savefig(gcf,strcat(saving.directory,'PandFid_heatmap_using_min_nlogp.fig'),'compact');
    
    
    %now make the same plots but using average parameter values
    figure('Name','Point Matching with avg -log(posterior) vals')

    Afinal=reshape(samples.avg_A,[3,3])';
    bfinal=samples.avg_b';
    
    if samp_def==1
        pfinal=reshape(samples.avg_p',[3*data.n2,1]);%reshape(MAP_p',[3*data.n2,1]);
        y0=[single(reshape(data.Y2',3*data.n2,1));single(pfinal)];
        [~,yout]=ode45(@(t,y) p_deform_mex(t,y,single(prior.kernel),single(data.n)),single([0,1]),y0);  %tspan=[0,1];
        yfinal=yout(end,:)';
        Y2d=double(reshape(yfinal(1:(3*data.n2)),data.n2,3)');
        %now apply affine
        transY2=(Afinal*Y2d)+bfinal;
    else
        transY2=(Afinal*data.Y2)+bfinal;
    end  
    
    subplot(1,2,1)
    scatter3(data.Y1(1,:),data.Y1(2,:),data.Y1(3,:),60,'MarkerFaceColor',[0.4940 0.1840 0.5560])
    hold on
    scatter3(transY2(1,:),transY2(2,:),transY2(3,:),140,'MarkerEdgeColor',[0.9290 0.6940 0.1250])
    hold on
    scatter3(data.Y2(1,:),data.Y2(2,:),data.Y2(3,:),60,'MarkerFaceColor',[0.9290 0.6940 0.1250])
    hold on
    legend('Y1','Transformed Y2','Y2','Location','southoutside')
    title('Original and Matched Points')
    axis equal    
    subplot(1,2,2)
    colours=linspace(0,1,data.n2);
    %[~,finalP]=max(matches,[],2);
    finalP=samples.min_P;
    for j=1:data.n1
        scatter3(data.Y1(1,j),data.Y1(2,j),data.Y1(3,j),140,colours(j)','filled')
        hold on
        scatter3(transY2(1,finalP(j)),transY2(2,finalP(j)),transY2(3,finalP(j)),140,colours(j));
        hold on
    end
    %now plot remaining 'extra' cells if any present
    lastcells=logical(1-ismember(linspace(1,data.n2,data.n2),finalP)); %extra cells in the case when n1!=n2 ones=cells to plot
    scatter3(transY2(1,lastcells),transY2(2,lastcells),transY2(3,lastcells),60,'k','filled');
    legend('Y1','Transformed Y2','Location','southoutside')
    hold off    
    axis equal
    title(sprintf('Matched Points : avg nlogPP=%.4f',samples.avg_nlog_posterior))
    savefig(gcf,strcat(saving.directory,'matched_points_using_avglogp.fig'),'compact');
    
    figure('Name','Permutation Sampling Heatmap and Fidelity Parameters (avg nlogP)')
    subplot(1,11,[1,2,3,4,5,6,7,8,9])
    h=heatmap(round(matches(1:data.n1,1:data.n2)./general.N1,2));
    h.Colormap=parula;
    h.XLabel='Cell in Y2';
    h.YLabel='Cell in Y1';
    h.Title='Probability of Matching';
    caxis([0,1])
    subplot(1,11,11)
    h2=heatmap(round(samples.avg_fid,2));%heatmap(round(MAP_fidelity',2));
    h2.Colormap=gray;
    h2.Title='Fidelity';
    caxis([0,1])
    savefig(gcf,strcat(saving.directory,'PandFid_heatmap_using_avglogp.fig'),'compact');
    
    
    %minimum nlog posterior to file along with associated values. 
    f=fopen(filename,'a');
    fwrite(f,sprintf('\n\nMinimum negative log posterior: %f \n', samples.min_nlog_posterior));
    fwrite(f,sprintf('Corresponding P:\n'));
    fwrite(f,sprintf('%d,', samples.min_P));
    fwrite(f,sprintf('\nCorresponding A:\n'));
    fwrite(f,sprintf('%f,', samples.min_A));
    fwrite(f,sprintf('\nCorresponding b:\n'));
    fwrite(f,sprintf('%f,', samples.min_b));
    fwrite(f,sprintf('\n'));
    if samp_def==1
        fwrite(f,sprintf('\nCorresponding p:\npx:'));
        fwrite(f,sprintf('%f,', samples.min_p(1,:)));
        fwrite(f,sprintf('\npy:'));
        fwrite(f,sprintf('%f,', samples.min_p(2,:)));
        fwrite(f,sprintf('\npz:'));
        fwrite(f,sprintf('%f,', samples.min_p(3,:)));
        fwrite(f,sprintf('\n'));
    end
    if samp_fid==1
        fwrite(f,sprintf('Corresponding fidelity terms:\n'));
        fwrite(f,sprintf('%f,', samples.min_fid'));
    end
    fclose(f);
    
    %now write down average nlog posterior and associated averages
    f=fopen(filename,'a');
    fwrite(f,sprintf('\n \n \nAverage negative log posterior: %f ', samples.avg_nlog_posterior));

    fwrite(f,sprintf('\nCorresponding A:\n'));
    fwrite(f,sprintf('%f,', samples.avg_A));
    fwrite(f,sprintf('\nCorresponding b:\n'));
    fwrite(f,sprintf('%f,', samples.avg_b));
    fwrite(f,sprintf('\n'));
    if samp_def==1
        fwrite(f,sprintf('\nCorresponding p:\npx:'));
        fwrite(f,sprintf('%f,', samples.avg_p(1,:)));
        fwrite(f,sprintf('\npy:'));
        fwrite(f,sprintf('%f,', samples.avg_p(2,:)));
        fwrite(f,sprintf('\npz:'));
        fwrite(f,sprintf('%f,', samples.avg_p(3,:)));
        fwrite(f,sprintf('\n'));
    end
    if samp_fid==1
        fwrite(f,sprintf('Corresponding fidelity terms:\n'));
        fwrite(f,sprintf('%f,', samples.avg_fid'));
    end
    fclose(f);
    
    figure()
    %plotcell to match distancehistograms
    for i=1:data.n1
        histogram(celldist_array(burn_in:end,i),101,'Normalization','pdf', 'DisplayStyle', 'stairs');
        hold on
    end
    xlabel('Cell to match Distance')
    ylabel('P(d)')
    savefig(gcf,strcat(saving.directory,'celltomatchdist.fig'),'compact');
    
    
end

end