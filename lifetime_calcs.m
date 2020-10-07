clearvars; clc;

while true
    close all
    promptOne = 'Input temperature = ';
    temperature = input(promptOne,'s');
    if temperature
        
        tempCat = strcat('*K',temperature,'.csv');
        filenames = dir(tempCat);
        data = csvread(filenames.name);
        t = data(:,1); a = data(:,2);
        amin = min(a);
        for i = 1:length(a)
            a(i,1) = a(i,1) - amin;
        end
        amax = max(a);
        for i = 1:length(a)
            a(i,1) = a(i,1) / amax;
        end
        
        %Interactively selects points
%         figure(1) 
        figureFullScreen('Name','Full screen figure size')
%         plot(t,a)
%         pause
%         semilogy(t,a)
        scatter(t,a,'.','SizeData',15)
        set(gca,'yscale','log');
        ylim([0, 10]);
        xlim([-.001,t(end,1)]);
        pause
        [X1(1,1), Y1(1,1)] = getpts(gcf);
%         [X1(2,1), Y1(2,1)] = getpts(gcf);
        close(figure(1),figure(2));
        
        %correcting for t = 0
        tlen = length(t);
        tbli = find(X1(1,1) <= t(:,1));
        for i = tbli(1,1):tlen
            t(i,1) = t(i,1) - t(tbli(1,1),1);
        end
        tbl = t(tbli(1,1):tbli(end,1),1);
        abl = a(tbli(1,1):tbli(end,1),1);
        
        figureFullScreen('Name','Full screen figure size')
        scatter(t,a,'.','SizeData',15)
        set(gca,'yscale','log');
        [X1(2,1),Y1(2,1)] = getpts(gcf);
        abl = abl(tbl <= X1(2,1));
        tbl = tbl(tbl <= X1(2,1));
        close(figure(1))
        
        figureFullScreen('Name','Full screen figure size')
        scatter(t,a,'.','SizeData',15)
        set(gca,'yscale','log');
        [X1(3,1),Y1(3,1)] = getpts(gcf);
        tbl = tbl(abl >= Y1(3,1));
        abl = abl(abl >= Y1(3,1));
        close(figure(1))
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% lsqcurvefit %%%%%%%%%%%%%%%%%%%%%%%%%%
        tblLen = length(tbl);
        
        options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
        lb = [];
        ub = [];
        
        x0 = [1,0.5,0.9615];
        z0 = [1,0.1009,0.9106,0.8001,0.9615];
        
        fun1 = @(x,tbl)paramsfun(x,tbl,tblLen);
        [x,resnormx,residualx,exitflagx,outputx,jacobian] = lsqcurvefit(fun1,x0,tbl,abl,lb,ub,options);
        
        abln = abl - x(3);
        abln_ind = find(abln >= 0);
        abln = abln(abln_ind);
        tbln = tbl(abln_ind);
        abln = abln - min(abln);
        abln = abln ./ max(abln);
        
        [f3,c3] = fit(tbln,abln,'exp1');
        fpar3 = coeffvalues(f3);
        

        
        fun2 = @(z,tbl)dparamsfun(z,tbl,tblLen);
        [z,resnormz,residualz,exitflagz,outputz] = lsqcurvefit(fun2,z0,tbl,abl,lb,ub,options);
        
%         ablnn = abl - z(5);
%         ablnn_ind = find(ablnn >= 0);
%         ablnn = ablnn(ablnn_ind);
%         tblnn = tbl(ablnn_ind);
%         ablnn = ablnn - min(ablnn);
%         ablnn = ablnn ./ max(ablnn);
        
%         [f4,c4] = fit(tblnn,ablnn,'exp2');
%         fpar4 = coeffvalues(f4);
          
        ablnn = abl;
        tblnn = tbl;
        fpar4(1) = z(1);
        fpar4(2) = z(2);
        fpar4(3) = z(3);
        fpar4(4) = z(4);
        fpar4(5) = z(5);
        f4 = fpar4(1).*exp(-fpar4(2).*tblnn) + fpar4(3).*exp(-fpar4(4).*tblnn) + fpar4(5);
        
        %%%%%%%%%%%%%%%%%%%%%%%%% lsqcurvefit END %%%%%%%%%%%%%%%%%%%%%%%
        
        
        %calculating transition rates and lifetimes

        transRate3 = -fpar3(2);
        lifeTime3 = 1 / transRate3 * 1E6;
        transRate41 = fpar4(2);
        lifeTime41 = 1 / transRate41 * 1E6;
        transRate42 = fpar4(4);
        lifeTime42 = 1 / transRate42 * 1E3;
        
        
        %Figures
       
        finfig1 = figure;
%         figureFullScreen(finfig1,'Full screen figure size')
        
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        
        box on
%         
        scatter(tbln,abln,'.','k','SizeData',15)
        set(gca,'yscale','log')
         
        xlabel('Time (Sec)')
        ylabel('Voltage (V)')
        dim = [0.2,0.2,0.258928564244083,0.140476186999253];        
        str = {['A_T_1 = ' num2str(transRate3),' s^{-1}'],['\tau_1 = ' num2str(lifeTime3),' \mus']};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','FontWeight','bold','EdgeColor','black')
        
        leg = findobj(finfig1,'type','legend');
        set(leg,'Visible','off')
        
        grid off
        hold on
%         semilogy(tbl,y3,'r')
        plot(f3,'-k')
        L = findobj(gcf,'type','line');
        set(L,'LineWidth',1.5)
        hold off
        
        pause
        
        
        
        finfig2 = figure;
%         figureFullScreen(finfig2,'Full screen figure size')
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        box on
       
        
        scatter(tblnn,ablnn,'.','k','SizeData',15)
        set(gca,'yscale','log')

        xlabel('Time (Sec)')
        ylabel('Voltage (V)')
        
        
        
        dim = [0.2,0.2,0.258928564244083,0.257142850188982];
        str = {['A_T_1 = ' num2str(transRate41),' s^{-1}'],['\tau_1    = ' num2str(lifeTime41),' \mus'],['A_T_2  = ' num2str(transRate42),' s^{-1}'],['\tau_2  = ' num2str(lifeTime42),' ms']};
        annotation('textbox',dim,'String',str,'FitBoxToText','on','FontWeight','bold','EdgeColor','black')
        
        leg = findobj(finfig2,'type','legend');
        set(leg,'Visible','off')
        
        grid off
        hold on
%         semilogy(tbl,y4,'r')
        plot(f4,'-k')
        L = findobj(gcf,'type','line');
        set(L,'LineWidth',1.5)
        hold off
        

        pause
       
        
        %saving figures
        promptThr = 'Do you want to save figures? ';
        repThr = input(promptThr,'s');
        if repThr
            figfname1 = strcat('NdYSO_site1_lifetime_singexp',temperature,'K.fig');
%             figfname2 = strcat('NdYSO_site1_lifetime_singexp',temperature,'K.png');
            figfname3 = strcat('NdYSO_site1_lifetime_doubexp',temperature,'K.fig');
%             figfname4 = strcat('NdYSO_site1_lifetime_doubexp',temperature,'K.png');
            set(finfig1, 'PaperPositionMode', 'auto');
            set(finfig2, 'PaperPositionMode', 'auto');
            saveas(finfig1,figfname1)
%             saveas(finfig1,figfname2)
            saveas(finfig2,figfname3)
%             saveas(finfig2,figfname4)
            
        end
        
        %saving transition rates vs. temperature for A - T graph and fit
        
        trratesuf = [];
        trrateduf = [];
        trratedff = [];
        promptTwo = 'Do you want to save it to file? ';
        repTwo = input(promptTwo,'s');
        if repTwo
           tempval = str2double(temperature);

           trratesuf(1,1) = tempval;
           trratesuf(1,2) = transRate3;
           dlmwrite('transrate_temp_singexp_nofilt.csv',trratesuf,'delimiter',',','-append');
           trrateduf(1,1) = tempval;
           trrateduf(1,2) = transRate41;
           dlmwrite('transrate_temp_doubexp1_nofilt.csv',trrateduf,'delimiter',',','-append');
           trratedff(1,1) = tempval;
           trratedff(1,2) = transRate42;
           dlmwrite('transrate_temp_doubexp2_nofilt.csv',trratedff,'delimiter',',','-append');
        end
        
    else
        break
    end
end