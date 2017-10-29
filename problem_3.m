function problem_3
    function main
        data = gen_data;
        corbanFigureDefaults;
        fits = run_fit(data);
        f2 = plot_datafit_fig(data, fits, 2);
        saveFigure(f2, 'prb3-');
        stats = calc_fit_stats(fits);
        ftest_out = ftest(stats);
        print_stats("No Noise", stats, ftest_out);
        
        data = add_noise(data);
        fits = run_fit(data);
        f3 = plot_datafit_fig(data, fits, 3);
        saveFigure(f3, 'prb3-Noise');
        stats = calc_fit_stats(fits);
        ftest_out = ftest(stats);
        print_stats("Noise", stats, ftest_out);
    end

%% Global Variables
ts = (0:100); % s
ks = {[1, 10]; [10, 20]; [10, 100]}; % 1 / s
n = length(ks);
% for i__ = 1:n; ks{i__} = ks{i__} .^ -1; end
monoexp_model = @(c,x) c(1) .* exp(-c(2) .* x);
biexp_model = @(c,x) c(1) .* exp(-c(2) .* x) ...
    + c(3) .* exp(-c(4) .* x);
fit_opt = 2; % 1 for nlinfit(...) 2 for fit(...)

% Noise
S = load('Noise.mat','All');
noise = S.All';

% Figure GLobal Var
num_runs = 0;

%% Run Script
clc;
diary('problem_3_diary.txt');
disp(datestr(now));
fprintf('vvvvvvvvvvvvvvvv BEGINNING SCRIPT vvvvvvvvvvvvvvvvvvvv\n\n');
main;
fprintf('\n\n^^^^^^^^^^^^^^^^^^ ENDING SCRIPT ^^^^^^^^^^^^^^^^^^^^^\n\n');
diary off;

%% Implementation

%%% Data Generation
    function data = gen_data
        exp_fun = @(k,x) 0.5 .* exp(-k(1) .* x) ...
            + 0.5 .* exp(-k(2) .* x);
        data = zeros(n,length(ts));
        for i = 1:n
            data(i,:) = exp_fun(ks{i}, ts);
        end
    end

%%% Addding Noise
    function noisy_data = add_noise(data)
        noisy_data = data + noise;
    end

%%% Curve Fitting
    function [fits] = run_fit(data)
        fits = cell(2,n);
                
        % Fit Paramerers
        if fit_opt == 1
            c_mono_0 = cell(1,n);
            c_bi_0 = cell(1,n);
            for i = 1:n
                c_mono_0{i} = [1, min(ks{i})];
                c_bi_0{i} = [1, ks{i}(1), 1, ks{i}(2)];
            end
        else
            monoexp_model_v2 = @(c1, c2, x) ...
                monoexp_model([c1, c2], x);
            biexp_model_v2 = @(c1, c2, c3, c4, x) ...
                biexp_model([c1, c2, c3, c4], x);
            fo1 = fitoptions('Method','NonlinearLeastSquares', ...
                'StartPoint', ones(1,2), ...
                'Lower', [-Inf, 0]);
            fo2 = fitoptions('Method','NonlinearLeastSquares', ...
                'StartPoint', ones(1,4), ...
                'Lower', [-Inf, 0, -Inf, 0]);
        end
        s1 = struct; s2 = struct;
        x = ts;
        for i = 1:n
            y = data(i,:);
            if fit_opt == 1
                [s1.params, s1.R, s1.J, ~, ~, ~] ...
                    = nlinfit(x,y, monoexp_model, c_mono_0{i});
                [s2.params, s2.R, s2.J, ~, ~, ~] ...
                    = nlinfit(x,y, biexp_model, c_bi_0{i});
            else
                [s1.fitobj, s1.gof, s1.output] ...
                    = fit(x', y', monoexp_model_v2, fo1); 
                s1.params = coeffvalues(s1.fitobj);
                [s2.fitobj, s2.gof, s2.output] ...
                    = fit(x', y', biexp_model_v2, fo2);
                s2.params = coeffvalues(s2.fitobj);
            end
            fits{1,i} = s1;
            fits{2,i} = s2;
            
        end 
    end

%%% Determining Fit Statistics
    function stats = calc_fit_stats(fits)
        stats = cell(2,n);
        for i = 1:2
            for j = 1:n
                if fit_opt == 2
                    s = struct;
                    fit_struct = fits{i, j};
                    fitobj = fit_struct.fitobj;
                    output = fit_struct.output;
                    gof = fit_struct.gof;
                    % Parameters
                    s.params = coeffvalues(fitobj);
                    s.conf = confint(fitobj);
                    s.rsq = gof.rsquare;
                    s.res = output.residuals;
                    s.ss = gof.sse;
                    s.df = output.numobs - output.numparam;
                    stats{i,j} = s;
                end
            end
        end
    end

%%% F-Test
    function results = ftest(stats)
        results = cell(1,n);
        for i = 1:n
            s = struct;
            ss_null = stats{1,i}.ss;
            ss_alt = stats{2,i}.ss;
            df_null = stats{1,i}.df;
            df_alt = stats{2,i}.df;
            s.f_stat = ((ss_null - ss_alt) ./ ss_alt) ...
                ./ ((df_null - df_alt) ./ df_alt);
            s.p_val = fcdf(s.f_stat,df_null,df_alt);
            results{i} = s;
        end
    end

%%% Displaying Values
    function print_stats(descr, stats, ftest_out)
        disp(descr);
        for i = 1:n
            [sm, sb] = stats{:,i};
            ft = ftest_out{i};
            fprintf('\n**************\n');
            fprintf('Model #1: k_1 = %d, k_2 = %d\n', ks{i});
            
            mono = cell(3,1);
            mono{1} = sprintf('%.2f, %.2f',sm.params);
            mono{2} = sprintf('(%.2e, %.2e), (%.2e, %.2e)', ...
                sm.conf(:,1), sm.conf(:,2));
            mono{3} = sprintf('%.2f',sm.rsq);
            
            bi = cell(3,1);
            bi{1} = sprintf('%.2f, %.2f, %.2f, %.2f',sb.params);
            bi{2} = sprintf(['(%.2e, %.2e), (%.2e, %.2e),' ...
                ' (%.2e, %.2e), (%.2e, %.2e)'], ...
                sb.conf(:,1), sb.conf(:,2),...
                sb.conf(:,3), sb.conf(:,4));
            bi{3} = sprintf('%.2f',sb.rsq);
            T = table(mono,bi);
            T.Properties.RowNames = {'Params'; 'Confidence Intervals'; ...
                'R^2'};
            T.Properties.VariableNames = {'Monoexp_Fit';'Biexp_Fit'};
            disp(T)
            fprintf('F-test: %.2f, p: %2f\n\n',ft.f_stat, ft.p_val);
            
            %{
            
            Model 1: k_1 = ___, k_2 = ___
            
            Model Sampled at 1/s:
            
                                     MONO    |   BI
            
            PPARAMETERS                      |
            CONFIDENCE INTERVALS
            R^2
            
            F Ration Test:
            F VALUE
            P VALUE

            %}
            
            
        end
    end

%% Figures  
    function fig = plot_fig1(data)
        same_screen_pos = [666 551 957 369];
        diff_screen_pos = [-1044 566 957 369];
        fig = setupFigure(1,'Comparison Figure',...
            same_screen_pos);
        plot_data_points(data)
        ax = gca;
        text(ax,'Interpreter','latex','String',['$$y\left(t\right)='...
            '0.5\exp{\left\{-k_1t\right\}}+'...
            '0.5\exp{\left\{-k_2t\right\}}$$'],...
            'Units','normalized','Position',[0.5 0.85],...
            'HorizontalAlignment','center','FontSize',20)
    end
    
    function fig = plot_datafit_fig(data, fits, num)
        num_runs = num_runs + 1;
        locs = {[44 409 959 532]; [679 58 959 532]};
        fig = setupFigure(num, 'Data With Fit', locs{num_runs});
        plot_data_points(data);
        plot(-10, -10, '-k','DisplayName','Mono-exponential Fit');
        plot(-10, -10, '--k','DisplayName','Bi-exponential Fit');
        legend('AutoUpdate', 'off');
        ax = gca;
        if true; x = -10:0.1:100; end
        for i = 1:n
            [mf, bf] = fits{:,i};
            if true
                ax.ColorOrderIndex = i;
                y = monoexp_model(mf.params,x);
                plot(x,y,'-');
%                 disp_eqn(eqn_str,[0.1]) 
                ax.ColorOrderIndex = i;
                y = biexp_model(bf.params,x);
                plot(x,y,'--');
            else
                ax.ColorOrderIndex = i;
                plot(mf.fitobj,'-');
%                 disp_eqn(eqn_str,[0.1 
                ax.ColorOrderIndex = i;
                plot(bf.fitobj,'--');
            end
        end

        
    end

    function plot_data_points(data)
        plot([-100 100],[0 0],'-k','LineWidth',0.5);
        plot([0 0],[-1 1],'-k','LineWidth',0.5);
        p1 = plot(ts,data,'o','MarkerSize',5,'LineWidth',1.5);
        ylabel('Fluorescence, unitless');
        xlabel('Time, s')
        leg_labels = compose('k_{1}, k_{2} = %.2f, %.2f', ...
            cell2mat(ks));
        legend(p1, leg_labels);
        ylim([-0.3 1]);
        xlim([-3.5 100]);
    end

    function fig = plot_comparison_fig(data)
        fig = setupFigure(1,'Comparison Figure',...
            [-1075 230 988 705]);
        ax1 = subplot(2,1,1);
        hold on;
        plot([0 100],[0 0],'-k');
        p1 = plot(ts,data(1:3,:),'o-','MarkerSize',5,'LineWidth',1);
        ylabel('Fluorescence, unitless');
        leg_labels = compose('\\tau_{1}, \\tau_{2} = %d, %d', ...
            cell2mat(ks));
        text(ax1,'Interpreter','latex','String',['$$y\left(t\right)='...
            '0.5\exp{\left\{-\frac{t}{\tau_1}\right\}}+'...
            '0.5\exp{\left\{-\frac{t}{\tau_2}\right\}}$$'],...
            'Units','normalized','Position',[0.5 0.8],...
            'HorizontalAlignment','center','FontSize',25)
        legend(p1, leg_labels)
        
        ax2 = subplot(2,1,2);
        hold on;
        plot([0 100],[0 0],'-k');
        p2 = plot(ts,data(4:6,:),'o-','MarkerSize',5,'LineWidth',1);
        ylabel('Fluorescence, unitless');
        xlabel('Time, s')
        leg_labels = compose('k_{1}, k_{2} = %d, %d', ...
            cell2mat(ks));
        text(ax2,'Interpreter','latex','String',['$$y\left(t\right)='...
            '0.5\exp{\left\{-k_1t\right\}}+'...
            '0.5\exp{\left\{-k_2t\right\}}$$'],...
            'Units','normalized','Position',[0.5 0.8],...
            'HorizontalAlignment','center','FontSize',25)
        legend(p2, leg_labels);
        
        linkaxes([ax1,ax2]);
        ylim([-0.1 1]);
    end

end