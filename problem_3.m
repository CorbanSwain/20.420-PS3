function problem_3
    function main
        data = gen_data;
        corbanFigureDefaults;
        f1 = plot_fig1(data);
        [mf, bf] = run_fit(data);
    
        % saveAllFigures('prb3');
    end

%% Constants
exp_fun = @(k,x) 0.5 .* exp(-k(1) .* x) ...
    + 0.5 .* exp(-k(2) .* x);
ts = (0:100); % s
ks = {[1, 10]; [10, 20]; [10, 100]}; % 1 / s
n = length(ks);
monoexp_model = @(c1,c2,x) c1 .* exp(-c2 .* x);
biexp_model = @(c1,c2,c3,c4,x) c1 .* exp(-c2 .* x) ...
    + c3 .* exp(-c4 .* x);
c_mono_0 = [0.5, 10];
c_bi_0 = [0.5, 10, 0.5, 10];
S = load('Noise.mat','All');
noise = S.All';
fit_opt = 2;

clc;
main;

%% Implementation

%%% Data Generation
    function data = gen_data
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
    function [mono_fits, bi_fits] = run_fit(data)
        mono_fits = cell(1,n);
        bi_fits = cell(1,n);
        s1 = struct; s2 = struct;
        x = ts;
        
        fo1 = fitoptions('Method','NonlinearLeastSquares', ...
            'StartPoint', c_mono_0, ...
            'Lower', [-Inf, 0, -Inf, 0]);
        fo2 = fitoptions('Method','NonlinearLeastSquares', ...
            'StartPoint', c_bi_0, ...
            'Lower', [-Inf, 0, -Inf, 0]);
        
        for i = 1:n
            y = data(1,:);
            if fit_opt == 1
                [s1.params, s1.R, ~, ~, ~, ~] = ...
                    nlinfit(x,y, monoexp_model, c_mono_0);
                [s2.params, s2.R, ~, ~, ~, ~] = ...
                    nlinfit(x,y, biexp_model, c_bi_0);
            else
                % TODO - implement again with fit function 
                s1.fitobj = fit(x', y', monoexp_model, fo1); 
                s2.fitobj = fit(x', y', biexp_model, fo2);
            end
            mono_fits{i} = s1;
            bi_fits{i} = s2;
        end 
    end


%% Figures
        function fig = plot_fig1(data)
        fig = setupFigure(1,'Comparison Figure',...
            [-1044 566 957 369]);
        hold on;
        plot([-100 100],[0 0],'-k','LineWidth',0.5);
        plot([0 0],[-1 1],'-k','LineWidth',0.5);
        p1 = plot(ts,data,'o','MarkerSize',5,'LineWidth',1.5);
        ylabel('Fluorescence, unitless');
        xlabel('Time, s')
        leg_labels = compose('k_{1}, k_{2} = %d, %d', ...
            cell2mat(ks));
        ax = gca;
        text(ax,'Interpreter','latex','String',['$$y\left(t\right)='...
            '0.5\exp{\left\{-k_1t\right\}}+'...
            '0.5\exp{\left\{-k_2t\right\}}$$'],...
            'Units','normalized','Position',[0.5 0.85],...
            'HorizontalAlignment','center','FontSize',20)
        legend(p1, leg_labels);
        ylim([-0.2 1]);
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