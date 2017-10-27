function problem_3
    function main
        data = gen_data;
        corbanFigureDefaults;
        plot_fig1(data);
        saveAllFigures('prb3');
    end

%% Constants
exp_fun = @(k,t) 0.5 .* exp(-k(1) .* t) ...
    + 0.5 .* exp(-k(2) .* t);
ts = (0:100); % s
ks = {[1, 10]; [10, 20]; [10, 100]}; % 1 / s
n = length(ks);
monoexp_model = @(c,t) c(1) .* exp(-c(2) .* t);
biexp_model = @(c,t) c(1) .* exp(-c(2) .* t) ...
    + c(3) .* exp(-c(4) .* t);

main

%% Implementation

%%% Data Generation
    function data = gen_data
      data = zeros(n*2,length(ts));
      for i = 1:n
          data(i,:) = exp_fun(ks{i} .^ -1, ts);
          data(i+3,:) = exp_fun(ks{i}, ts);
      end
    end

%%% Curve Fitting


%% Figures

    function plot_fig1(data)
        fig = setupFigure(1,'Initial Data Timecourses');
        ax1 = subplot(2,1,1);
        hold on;
        plot([0 100],[0 0],'-k');
        p1 = plot(ts,data(1:3,:),'o-','MarkerSize',5,'LineWidth',1);
        ylabel('Fluorescence, unitless');
        leg_labels = compose('\\tau_{1}, \\tau_{2} = %d, %d', ...
            cell2mat(ks));
        text(ax1,'Interpreter','latex','String',['$$y\left(t\right)='...
            '0.5\exp{\left\{\frac{t}{\tau_1}\right\}}+'...
            '0.5\exp{\left\{\frac{t}{\tau_2}\right\}}$$'],...
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
            '0.5\exp{\left\{k_1t\right\}}+'...
            '0.5\exp{\left\{k_2t\right\}}$$'],...
            'Units','normalized','Position',[0.5 0.8],...
            'HorizontalAlignment','center','FontSize',25)
        legend(p2, leg_labels);
        
        linkaxes([ax1,ax2]);
        ylim([-0.1 1]);
    end
end