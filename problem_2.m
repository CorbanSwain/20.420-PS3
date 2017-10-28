function problem_2
%{
Problem Statement 

2. You are interested in experimentally determining the ligand binding affinity of a cell
surface receptor with Kd =10 pM and association rate constant kon =10^5 M?1s?1. You incubate
nine separate cell suspensions with the following concentrations of fluorescently labeled
ligand: 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, and 1,000.0 pM. Using Matlab OR
Python (your choice), simulate the binding isotherm (binding curve) you would obtain if
you measured cell-bound fluorescence at the following times: 15 minutes; 1 hour; 4
hours; 16 hours; 64 hours. You may assume pseudo-first order kinetics (i.e. the ligand is
not used up by binding). Curve-fit each of these five isotherms to obtain an estimate of Kd.
Is it necessary for every sample to reach equilibrium in order to obtain an accurate Kd
estimate?
%}
    function main
        data = run_simulation;
        
        corbanFigureDefaults;
        % plot_fig1(data);      
        % plot_fig2(data);
        
        fits = curve_fit(data);
        f3 = plot_fig3(data,fits);
        % FIXME - Add residual plot
        saveAllFigures('prb2');

    end

%% Constants
lig_concs = [ ...
    0.1, 0.3, 1, ...
    3, 10, 30, ...
    100, 300, 1000]; % pM
timepoints = [ ...
    0.25, 1, 4, 16, 64]; % hours
Kd = 10; % pM
k_on = 1e5 * 1e-12 * 60 * 60; % 1 / (pM h)
binding_func = @(p,L) L ./ (L + p);
n = length(timepoints);

main;

%% Data Generation
    function data = run_simulation
        [m_lig_concs, m_timepoints] = ...
            meshgrid(lig_concs, timepoints);
        c1 = m_lig_concs ./ (Kd + m_lig_concs);
        c2 = (Kd + m_lig_concs) .* k_on;
        data = c1 .* (1 - exp(-c2 .* m_timepoints));
    end

%% Curve Fitting
    function fits = curve_fit(data)
        p_0 = 1;
        fits = cell(1,n);
        for i = 1:n
            fits{i} = nlinfit(lig_concs,data(i,:),binding_func,p_0);
            % fits{i} = fit(lig_concs,data(1,:
        end
    end

%% Plotting

    function plot_data_pts(data)
        data = flipud(data);
        p = plot(log10(lig_concs),data,'.');
        xlabel('log[L], log(pM)');
        ylabel('Relative Fluorescence, unitless');
        legend_vals = strcat(num2str(timepoints')," Hours");
        legend_vals(1) = "15 Minutes";
        legend_vals = flipud(legend_vals);
        legend(p, legend_vals,'Location','northwest',...
            'AutoUpdate','off');
        
    end

    function fig = plot_fig1(data)
        fig = etupFigure(1,'Binding Isotherms');
        plot_data_pts(data);

    end

    function fig = plot_fig2(data)
        fig = setupFigure(2,'Binding Timecourses');
        p = plot(timepoints,data,'.-');
        xlabel('Time, hours');
        ylabel('Relative Fluoroesence, unitless')
        legend_vals = strcat(num2str(lig_concs')," pM");
        legend(legend_vals,'Location','northoutside',...
            'Orientation','horizontal');
    end

    function fig = plot_fig3(data,fits)
        fig = setupFigure(3,"Overlayed Fits",[-1691 489 970 396]);
        plot_data_pts(data);
        ax = gca;
        ax.ColorOrderIndex = 1;
        x = logspace(log10(lig_concs(1)),...
            log10(lig_concs(end)), 100);
        text_x_posns = [2.75, 2.75, 2.75, 2.75, 2.75]-0.25;
        y_fudge = [0.06, 0.14, 0.2, 0.13, 0.08];
        for i = 1:n
            ind = n - i + 1;
            fit_Kd = fits{ind};
            y = binding_func(fit_Kd,x);
            p = plot(log10(x),y,'-','LineWidth',1.5);
            text_x_pos = text_x_posns(ind);
            text_y_pos = binding_func(fit_Kd, 10 .^ text_x_pos);
            del_x = 0.5;
            del_y = -text_y_pos ...
                + binding_func(fit_Kd, 10 .^ (text_x_pos + del_x));
            angle = atan(del_y ./ del_x) .* 180 ./ pi;
            exp_val = floor(log10(fit_Kd) ./ 3) .* 3;
            bas_val = fit_Kd ./ (10 .^ exp_val);
            if eq(exp_val, 0)
                Kd_str = sprintf('%.0f',bas_val);
            else
                Kd_str = sprintf('%.0f\\times10^{%d}',bas_val,exp_val);
            end
            if ind == 5
                Kd_str = "K_{d} = " + Kd_str;
            end
            text(text_x_pos + 0.25,text_y_pos + y_fudge(ind), Kd_str,...
                'HorizontalAlignment','center',...
                'Color',p.Color,...
                'FontWeight','bold',...
                'FontSize',16,...
                'Rotation',angle * 1.75);
        end
    end

    
end
