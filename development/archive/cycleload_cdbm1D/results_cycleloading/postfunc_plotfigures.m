function [outputArg1,outputArg2] = postfunc_plotfigures(x, y, plot_params, save_filename)
    % This function handles plotting of stress-strain curves
    
    % Create the figure
    figure();
    
    % Get number of datasets
    num_datasets = length(x);
    
    % Loop through the datasets and plot them
    for i = 1:num_datasets
        plot(x{i}, y{i}, plot_params.line_styles{i}, 'LineWidth', plot_params.linewidth); 
        hold on;
    end
    
    % Set axis limits
    axis(plot_params.framelimits);
    
    % Set labels, legend, and title
    xlabel(plot_params.xlabel_str, 'FontSize', plot_params.frontsize_titles, 'Interpreter', "latex");
    ylabel(plot_params.ylabel_str, 'FontSize', plot_params.frontsize_titles, 'Interpreter', "latex");
    legend(plot_params.legend_labels, 'Location', "best", 'FontSize', plot_params.frontsize_legend);
    title(plot_params.plot_title, 'FontSize', plot_params.frontsize_titles, 'Interpreter', "latex");
    
    % Set tick label sizes
    ax = gca; % Get the current axes
    ax.XAxis.FontSize = plot_params.frontsize_ticklabels; 
    ax.YAxis.FontSize = plot_params.frontsize_ticklabels;
    
    % Save the figure to a file
    saveas(gcf, save_filename);
end