% function generategifs(x,y, plot_params)
% 
% % Create a figure for plotting
% figure;
% 
% % Loop through the data points to animate the point along the curve
% for i = 1:length(x)
%     plot(x, y, 'b-', 'LineWidth', 2); % Plot the sine curve
%     hold on;
% 
%     % Plot the moving point
%     plot(x(i), y(i), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% 
%     % Set axis limits
%     axis(plot_params.framelimits);
% 
%     % Add labels and title
%     xlabel(plot_params.xlabel_str,FontSize=plot_params.frontsize_titles,Interpreter="latex");
%     ylabel(plot_params.ylabel_str,FontSize=plot_params.frontsize_titles,Interpreter="latex");
%     title(plot_params.plot_title,FontSize=plot_params.frontsize_titles,Interpreter="latex");
%     ax = gca; % Get the current axes
%     ax.XAxis.FontSize = plot_params.frontsize_ticklabels; % Set the font size of the X-axis
%     ax.YAxis.FontSize = plot_params.frontsize_ticklabels; % Set the font size of the X-axis
% 
%     % Capture the plot as a frame
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind, cm] = rgb2ind(im, 256);
% 
%     % Write to the GIF file
%     if i == 1
%         imwrite(imind, cm, plot_params.filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.02);
%     else
%         imwrite(imind, cm, plot_params.filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.02);
%     end
% 
%     hold off;
% end
% end

function postfunc_generategifs(x, y, plot_params, filename)
    % x and y should be cell arrays containing the datasets
    % Example: x = {x1, x2, x3}, y = {y1, y2, y3}

    % Get the number of datasets
    num_sets = length(x);

    % Find the maximum length of all data sets to loop over
    max_len = max(cellfun(@length, x));

    % Create a figure for plotting
    figure;

    % Plot each dataset's static curve (lines) once
    hold on;
    for j = 1:num_sets
        % Plot the full curve for each dataset
        plot(x{j}, y{j}, plot_params.line_styles{j}, 'LineWidth', 2); 
    end
    hold off;

    % Set axis limits once
    axis(plot_params.framelimits);

    % Add labels and title once
    xlabel(plot_params.xlabel_str, 'FontSize', plot_params.frontsize_titles, 'Interpreter', "latex");
    ylabel(plot_params.ylabel_str, 'FontSize', plot_params.frontsize_titles, 'Interpreter', "latex");
    title(plot_params.plot_title, 'FontSize', plot_params.frontsize_titles, 'Interpreter', "latex");

    % Adjust tick label font sizes
    ax = gca; % Get the current axes
    ax.XAxis.FontSize = plot_params.frontsize_ticklabels;
    ax.YAxis.FontSize = plot_params.frontsize_ticklabels;

    % Loop through the maximum number of data points
    for i = 1:max_len
        % Only plot moving points every 50 steps
        if mod(i, plot_params.ploteveynsteps) == 0 || i == max_len || i == 1 % Plot at every 50th step or at the last step
            % Re-plot the static curves and moving points for each dataset
            hold on;
            for j = 1:num_sets
                % Plot the moving point for the current set only if i is within the dataset's range
                if i <= length(x{j})
                    plot(x{j}(i), y{j}(i), plot_params.marker_styles{j}, 'MarkerSize', 10, 'MarkerFaceColor', plot_params.marker_facecolors{j});
                end
            end
            hold off;

            % Capture the plot as a frame
            frame = getframe(gcf);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);

            % Write to the GIF file
            if i == 1
                imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.02);
            else
                imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.02);
            end

            % Clear only the moving points (keep the static curves)
            if i < max_len
                cla(ax); % Clear only the current frame
                hold on;
                % Re-plot the static curves to persist
                for j = 1:num_sets
                    plot(x{j}, y{j}, plot_params.line_styles{j}, 'LineWidth', 2);
                end
                hold off;
            end
        end
    end
end