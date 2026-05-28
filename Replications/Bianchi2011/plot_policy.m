function plot_policy(x_data, y_data, states, target_state, color, label)
    % Find all times this specific shock happened
    idx = (states == target_state);
    
    if sum(idx) > 5
        x_pts = x_data(idx);
        y_pts = y_data(idx);
        
        % Sort them to draw a line
        [x_pts, sort_i] = sort(x_pts);
        y_pts = y_pts(sort_i);
        
        plot(x_pts, y_pts, 'o-', 'Color', color, 'MarkerSize', 4, 'DisplayName', label);
    end
end