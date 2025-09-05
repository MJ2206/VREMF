function [CS] = knee(weights)
[sorted_weights, original_indices] = sort(weights, 'descend');
n = length(sorted_weights);
x = 1:n;
y = sorted_weights;
start_point = [1, sorted_weights(1)];
end_point = [n, sorted_weights(end)];
slope = (sorted_weights(end) - sorted_weights(1)) / (n - 1);
intercept = sorted_weights(1) - slope;
distances = abs(slope * x - y + intercept) / sqrt(slope^2 + 1);
[~, knee_index] = max(distances);
knee_point_weight = sorted_weights(knee_index);
selected_indices_sorted = sorted_weights > knee_point_weight;
CS = original_indices(selected_indices_sorted);
end