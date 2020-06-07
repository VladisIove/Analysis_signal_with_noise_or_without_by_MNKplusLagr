function y = function_add_grey_noise(m, n)
    % function: y = function_add_blue_noise(m, n)
    % m - number of matrix rows
    % n - number of matrix columns
    % y - matrix with grey noise samples 
    % The function generates a matrix of grey noise samples
    % (columnwise). In terms grey noise it is sum of two noise:
    % broun noise and violet noise
    y = function_add_broun_noise(m,n) + function_add_violet_noise(m,n);
end