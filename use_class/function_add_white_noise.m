function y = function_add_white_noise(matrixX, matrixY)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   White Noise Generation with MATLAB                 %
    %                                                      %
    %   y - function without noise                         %
    %   matrixX - size matrix on X-axis                    %
    %   matrixY - size matrix on Y-axis                    %
    %                                                      %
    %   return:                                            %
    %   y - function with white noise                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = randn(matrixX,matrixY);
    y = x;
end