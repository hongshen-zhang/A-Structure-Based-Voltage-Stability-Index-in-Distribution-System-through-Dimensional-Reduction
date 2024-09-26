%% Easily calculate the error for real and predict list
% date : 20210824
% author : hongshen zhang
function [average_error,median_error] = average_median_error(real,predict)
average_error = mean(100 * abs(real - predict) ./ real);
median_error = median(100 * abs(real - predict) ./ real);
end