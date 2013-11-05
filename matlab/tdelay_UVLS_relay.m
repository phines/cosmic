function [tdelay_UVLS_relay] = tdelay_UVLS_relay(zero_crossing)
% usage: [tdelay_temp_relay] = tdelay_UVLS_relay(zero_crossing)
% obtain the particular time delay based on how much the observed value
% corsses threshold 

% Please adjust the x range and time_delay values
t_temp = @(x)(x>-100 & x<=1e-4).*20 + (x>-200 & x<=-100).*10 + ...
        (x>-300 & x<=-200).*5 + (x>-400 & x<=-300).*2 + (x<=-400).*0;

tdelay_UVLS_relay = t_temp(zero_crossing);