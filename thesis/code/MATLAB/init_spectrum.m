function plot_data = init_spectrum
thesis_data;

if VORTEX
%   plot_data = sort(real([s_freq;-1*s_freq;u_freq]));
   plot_data = freq_c;
else
%   plot_data = sort(real(diag(e_value)));
   plot_data = freq_d;
end

return