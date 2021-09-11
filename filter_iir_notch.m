function [filtered] = filter_iir_notch(signal, w0, bw, dim)
%Applying IIR Notch Filter (band-stop) at frequency w0, with width bw (3dB)
%   for multiple channel signals dim will specify the dimension to apply the filter
%   for singal channel dim is unused
[a, b] = iirnotch(w0, bw);

if ndims(signal)==1
    filtered=filter(b, a, signal);
elseif ndims(signal)==2
    filtered=filter(b, a, signal, [], 2);
end
end

