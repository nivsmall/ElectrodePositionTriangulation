function [V_conjugates] = continuousStream2triggered(cont_stream, trig_indices, exciting_channels, num_trials)
% converts 2D continuous stream to 3D, folding data @ triggers (in time)
L=length(trig_indices);
W = floor((trig_indices(L)-trig_indices(1))/L);	% Mat width - #samples in excite period
H = length(cont_stream(:, 1));                  % Mat height - #channels
D = exciting_channels;                          % Mat depth  - #exciting channels
contM = zeros(H, W, D, num_trials);
V_conjugates = zeros(H, D);
for experiment = 1:num_trials
    for i = 1:D
        %disp(trig_indices(i)+W-11)
        contM(:, :, i, experiment) = cont_stream(:, trig_indices(i+D*(experiment-1)): trig_indices(i+D*(experiment-1))+W-1); % cont_stream(:, (i-1)*W+1: i*W);
    end
end
M = min(contM, [],2);
M = reshape(M, [H,D,num_trials]);

for i = 1:H
    for j = 1:D
        [Vhat,Vpci] = mle(reshape(M(i,j,:), [num_trials,1]))
        %disp(reshape(M(i,j,:), [num_trials,1]))
        %disp(V);
        V_conjugates(i,j) = Vhat(1);    % maximum likelihood estimate over trials.
    end
end
end

%W = (trig_indices(2)-trig_indices(1))*exciting_channels;    % Mat width - excitation signal period (in the sense of #samples)
%H = length(cont_stream(:, 1));                              % Mat height - #channels
%D = length(trig_indices)/exciting_channels;                 % Mat depth  - #channels*#trials
%T = num_trials
%assert(floor(D)==D, 'set up of triggers and channel data is incorrect!');