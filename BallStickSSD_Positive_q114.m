function [sumRes, S] = BallStickSSD_Positive_q114(x, Avox, bvals, qhat)
% Extract the parameters
S0 = x(1).^2; 

diff = x(2).^2;

f = Siggy(x(3)); %  Must be between [0 1]

theta = x(4);

phi = x(5);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');

S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));

change = Avox';

% Compute the sum of square differences
sumRes = sum((change - S').^2);

