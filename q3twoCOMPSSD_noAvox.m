function S = q3twoCOMPSSD_noAvox(x, bvals, qhat)
% Extract the parameters
S0 = x(1); 

diff = x(2);

f = Siggy(x(3)); %  Must be between [0 1]

theta = x(4);

phi = x(5);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

% fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');

fibdotgrad_q3 = (sum(qhat.*repmat(fibdir, [length(qhat) 1])')).^2;

Si = exp(-bvals*diff.*fibdotgrad_q3);

Se = exp(-bvals.*diff);

S = S0.*(f.*Si + (1-f).*Se);

