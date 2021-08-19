% Copyright (c) 2021 Jason Pereira jason.pereira@york.ac.uk
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software isfurnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.


% Choose transmissivity values for the AD channels.
eta1 = 0.4; eta2 = 0.45;

pauliX = [0,1;1,0]; pauliZ = [1,0;0,-1];
count = 0;

for i = 1:100
    % Generate random pure 2-qubit states.
    params = rand(1,8)*pi/2;
    w1 = params(1); r1a = params(2); r1b = params(3); r1c = params(4);
    w2 = params(5); r2a = params(6); r2b = params(7); r2c = params(8);
    rho1 = [sin(w1)^2,0,0,cos(w1)*sin(w1);0,0,0,0;0,0,0,0;...
        cos(w1)*sin(w1),0,0,cos(w1)^2];
    rho2 = [sin(w2)^2,0,0,cos(w2)*sin(w2);0,0,0,0;0,0,0,0;...
        cos(w2)*sin(w2),0,0,cos(w2)^2];
    u1 = kron(eye(2),pauliZ^r1a)*kron(eye(2),pauliX^r1b)*...
        kron(eye(2),pauliZ^r1c);
    u2 = kron(eye(2),pauliZ^r2a)*kron(eye(2),pauliX^r2b)*...
        kron(eye(2),pauliZ^r2c);
    rho1A = u1*rho1*u1'; rho2 = u2*rho2*u2';
    
    % Find the channel outputs.
    out1A = ADoutput(rho1A,eta1,2); out2 = ADoutput(rho2,eta2,2);
    
    % Generate a new pure 2-qubit state.
    params = rand(1,4)*pi/2;
    w1 = params(1); r1a = params(2); r1b = params(3); r1c = params(4);
    rho1 = [sin(w2)^2,0,0,cos(w2)*sin(w2);0,0,0,0;0,0,0,0;...
        cos(w2)*sin(w2),0,0,cos(w2)^2];
    u1 = kron(eye(2),pauliZ^r1a)*kron(eye(2),pauliX^r1b)*...
        kron(eye(2),pauliZ^r1c);
    rho1B = u1*rho1*u1';
    
    % Find the channel output.
    out1B = ADoutput(rho1B,eta1,2);
    
    % Generate a random intermediate state.
    p = rand;
    rho1C = (p*rho1A + (1-p)*rho1B);
    
    % Find the channel output.
    out1C = ADoutput(rho1C,eta1,2);
    
    % Calculate the relative fidelities.
    try
        warning('off','MATLAB:sqrtm:SingularMatrix')
        fidIn = sum(sqrt(abs(eigs(sqrtm(rho1A)*rho2*sqrtm(rho1A)))));
        fidOut = sum(sqrt(abs(eigs(sqrtm(out1A)*out2*sqrtm(out1A)))));
        relFidA = fidOut/fidIn;
        fidIn = sum(sqrt(abs(eigs(sqrtm(rho1B)*rho2*sqrtm(rho1B)))));
        fidOut = sum(sqrt(abs(eigs(sqrtm(out1B)*out2*sqrtm(out1B)))));
        relFidB = fidOut/fidIn;
        fidIn = sum(sqrt(abs(eigs(sqrtm(rho1C)*rho2*sqrtm(rho1C)))));
        fidOut = sum(sqrt(abs(eigs(sqrtm(out1C)*out2*sqrtm(out1C)))));
        relFidC = fidOut/fidIn;
    catch
        warning('on','MATLAB:sqrtm:SingularMatrix')
    end
    
    % Check if superadditivity holds.
    if relFidC < (p*relFidA + (1-p)*relFidB)
        count = count + 1;
    end
end

% Number of violations of superadditivity.
count