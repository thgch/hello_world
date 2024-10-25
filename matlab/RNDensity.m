% This is material illustrating the methods from the book
% Financial Modelling  - Theory, Implementation and Practice with Matlab
% source
% Wiley Finance Series, Material for 2nd edition
% ISBN 978-0-470-74489-5
%
% Date: 31.03.2015
%
% Authors:  Joerg Kienitz
%           Daniel Wetterau
%
% Please send comments, suggestions, bugs, code etc. to
% kienitzwetterau_FinModelling@gmx.de
%
% (C) Joerg Kienitz, Daniel Wetterau
% 
% Since this piece of code is distributed via the mathworks file-exchange
% it is covered by the BSD license 
%
% This code is being provided solely for information and general 
% illustrative purposes. The authors will not be responsible for the 
% consequences of reliance upon using the code or for numbers produced 
% from using the code.

function [out1,out2] = RNDensity( K,C )
    
    delta = K(2:end) - K(1:end-1);
    dC = C(2:end) - C(1:end-1);
    dC = dC ./ delta;
    
    delta = K(3:end)-K(2:end-1);
    dC2 = dC(2:end) - dC(1:end-1);
    dC2 = dC2./delta;
    
    out2 = dC2;
    out1 = K(2:end-1);
end

