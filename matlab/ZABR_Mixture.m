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


clc
clear
close all

alpha = 0.0873; beta = .7;                          % SABR params
rho = -0.48; epsi = 0.47;                           % SABR params

S = 0.0325; T = 10;                                 % asset params


Kl = 0.00; Ku = 0.12; delta = 0.0025;                % params for k range
K = [(Kl:delta:S-delta),S,(S+delta:delta:Ku)]';     % strike range




sigS = @(S)alpha*S.^beta;                           % vol function

gamma = [0,0.5,1,1.5,1.7];                          % ZABR parameter

nG = length(gamma); nK = length(K);                 % no of gammas and Ks
% init arrays and options for later use
nubar = zeros(nK,nG);
nu = zeros(nK,nG);
impvol = zeros(nK,nG);
callMat = zeros(nK,nG); callMat2 = callMat;         % call prices
callnubar = callMat;
options = odeset('RelTol',1e-4,'AbsTol',1e-6);      % options for imp vol calc
fy = zeros(size(K));
lognormalvolv = zeros(nG,nK); lognormalvolv2 = lognormalvolv;

figure;

for j = 1:nG
    
    if abs(gamma(j) - 1.0) < 1e-06 
        adjepsi = epsi;
    else
        adjepsi = epsi*alpha^(1-gamma(j));
    end
    % % approximated integral value
    % y = zeros(size(K));
    % for i = 1:length(K)
    %     y(i) = integral(@(S)1./sigS(S),K(i),S)*alpha^(gamma-2);
    % end

    % analytical integral value
    y = (S^(1-beta)-K.^(1-beta))/(1-beta)*alpha^(gamma(j)-2);

    A = @(y)1+(gamma(j)-2)^2*adjepsi^2*y.^2+2*rho*(gamma(j)-2)*adjepsi*y;
    B = @(y)2*rho*(1-gamma(j))*adjepsi+2*(1-gamma(j))*(gamma(j)-2)*adjepsi^2*y;
    C = (1-gamma(j))^2*adjepsi^2;

    f0 = 0; 
    y0 = 0;
    for i = 1:nK
        if K(i)==S
            fy(i) = 0;
        else
            odefun = @(t,f)0.5*(-B(t).*f+sqrt(B(t).^2.*f.^2-4*A(t).*(C*f.^2-1)))./A(t);
            [tt,yy] = ode45(odefun,[y0,y(i)],f0,options);
            fy(i) = yy(end);
            f0 = fy(i);
            y0 = y(i);
        end
    end

    x = alpha^(1-gamma(j))*fy;
    % implied lognormal volatility
    nubar(:,j) = log(S./K)./x; nubar(K==S,j) = alpha*S^(beta-1);
    % implied normal volatility
    nu(:,j) = (S-K)./x; nu(K==S,j) = S^beta*alpha;
    
    % transfers the implied normal vols in implied black vols
    c1 = nu(:,j)*sqrt(T);
    c2 = (S-K)./c1;
    callBachelier = (S-K).*normcdf(c2)+c1.*normpdf(c2);
    impvol(:,j) = blsimpv(S,K,0,T,callBachelier,2,0,1e-4,{'Call'});

    
    phi = sigS(K)./odefun(y,fy);                % nu(k), eq (9) of paper
   
    % constants used for approximation of normcdf
    a1 = .319381530;     a2 = -.356563782;
    a3 = 1.781477937;    a4 = -1.821255978;
    a5 = 1.330274429;    p = .2316419;
    
    xi = abs(x)/sqrt(T); % should it be xi = abs(fy)/sqrt(T)?;
    Px = 2*(1-xi.*(normcdf(-xi)./normpdf(xi)));
    
    % approx of normcdf(-xi) from Abramowitz, Stegun
    tx = 1./(1-p*xi);
    Px2 = 1-normpdf(-xi).*tx.*(a1+tx.*(a2+tx.*(a3+tx.*(a4+tx*a5))));
    Px2 = 2*(1-xi.*Px2./normpdf(xi));
    
    theta2 = phi.^2.*Px2; %  theta2 = phi.^2.*Px2;  % exact or approx
    
    callMat(1,j) = callBachelier(1); 
    callMat2(1,j) = callBachelier(1);
    callMat(end,j) = callBachelier(end); 
    callMat2(end,j) = callBachelier(end);
    
    %zj = .5*T*theta2/delta^2;
    lhs = eye(length(theta2(2:end-1)))-0.5*T*(diag(theta2(3:end-1),-1)-2*diag(theta2(2:end-1))+diag(theta2(2:end-2),1))/delta.^2;
    rhs = max(S-K(2:end-1),0);
    BoundaryCond = zeros(size(rhs));
    BoundaryCond(end)= (theta2(end-1)/delta.^2) * callBachelier(end);
    BoundaryCond(1)= (theta2(2)/delta.^2) * callBachelier(1);
    rhs = rhs + 0.5 * T * BoundaryCond;
    callMat(2:end-1,j) = lhs\rhs;
	callMat(1,j) = callBachelier(1); 
    
    % without adjustments
    theta2 = phi.^2;
    lhs = eye(length(theta2(2:end-1)))-0.5*T*(diag(theta2(3:end-1),-1)-2*diag(theta2(2:end-1))+diag(theta2(2:end-2),1))/delta.^2;
    rhs = max(S-K(2:end-1),0);
    BoundaryCond = zeros(size(rhs));
    BoundaryCond(end)= (theta2(end-1)/delta.^2) * callBachelier(end);
    BoundaryCond(1)= (theta2(2)/delta.^2) * callBachelier(1);
    rhs = rhs + 0.5 * T * BoundaryCond;
    callMat2(2:end-1,j) = lhs\rhs;
	callMat2(1,j) = callBachelier(1); 
	callMat2(end,j) = callBachelier(end);
    
    callnubar(:,j) = blsprice(S,K,0,T,nubar(:,j));
    %figure
    subplot(nG,3,3*(j-1)+1);
    plot(K,callMat(:,j),'-.',K,callnubar(:,j),'.',K,callMat2(:,j),'-','LineWidth',2); title(['Call Option Prices adjusted and unadjusted \beta_v=',num2str(gamma(j))]);
    legend('adjusted','standard','raw')
    
    [a,b] = RNDensity(K,callMat(:,j)); [aa,bb] = RNDensity(K,callnubar(:,j)); [aaa,bbb] = RNDensity(K,callMat2(:,j));
    
    %figure;
    subplot(nG,3,3*(j-1)+2);
    plot(a,b,'-.',aa,bb,'.',aaa,bbb,'-','LineWidth',2); hold off; title(['Densities adjusted and unadjusted \beta_v=',num2str(gamma(j))]);
    legend('adjusted','standard','raw')

    lognormalvolv(j,:) = blsimpv(S,K,0,T,callMat(:,j),2,0,1e-4,{'Call'}); lognormalvolv2(j,:) = blsimpv(S,K,0,T,callMat2(:,j),2,0,1e-4,{'Call'});
    lnvol = blsimpv(S,K,0,T,callnubar(:,j),2,0,1e-4,{'Call'});
    %for k = 1:length(K)
    %    lognormalvolv(j,k) =  blsimpv(S,K,0,T,callMat(:,j),2,0,1e-4,{'Call'});%ImpliedVolLN(S,K(k),0,0,T,callMat(k,j),1,1e-4,0,2,50);
    %end 
    %figure;
    subplot(nG,3,3*(j-1)+3);
    plot( K(2:end-1),lognormalvolv(j,2:end-1),'-.', K(2:end-1),lnvol(2:end-1),'.', K(2:end-1),lognormalvolv2(j,2:end-1),'-','LineWidth',2); title(['Lognormal Implied Volatilities \beta_v=',num2str(gamma(j))]);
    legend('adjusted','standard','raw')


end

% figure 3 of Andreasen
figure
plot(K,impvol);
figure;
plot(K,nubar)
grid on
legend('gam = 0','gam = 0.5','gam = 1','gam = 1.5','gam = 1.7','Location','NorthEastOutside')


% mixture of ZABR
% left from S we use mixture mix11 / mix12
mix11 = 0.9; mix12 = 0.1; 
% right from S we use mixture mix21 / mix22 
mix21 = 0.5; mix22 = 0.5;

% Calculate densities and mix 
[a,b] = RNDensity(K(3:end),callMat(3:end,1)); 
[aa,bb] = RNDensity(K(3:end),callMat(3:end,end));
aaa = a; bbb=b;
bbb(a<S) = mix11 * b(a<S) + mix12 * bb(a<S);
bbb(a>=S) = mix21 * b(a>=S) + mix22 * bb(a>=S);

figure;
hold on; plot(a,b,'b'); plot(aa,bb,'r'); plot(aaa,bbb,'g'); hold off;
legend('gam = 0','gam = 1.7','gam = mix','Location','NorthEastOutside')
title('Mixture of ZABRs')