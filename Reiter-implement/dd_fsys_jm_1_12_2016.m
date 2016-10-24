function [F, FLHS, FRHS] = dd_fsys_jm_12015(X, Xmin1, eta, eps) 

global VDind  VNDind Vind XIind distind Ckind Nkind ...
    Zkind  Lambda_t_t1ind Brkind Bind  Qind Lambdaind  ...
    Bgind Wind Tkind Iind Aind Yind Dgind Mind ...
    Piind Dkind Rind Lind Kind Gind deltaind numstates numeps ... 
    etac xik Delta eta1 eta2 rho_zk rho_a sigma_zk ... 
    alpha Tkbar gamma_t bybar ibar phi_i gybar gamma_g znum bnum ...
    dsnum bdensenum betab deltak betak nmaxit ntol maxvfit vftol ...
    bmaxit gsstol gssphi psi Ikind defls b0ind pr_mat_z pr_mat_ds ...
    bdense0ind etaVDind etaVNDind etaBRind etaBind etaKind xibar ...
    b0 bdense0 z0 epskzind epsiind epsaind indmat denseindmat; 


%%%%%%%%%%%%%%%%%%%%%%%
%CREATES THE SYSTEM OF NONLINEAR EQUATIONS TO BE PERTURBED
%
%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%

F       = zeros(numstates, 1); 
FLHS    = zeros(numstates, 1); 
FRHS    = zeros(numstates, 1); 

%%
%READ IN DATA: NOTE THAT ANYTHING INDEXED "MIN1" is "t" IN THE MODEL
%WRITEUP
CKmin1  = Xmin1(Ckind);
NKmin1  = Xmin1(Nkind);
ZKmin1  = Xmin1(Zkind);
%Lambda_t_t1min1     = Xmin1(Lambda_t_t1ind);
Brkmin1 = Xmin1(Brkind);    
Bmin1   = Xmin1(Bind);
Qmin1val   = Xmin1(Qind);
Lambdamin1  = Xmin1(Lambdaind);
Bgmin1  = Xmin1(Bgind);
Wmin1val   = Xmin1(Wind);
Tkmin1  = Xmin1(Tkind);
Imin1  = Xmin1(Iind);
Amin1   = Xmin1(Aind);
Ymin1   = Xmin1(Yind);
Dgmin1  = Xmin1(Dgind);
Mmin1   = Xmin1(Mind);
Pimin1    = Xmin1(Piind);
Dkmin1  = Xmin1(Dkind);
Lmin1   = Xmin1(Lind);
Kmin1   = Xmin1(Kind);
Gmin1   = Xmin1(Gind);
deltamin1   = Xmin1(deltaind);
Rmin1       = Xmin1(Rind);
Ikmin1      = Xmin1(Ikind);

Ck  = X(Ckind);
NK  = X(Nkind);
ZK  = X(Zkind);
%Lambda_t_t1     = X(Lambda_t_t1ind);
Brk = X(Brkind);
B   = X(Bind);
Q   = X(Qind);
Lambda  = X(Lambdaind);
Bg  = X(Bgind);
W   = X(Wind);
Tk  = X(Tkind);
I  = X(Iind);
A   = X(Aind);
Y   = X(Yind);
Dg  = X(Dgind);
M   = X(Mind);
Pi    = X(Piind);
Dk  = X(Dkind);
L   = X(Lind);
K   = X(Kind);
G   = X(Gind);
delta   = X(deltaind);
R   = X(Rind);
IK   = X(Ikind);

%%
%EVALUATING THE SYSTEM
FLHS(Ckind)  = CKmin1;
FLHS(Nkind)  = NKmin1;
FLHS(Zkind)  = log(ZK);
%FLHS(Lambda_t_t1ind)   = Lambda_t_t1min1;
FLHS(Brkind) = Brkmin1;
% FLHS(Bind)  = - Bmin1 - Qmin1val * Brkmin1 * etac; 
FLHS(Bind)  =  Bmin1; %Bmin1 and B track TOTAL debt, while Bg and Br track govt and risky debt, respectively
FLHS(Qind)   = Qmin1val;
FLHS(Lambdaind)  = Lambdamin1;
FLHS(Bgind)  = Bg;
FLHS(Wind)   = Wmin1val;
FLHS(Tkind)  = Tkmin1;
FLHS(Iind)  = Imin1;
FLHS(Aind)   = log(A);
FLHS(Yind)   = Ymin1;
FLHS(Dgind)  = Dgmin1;
FLHS(Mind)   = Mmin1;
FLHS(Piind)    = Pi;
FLHS(Dkind)   = Dkmin1;
FLHS(Lind)   = Lmin1;
FLHS(Kind)   = 1 + xik*(Ikmin1/Kmin1 - deltak);
FLHS(Gind)   = log(Gmin1);
FLHS(deltaind)   = deltamin1;
FLHS(Rind)   = Rmin1;
FLHS(Ikind)   = Ikmin1;


%%
SSconstants    = importdata('SSconstants.txt');
NKss           = SSconstants(1);
Brkss          = SSconstants(2);
Qss         = SSconstants(3);
Lambdass    = SSconstants(4);
Tkss    = SSconstants(5);
Gss     = SSconstants(6);
Bgss    = SSconstants(7);
Wss     = SSconstants(8);
Iss     = SSconstants(9);
Ass     = SSconstants(10);
Yss     = SSconstants(11);
Kss     = SSconstants(12);
Lss     = SSconstants(13);
Mss     = SSconstants(14);
Dgss    = SSconstants(15);
Piss    = SSconstants(16);
Dkss    = SSconstants(17);
Rss     = SSconstants(18);
deltass = SSconstants(19);
Ikss    = SSconstants(20);
Bss     = SSconstants(21);
Ckss    = SSconstants(22);
Zkss    = SSconstants(23);
Lambda_t_tss = SSconstants(24);

CB = (Yss - Ckss - Gss - deltak * Kss + Delta* Brkss)/etac;
NB =  (Lss - Zkss * NKss ) / etac ;
Br = Brkss/etac;

%Fraction of debt defaulted on
Ckmin1       = Ymin1 - etac * CB - Gmin1 - K + Kmin1 - deltak * Kmin1 - xik * ( ( K - Kmin1) / Kmin1)^2 * Kmin1 / 2 + Delta  * Brkmin1; 
FRHS(Ckind)  = Ckmin1;
FRHS(Nkind)  = ( (Wmin1val * ZKmin1) / (eta1 * Ckmin1) )^(1 / eta2);
FRHS(Zkind)  = rho_zk * log(ZKmin1) + eps(epskzind);
Lambda_t_t1min1 =  betak * (Ckmin1 / Ck);
FRHS(Brkind) = etac*Br;
FRHS(Bind) = Bgmin1 + Qmin1val*Brkmin1;
FRHS(Qind)     = Lambda_t_t1min1 * (1 - delta - Delta) / (1 + Pi) + eta(etaBRind);
FRHS(Lambdaind)  = betak * Lambda * (1 + I)/(1 + Pi) + eta(etaBind);
FRHS(Bgind)  = -G + ((1 + I)/(1 + Pi)) * Bgmin1 + Tk;
FRHS(Wind)   = Mmin1 * (1 - alpha) * Amin1 * Kmin1^alpha * Lmin1^(- alpha);
FRHS(Tkind)  = Tkbar + gamma_t * log( ( Bgmin1 / Ymin1) / bybar); 
FRHS(Iind)   = ibar + phi_i * log(1 + Pimin1) + eps(epsiind);
FRHS(Aind)   = rho_a* log(Amin1) + eps(epsaind);
FRHS(Yind)   = Amin1 * Kmin1^(alpha) * Lmin1^(1 - alpha); 
FRHS(Dgind)  = Ymin1 - Mmin1 *(Amin1 * Kmin1^(alpha) * Lmin1^(1 - alpha) ); 
FRHS(Mind)   = (1 / Amin1)*(((Rmin1 + deltak)/alpha)^(alpha)* (Wmin1val/(1 - alpha))^(1-alpha)); 
FRHS(Piind)  = M / Mmin1-1; 
FRHS(Dkind)  = (1 + Rmin1)*Kmin1 - K - ((xik /2) * ( (K - Kmin1)/Kmin1)^2) * Kmin1; 
FRHS(Rind)    = Mmin1 * alpha * Amin1 * Kmin1^(alpha - 1) * Lmin1^(1 - alpha) - deltak; 
FRHS(Lind)   = ZKmin1 * NKmin1 + etac * NB; 
FRHS(Kind)   = Lambda_t_t1min1 * ( R + 1 + xik * ( IK / K - deltak) * (IK / K - deltak + 1) - (xik/2) * (IK / K - deltak)^2) + eta(etaKind); 
FRHS(deltaind)   = DeltaB;
FRHS(Gind)   = log(Ymin1) + log(gybar) - gamma_g * log((Bgmin1 / Ymin1)/bybar);
FRHS(Ikind)   = K - (1 - deltak) * Kmin1; 
F = FLHS - FRHS; 

end