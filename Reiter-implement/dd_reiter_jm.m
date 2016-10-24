%%%%%%%%%%%%%%%%%%%%%%%
%MASTER FILE TO PERTURB DD MODEL AROUND SS
%
%
%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%

clear 

%%
%GLOBALS
global VDind  VNDind Vind XIind distind Ckind Nkind ...
    Zkind  Lambda_t_t1ind Brkind Bind  Qind Lambdaind  ...
    Bgind Wind Tkind Iind Aind Yind Dgind Mind ...
    Piind Dkind Rind Lind Kind Gind deltaind numstates numeps ... 
    etac delta xik Delta eta1 eta2 rho_zk rho_a sigma_zk ... 
    alpha Tkbar gamma_t bybar ibar phi_i gybar gamma_g znum bnum ...
    dsnum bdensenum betab deltak betak nmaxit ntol maxvfit vftol ...
    bmaxit gsstol gssphi psi Ikind defls b0ind pr_mat_z pr_mat_ds ...
    bdense0ind etaVDind etaVNDind etaBRind etaBind etaKind xibar ...
    b0 bdense0 z0 epskzind epsiind epsaind; 

%%
%DIMENSION DECLARATIONS
constantvec = importdata('constantvec.txt'); 
znum        = constantvec(1); 
bnum        = constantvec(2); 
dsnum       = constantvec(3); 
bdensenum   = constantvec(4); 
betab       = constantvec(5); 
eta1        = constantvec(6);
eta2        = constantvec(7); 
etac         = constantvec(8); 
Delta       = constantvec(9); 
alpha       = constantvec(10); 
deltak      = constantvec(11);
betak       = constantvec(12); 
nmaxit      = constantvec(13); 
ntol        = constantvec(14); 
maxvfit  = constantvec(15); 
vftol    = constantvec(16); 
bmaxit   = constantvec(17); 
gsstol   = constantvec(18); 
gssphi   = constantvec(19); 
psi      = constantvec(20); 
b0ind    = constantvec(21);
bdense0ind = constantvec(22); 
xibar   = constantvec(23); 
xik     = constantvec(24);
rho_zk  = constantvec(25);
rho_a   = constantvec(26);
ibar    = constantvec(27);
phi_i   = constantvec(28);
bybar   = constantvec(29);
gybar   = constantvec(30);
Tkbar   = constantvec(31);
gamma_t = constantvec(32);
gamma_g = constantvec(33);


b0 = importdata('b0.txt'); 
bdense0 = importdata('bdense0.txt'); 
z0 = importdata('z0.txt');

Vss     = importdata('V.txt');
VNDss     = importdata('VND.txt');
VDss     = importdata('VD.txt');

pr_mat_ds    = importdata('pr_mat_ds.txt');
pr_mat_z     = importdata('pr_mat_z.txt');

distss  = importdata('dist.txt');

XIss = importdata('xistarmat.txt');

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
defls   = SSconstants(25);

%%
%X VECTOR 
count = 0; 
VDind   = 1:znum;
count   = count+length(VDind); 
VNDind  = (count + 1) : (znum*bnum * dsnum + count); 
count   = count +length(VNDind); 
Vind    = (count + 1) : (znum*bnum * dsnum + count); 
count   = count + length(Vind);
%XIind   = (count + 1 ) : (znum*bnum*dsnum + count);
%count   = count + length(XIind);
distind   = (count + 1 ) : (znum*bdensenum*dsnum + count); 
count   = count + length(distind);
Ckind   = count+1;
count   = count + 1;
Nkind   = count+1;
count   = count + 1;
Zkind   = count+1;
count   = count + 1;
Lambda_t_t1ind   = count+1;
count   = count + 1;
Brkind   = count+1;
count   = count + 1;
Bind   = count+1;
count   = count + 1;
Qind   = count+1;
count   = count + 1;
Lambdaind   = count+1;
count   = count + 1;
Bgind   = count+1;
count   = count + 1;
Wind   = count+1;
count   = count + 1;
Tkind   = count+1;
count   = count + 1;
Iind   = count+1;
count   = count + 1;
Aind   = count+1;
count   = count + 1;
Yind   = count+1;
count   = count + 1;
Dgind   = count+1;
count   = count + 1;
Mind   = count+1;
count   = count + 1;
Piind   = count+1;
count   = count + 1;
Dkind   = count+1;
count   = count + 1;
Lind   = count+1;
count   = count + 1;
Kind   = count+1;
count   = count + 1;
Gind   = count+1;
count   = count + 1;
deltaind   = count+1;
count   = count + 1;
Rind    = count + 1;
count   = count + 1; 
Ikind   = count + 1; 
numstates = count + 1; 

%EPSILONS
epskzind = 1; 
epsiind  = 2; 
epsaind  = 3; 

epsss = zeros(3,1);

%EXPECTATIONAL ERRORS
counteta = 0; 
etaVDind = 1:znum;
counteta = counteta + length(etaVDind); 
etaVNDind  = (counteta + 1) : (znum*bnum * dsnum + counteta);
counteta = counteta + length(etaVNDind); 
etaBRind = counteta + 1;
counteta = counteta + 1; 
etaBind  = counteta + 1; 
counteta = etaBind + 1; 
etaKind  = counteta + 1; 

numetas = counteta+1;

etass = zeros(numetas,1);

%%
%INPUT STEADY STATE VECTOR
Xss = zeros(numstates,1);
Xss(VDind) = VDss; 
Xss(VNDind) = VNDss; 
Xss(Vind)   = Vss; 
%Xss(XIind)   = XIss; 
Xss(distind)   = distss; 
Xss(Ckind)   = Ckss; 
Xss(Nkind)   = NKss;
Xss(Zkind)   = Zkss;
Xss(Lambda_t_t1ind)   = Lambda_t_tss;
Xss(Brkind)   = Brkss;
Xss(Bind)       = Bss;
Xss(Qind)       = Qss;
Xss(Lambdaind)  = Lambdass;  
Xss(Bgind)      = Bgss;
Xss(Wind)       = Wss;
Xss(Tkind)      = Tkss;
Xss(Iind)       = Iss;
Xss(Aind)   = Ass;
Xss(Yind)   = Yss;
Xss(Dgind)  = Dgss; 
Xss(Mind)   = Mss;
Xss(Piind)      = Piss; 
Xss(Dkind)      = Dkss;
Xss(Lind)        = Lss;
Xss(Kind)       = Kss;
Xss(Gind)       = Gss; 
Xss(deltaind)   = deltass;
Xss(Rind)       = Rss;
Xss(Ikind)      = Ikss; 


%CALL SOLUTION
disp('COMPUTING!!!');
[Fss,FLHss,FRHss] = dd_fsys_jm(Xss,Xss,etass,epsss); 



