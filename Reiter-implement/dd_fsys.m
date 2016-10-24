function [F, FLHS, FRHS] = dd_fsys(X, Xmin1, eta, eps) 

global VDind  VNDind Vind XIind distind Ckind NKind ...
    Zkind  Lambda_t_t1ind Brkind Bind  Qind Lambdaind  ...
    Bgind Wind Tkind Iind Aind Yind Dgind Mind ...
    Piind Dkind Rind Lind Kind Gind deltaind numstates numeps ... 
    etac delta xik Delta eta1 eta2 rho_zk rho_a sigma_zk ... 
    alpha Tkbar gamma_t bybar ibar phi_i gybar gamma_g znum bnum ...
    dsnum bdensenum betab deltak betak nmaxit ntol maxvfit vftol ...
    bmaxit gsstol gssphi psi Ikind defls b0ind pr_mat_z pr_mat_ds ...
    bdense0ind etaVDind etaVNDind etaBRind etaBind etaKind xibar ...
    b0 bdense0 z0 epskzind epsiind epsaind; 


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
%READ IN DATA
VDmin1  = Xmin1(VDind);
VNDmin1 = Xmin1(VNDind); 
Vmin1    = Xmin1(Vind);
%XImin1  = Xmin1(XIind);
distmin1  = Xmin1(distind);
CKmin1  = Xmin1(Ckind);
NKmin1  = Xmin1(NKind);
ZKmin1  = Xmin1(Zkind);
Lambda_t_t1min1     = Xmin1(Lambda_t_t1ind);
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

VD  = X(VDind);
VND = X(VNDind); 
V    = X(Vind);
%XI  = X(XIind);
dist  = X(distind);
CK  = X(Ckind);
NK  = X(NKind);
ZK  = X(Zkind);
Lambda_t_t1     = X(Lambda_t_t1ind);
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
FLHS(VDind)     = VDmin1;
FLHS(VNDind)    = VNDmin1;
FLHS(Vind)      = Vmin1;
%FLHS(XIind)  = XImin1;
FLHS(distind)  = distmin1;
FLHS(Ckind)  = CKmin1;
FLHS(NKind)  = NKmin1;
FLHS(Zkind)  = ZKmin1;
FLHS(Lambda_t_t1ind)   = Lambda_t_t1min1;
FLHS(Brkind) = Brkmin1;
FLHS(Bind)  = - Bmin1 - Qmin1val * Brkmin1;
FLHS(Qind)   = Qmin1val;
FLHS(Lambdaind)  = Lambdamin1;
FLHS(Bgind)  = Bgmin1;
FLHS(Wind)   = Wmin1val;
FLHS(Tkind)  = Tkmin1;
FLHS(Iind)  = Imin1;
FLHS(Aind)   = log(Amin1);
FLHS(Yind)   = Ymin1;
FLHS(Dgind)  = Dgmin1;
FLHS(Mind)   = Mmin1;
FLHS(Piind)    = Pimin1;
FLHS(Dkind)   = Dkmin1;
FLHS(Lind)   = Lmin1;
FLHS(Kind)   = 1 + xik*(Ikmin1/Kmin1 - delta);
FLHS(Gind)   = log(Gmin1);
FLHS(deltaind)   = deltamin1;
FLHS(Rind)   = Rmin1;
FLHS(Ikind)   = Ikmin1;




%%
%CONSTRUCTING VALUE FUNCTIONSs
V = reshape(V, [dsnum bnum znum]);
V2 = 0*V;

for dscount = 1:dsnum 
    for zcount = 1:znum 
        V2(dscount,:, zcount) = splinefunc(b0, V(dscount,:,zcount), 1e30, 1e30);
    end
end

%DEFAULT VALUE FUNCTIONS
dscount = 2; 
for zcount = 1:znum 
    zval = z0(zcount);
    cval = Wmin1val*zval*defls;
    FRHS(VDind(zcount)) = log(cval) - eta1 * (defls)^(1 + eta2)/(1 + eta2);    

    for zprimecount = 1:znum 
        for dsprimecount = 1:2
            FRHS(VDind(zcount)) = FRHS(VDind(zcount)) + betab * pr_mat_z(zcount, zprimecount) * ...
                            pr_mat_ds(dscount, dsprimecount) * ...
                            splintfunc(b0,V(dsprimecount, :, zprimecount), V2(dsprimecount, :, zprimecount), 0); 
        end        
    end    
    FRHS(VDind(zcount)) = FRHS(VDind(zcount)) + eta(etaVDind(zcount));
end

%NO DEFAULT VALUE FUNCTIONS
tempcount = 0;

for zcount = 1:znum 
    for bmin1count = 1:bnum 
        for dscount = 1:dsnum
            tempcount = tempcount + 1;
            
            %if in defaulted state value same as defaulted state value
            if dscount == 2
                FRHS(VNDind(tempcount)) = FRHS(VDind(zcount)) - eta(etaVDind(zcount)) + eta(etaVNDind(tempcount));
            end
            
            %if not in defaulted state
            if dscount == 1 
            
                zval = z0(zcount);
                bmin1val = b0(bmin1count);
                               
                %GSS 
                a = b0(1); 
                c = b0(bnum);
                b = gssphi*a + (1 - gssphi)*c;
                           
                for goldenit = 1:bmaxit
                    
                    if abs(a-c) < gsstol
                        break;
                    end
                    
                    
                    if goldenit == 1
                        d = a;
                    elseif goldenit == 2
                        d = b;
                    elseif goldenit == 3
                        d = c;
                    else
                        if abs(b-c) > abs(b-a)
                            d = b*gssphi + c*(1 - gssphi);
                            bcflag = 1;
                            abflag = 0;
                        else 
                            d = a*gssphi + b*(1 - gssphi);
                            bcflag = 0;
                            abflag = 1;
                        end

                    end
                    
                    
                    %evaluate fd
                    bval = d;
                    
                    %last period borrowing and current borrowing are determined now
                    nub = 10.0;
                    
                    %asymptote of the N equation
                    nlb = (Qmin1val*bval - bmin1val)./(Wmin1val*zval);
                    %nlb = max([nlb,zeros(bnum,1)],[],2)+nlbcons;
                    nlb = max([nlb 0]);

                    %will be negative
                    na = nlb;
                    nfa = eta1*(na^eta2) - (Wmin1val*zval/(Wmin1val*zval*na+bmin1val-Qmin1val*bval));

                    %will be positive
                    nb = nub;
                    nfb = eta1*(nb^eta2) - (Wmin1val*zval/(Wmin1val*zval*nb+bmin1val-Qmin1val*bval));

                    %actually perform bisection
                    for nct=1:nmaxit;

                        %test for convergence
                        if (max(abs(nb-na))<ntol); break; end;

                        %if not yet converged, choose midpoint candidate value
                        nx = (na+nb)/2;
                        fx = eta1*(nx^eta2) - (Wmin1val*zval/(Wmin1val*zval*nx+bmin1val-Qmin1val*bval));

                        %classify errors, noting that fx is strictly increasing
                        if fx >= 0 
                            nb = nx;
                        else
                            na = nx;
                        end
                        
                    end
                    
                    nval = nx;
                    nvalmat(zcount,bmin1count) = nval;                    
                    
                    cval = bmin1val + Wmin1val * zval * nval - Qmin1val*bval;
                    cvalmat(zcount,bmin1count) = cval;
                    
                    fd = log(cval) - eta1 * nval^(1 + eta2)/(1 + eta2);
                    
                    for zprimecount = 1:znum 
                        fd = fd + betab * ...
                            pr_mat_z(zcount, zprimecount) * ...
                             splintfunc(b0,V(1, :, zprimecount), V2(1, :, zprimecount), bval); 
                    end
                    
                    if goldenit == 1
                        fa = fd;
                    elseif goldenit == 2
                        fb = fd;
                    elseif goldenit == 3
                        fc = fd;
                    else
                                                              
                        if bcflag == 1 && fd > fb 
                            a = b;
                            b = d;
                            fa =fb;
                            fb = fd;

                        elseif bcflag == 1 && fd <= fb
                            c = d;
                            fc = fd;

                        elseif abflag == 1 && fd > fb
                            b = d;
                            fb = fd;

                            c = b;
                            fc = fb;

                        elseif abflag == 1 && fd <= fb 
                            a = d;
                            fa = fd;
                        end
                    end
                    
                                   
                end
                
                %borrower policy function
                bvalmat(zcount,bmin1count) = d;
                
                FRHS(VNDind(tempcount)) = fd + eta(etaVNDind(tempcount));
            end
            
            if dscount == 1 
               
                %default "policy function"
                xivalmat(zcount,bmin1count) = VDmin1(zcount) - VNDmin1(tempcount);
                %zcount, bmin1count,dscount);
                
                xival = xivalmat(zcount,bmin1count);
                
                if xival < 0 
                    FRHS(Vind(tempcount)) = VNDmin1(tempcount);
                elseif xival > xibar 
                    FRHS(Vind(tempcount)) = VDmin1(zcount);
                else
                    FRHS(Vind(tempcount)) = -xival^2 / (2*xibar) + (xival / xibar) * VDmin1(zcount) ...
                        + (1-xival/xibar) * VNDmin1(tempcount);                  
                end                               
            end
        end
    end
end


%%
%DISTRIBUTION PUSH

%policies defined on a dense grid
xival2mat = zeros(znum,bdensenum);
bval2mat    = zeros(znum,bdensenum);

%linear interpolation
for zcount = 1:znum
    for bcount = 1:bdensenum 
        
        %current debt position
        bvalt = bdense0(bcount);
        
        %getting interval and interpolation weights on "sparse" matrix
        bind = sum(bval >= b0);
        
        if bind < 1 
            bind =1;
            bweight =0;
        elseif bind >= bdensenum 
            bind = bdensenum-1;
            bweight = 1;
        else 
            bweight = (bval - bdense0(bind))/(bdense0(bind+1) - bdense0(bind));
        end
        
        %borrowing policy interpolated
        bval2mat(zcount,bcount) = (1- bweight) * bvalmat(zcount,bind) + bweight * bvalmat(zcount,bind+1);
        
        %default policy interpolated
        xival2mat(zcount,bcount) = (1 - bweight)*xivalmat(zcount,bind) + bweight * xivalmat(zcount, bind+1);
       
    end
end
    
%pushing distribution
tempcount = 0;
for zcount = 1:znum 
    for bmin1count = 1:bdensenum 
        for dscount = 1:dsnum
            
            tempcount = tempcount +1;
                        
            %defaulted mass
            if dscount == 2
                for zprimecount = 1:znum
                    for dsprimecount = 1:dsnum
                       temprimecount = (zprimecount-1)*dsnum *bdensenum + (bdense0ind - 1)*dsnum + dsprimecount;
                       FRHS(distind(temprimecount )) = FRHS(distind(temprimecount)) + pr_mat_z(zcount,zprimecount) * ...
                            pr_mat_ds(dscount,dsprimecount) * distmin1(tempcount);
                    end
                end
                
            else
                
                xival = xival2mat(zcount,bmin1count);
                bval = bval2mat(zcount,bcount);

                
                bind = sum(bval >= bdense0);
                
                if bind < 1 
                    bind =1;
                    bweight =0;
                elseif bind >= bdensenum 
                    bind = bdensenum-1;
                    bweight = 1;
                else 
                    bweight = (bval - bdense0(bind))/(bdense0(bind+1) - bdense0(bind));
                end
                
                for zprimecount = 1:znum
                    
                    temprimecount = (zprimecount-1)*dsnum *bdensenum + (bind - 1)*dsnum + 1;
                    
                    FRHS(distind(temprimecount)) = FRHS(distind(temprimecount)) + pr_mat_z(zcount,zprimecount) * ...
                        distmin1(tempcount) *(1 - bweight);
                                              
                    
                    %temprimecount = (zprimecount-1)*dsnum *bdensenum + (bind + 1)*dsnum + 1;
                    temprimecount = temprimecount + 1;  
                    
                    FRHS(distind(temprimecount)) = FRHS(distind(temprimecount)) + pr_mat_z(zcount,zprimecount) * ...
                       distmin1(tempcount) *bweight;
                            
                end
            end            
        end
    end
end

%%
%COMPUTE SOME AGGREGATES
%DENSE POLICY FUNCTIONS FOR LABOR SUPPLY AND CONSUMPTION

for zcount = 1:znum    
    
    zval = z0(zcount);
    
    for bmin1count = 1:bdensenum
    
        bval1 = bdense0(bmin1count);
        bval2 = bval2mat(zcount, bmin1count);
            
        xival = xival2mat(zcount, bmin1count);
                
        if xival < 0 
            pdefault = 0;
        elseif xival > xibar 
            pdefault = 1;
        else
            pdefault = xival / xibar;
        end

        %calculate labor supply with golden section search

        %set out bounds of labor supply search
        nlb = (Qmin1val*bval2 - bval1) / (Wmin1val*zval);
        nlb = max([nlb 0]);

        %initialize and do function evals
        na = nlb + 1e-4;
        nb = nlb + 10.0;

        %will be negative
        fa = eta1*(na^eta2) - (Wmin1val*zval/(Wmin1val*zval*na+bval1-Qmin1val*bval2));

        %will be positive
        fb = eta1*(nb^eta2) - (Wmin1val*zval/(Wmin1val*zval*nb+bval1-Qmin1val*bval2));

        %actually perform bisection
        for nct=1:nmaxit;

            %test for convergence
            if (max(abs(nb-na))<ntol); 
                break; 
            end;

            %if not yet converged, choose midpoint candidate value
            nx = (na+nb)/2;
            fx = eta1*(nx^eta2) - (Wmin1val*zval/(Wmin1val*zval*nx+bval1-Qmin1val*bval2));

            %classify errors, noting that fx is strictly increasing
            if fx >= 0 
                nb = nx;
            else
                na = nx;
            end
        end

        nval =  (na + nb) / 2;

        %dense policy functions
        ndensepol(zcount, bmin1count) = nval;
        cdensepol(zcount, bmin1count) = (bval1 + Wmin1val * zval * nval - Qmin1val * bval2 );
       
    end            
end

CB = 0;
NB = 0;
Br = 0;
DeltaB = 0;
tempcount = 0;

for zcount = 1:znum
    
    zval = z0(zcount);
    
    for bmin1count = 1:bdensenum
    
        bval1 = bdense0(bmin1count);
        bval2 = bval2mat(zcount, bmin1count);
            
        for dscount = 1:dsnum
            
            tempcount = tempcount+1;
            
            if dscount == 2
                
               CB = CB +  Wmin1val*zval*defls * distmin1(tempcount);
               NB = NB + defls*zval*distmin1(tempcount);
                                           
            else
                xival = xival2mat(zcount, bmin1count);
                
                if xival < 0 
                    pdefault = 0;
                elseif xival > xibar 
                    pdefault = 1;
                else
                    pdefault = xival / xibar;
                end
                
                %calculate labor supply with golden section search
                               
                nval =  ndensepol(zcount, bmin1count);        
                cval = cdensepol(zcount, bmin1count);
                
                CB = CB + cval * distmin1(tempcount) *(1 - pdefault);
                CB = CB + Wmin1val*zval*defls * distmin1(tempcount) *pdefault;
                
                NB = NB + zval*nval*distmin1(tempcount)*(1 - pdefault);
                NB = NB + defls*zval*distmin1(tempcount)*pdefault;
                
                bval = splintfunc(b0, bvalmat(zcount,:), bval2mat(zcount,:), bdense0(bmin1count));
                
                Br = Br + bval * distmin1(tempcount)*(1 - pdefault);
                DeltaB = DeltaB + bval * distmin1(tempcount)*(pdefault);
                
                
            end            
        end
    end
end

%Fraction of debt defaulted on
DeltaB = DeltaB/Br;
   
%FRHS(VDind)     = VDmin1;
%FRHS(VNDind)    = VNDmin1;
%FRHS(Vind)      = Vind1;
%FRHS(XIind)     = XIind1;
%FRHS(distind)     = distmin1;
Ckmin1       = Ymin1 - etac * CB - Gmin1 - K + Kmin1 - delta * xik * ( ( K - Kmin1) / Kmin1)^2 * Kmin1 / 2 + Delta * Brkmin1; 
FRHS(Ckind)  = Ckmin1;
FRHS(NKind)  = ( (Wmin1val * ZKmin1) / (eta1 * Ckmin1) )^(1 / eta2);
FRHS(Zkind)  = rho_zk * log(ZKmin1) + eps(epskzind);
FRHS(Lambda_t_t1ind)   = betak * CK / Ckmin1;
FRHS(Brkind) = Br;
FRHS(Bind)  = - (1 - DeltaB - Delta)/(1 + Pimin1) & Brkmin1 - (1 + Imin1) / (1 + Pimin1) * Bgmin1 + Wmin1val * ZKmin1 * NKmin1 + Dkmin1 + Dgmin1 - Tkmin1 - Ckmin1;
FRHS(Qind)   = Lambda_t_t1min1 * (1 - DeltaB - Delta) / (1 + Pimin1) + eta(etaBRind);
FRHS(Lambdaind)  = Lambdamin1 * (1 + Imin1)/(1 + Pimin1) + eta(etaBind);
FRHS(Bgind)  = -Gmin1 + (1 + Imin1)/(1 + Pimin1) * Bmin1 + Tkmin1;
FRHS(Wind)   = Mmin1 * (1 - alpha) * Amin1 * Kmin1^alpha * Lmin1^(- alpha);
FRHS(Tkind)  = Tkbar + gamma_t * log( ( Bgmin1 / Ymin1) / bybar); %changed bmin1 into Bgmin1
FRHS(Iind)   = ibar + phi_i * log(1 + Pimin1) + eps(epsiind);
FRHS(Aind)   = rho_a* log(Amin1) + eps(epsaind);
FRHS(Yind)   = Amin1 * Kmin1^(alpha) * Lmin1^(1 - alpha);
FRHS(Dgind)  = Ymin1 - Mmin1 * Ymin1;
FRHS(Mind)   = ((Rmin1 + delta)/alpha)^(alpha)* (Wmin1val/(1 - alpha))^(1-alpha);
FRHS(Piind)  = M / Mmin1;
FRHS(Dkind)  = (1 + Rmin1)*Kmin1 - K - xik * ( (K - Kmin1)/Kmin1)^2 * Kmin1 / 2;
FRHS(Lind)   = ZKmin1 * NKmin1 + etac * NB;
FRHS(Kind)   = Lambda_t_t1 * ( R + 1 + xik * ( IK / K - delta) * (IK / K - delta + 1) - xik/2 * (IK / K - delta)^2) + eta(etaKind);
FRHS(Gind)   = -log(Ymin1) + log(gybar) - gamma_g * log(Bgmin1 / Ymin1);
FRHS(deltaind)   = DeltaB;
FRHS(Rind)    = Mmin1 * alpha * Amin1 * Kmin1^(alpha - 1) * Lmin1^(1 - alpha) - delta;
FRHS(Ikind)   = K - (1 - delta) * Kmin1;


F = FLHS - FRHS; 

end