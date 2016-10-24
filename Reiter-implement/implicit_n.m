%implicit labor supply bisection

b0 = importdata('b0.txt');

nmaxit = 2500;
nlbcons = 1e-5;
nub = 10.0;
ntol = 1e-5;
eta1 = 21.6;
eta2 = 2.0;
w = 0.17091561025416174;
Q = 0.97;
bnum = 25;

%last period borrowing and current borrowing are determined now
bmin1vec = b0;
bvec = zeros(bnum,1);
bvec(:) = b0(13);
zvec = bvec;
zvec(:) = 1;

%asymptote of the N equation
nlb = (Q*bvec - bmin1vec)./(w*zvec);nlb = max([nlb,zeros(bnum,1)],[],2)+nlbcons;

na = zeros(bnum,1); fa = na; nb = na; fb = na; nx = na; fx = na;

%will be negative
na = nlb;
fa = eta1*(na.^eta2) - (w*zvec./(w*zvec.*na+bmin1vec-Q*bvec));

%will be positive
nb(:) = nub;
fb = eta1*(nb.^eta2) - (w*zvec./(w*zvec.*nb+bmin1vec-Q*bvec));

%actually perform bisection
for nct=1:nmaxit;
    
    %text for convergence
    if (max(abs(nb-na))<ntol); break; end;
    
    %if not yet converged, choose midpoint candidate value
    nx = (na+nb)/2;
    fx = eta1*(nx.^eta2) - (w*zvec./(w*zvec.*nx+bmin1vec-Q*bvec));
    
    %classify errors, noting that fx is strictly increasing
    pos = fx>=0;
    neg = fx<0;
    
    %update appropriate endpoints
    nb(pos) = nx(pos);
    na(neg) = nx(neg);
    
end;
n = nx;

