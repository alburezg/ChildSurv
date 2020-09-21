function out=kinship_function_DR(U,F)
%inputs
% U=subdiagonal survival matrix
% F = fertility matrix
% U and F refer to living ages
% Utilde and Ftilde, including ages at death, are constructed here,
%inside the function

%size of U
[om,om]=size(U);
%identity matrix
Iom=eye(om);

%create model with absorbing states for age at death

Mort=diag(1-sum(U));

%construct survival matrix

%must specify how deaths should be counted:
%   accumulated deaths up to age x of Focal or ...
%   experienced deaths at age x of Focal

% accumulated deaths 
type_death='cum';
Utilde=[U zeros(om,om);
    Mort Iom];

% %experienced deaths 
%    type_death='exp'
%     Utilde=[U zeros(om,om);
%         Mort zeros(om,om)];

%construct fertility matrix Ftilde
Ftilde=[F zeros(om,om);
    zeros(om,2*om)];

%identity matrix for Utilde
I2om=eye(2*om);

%find distribution of age of mothers in stable pop

Amat=U+F;
[w,d]=eig(Amat);
d=diag(d);
pick=find(d==max(d));
w=w(:,pick);
w=w/sum(w);

lambda=d(pick);

clear d

%distribution of ages of mothers of children
pi=F(1,:)'.*w;

%append zeros for the dead
pi=[pi;zeros(om,1)];
%normalize the vector
pi=pi/sum(pi);

%frequently used zero vector for initial condition
zvec=zeros(2*om,1);

%frequently used om-1 limit for iterations
omz=om-1;

%following code calculates age distribution, including the dead
%for each type of kin, and stores these as columns of an array
% e.g., a(x) = daughters at age x; A(:,x) contains a(x)

% a: daughters of focal

az=zvec;
A(:,1)=az;
for ix=1:omz
    A(:,ix+1)=Utilde*A(:,ix) + Ftilde*I2om(:,ix);
end % for ix


% b = granddaughters of Focal
b=zvec;
B(:,1)=b;
for ix=1:omz
    B(:,ix+1)=Utilde*B(:,ix) + Ftilde*A(:,ix);
end


% c = greatgranddaughters of Focal
c=zvec;
C(:,1)=c;
for ix=1:omz
    C(:,ix+1)=Utilde*C(:,ix) +Ftilde*B(:,ix);
end


% d = mothers of Focal
dzero=pi;
D(:,1)=dzero;
for ix=1:omz
    D(:,ix+1)=Utilde*D(:,ix);
end


% g = maternal grandmothers of Focal
gzero=(D*pi(1:om));
gzero(om+1:end)=0;
G(:,1)=gzero;
for ix=1:omz
    G(:,ix+1)=Utilde*G(:,ix);
end


% h = greatgrandmothers of Focal
hzero=G*pi(1:om);
hzero(om+1:end)=0;
H(:,1)=hzero;
for ix=1:omz
    H(:,ix+1)=Utilde*H(:,ix) + 0;
end


% m = older sisters of Focal
mzero=A*pi(1:om);
mzero(om+1:end)=0;
M(:,1)=mzero;
for ix=1:omz
    M(:,ix+1)=Utilde*M(:,ix) + 0;
end

% n = younger sisters
nzero=zvec;
N(:,1)=nzero;
for ix=1:omz
    N(:,ix+1)=Utilde*N(:,ix) + Ftilde*D(:,ix);
end


% p = nieces through older sisters
pzero=B*pi(1:om);
pzero(om+1:end)=0;
P(:,1)=pzero;
for ix=1:omz
    P(:,ix+1)=Utilde*P(:,ix) + Ftilde*M(:,ix);
end

% q = nieces through younger sisters
qzero=zvec;
Q(:,1)=qzero;
for ix=1:omz
    Q(:,ix+1)=Utilde*Q(:,ix) + Ftilde*N(:,ix);
end

% r = aunts older than mother
rzero=M*pi(1:om);
rzero(om+1:end)=0;
R(:,1)=rzero;
for ix=1:omz
    R(:,ix+1)=Utilde*R(:,ix) + 0;
end

% s = aunts younger than mother
szero=N*pi(1:om);
szero(om+1:end)=0;
S(:,1)=szero;
for ix=1:omz
    S(:,ix+1)=Utilde*S(:,ix) + Ftilde*G(:,ix);
end

% t = cousins from older aunts
tzero=P*pi(1:om);
tzero(om+1:end)=0;
T(:,1)=tzero;
for ix=1:omz
    T(:,ix+1)=Utilde*T(:,ix) + Ftilde*R(:,ix);
end


% v = cousins from aunts younger than mother
vzero=Q*pi(1:om);
vzero(om+1:end)=0;
V(:,1)=vzero;
for ix=1:omz
    V(:,ix+1)=Utilde*V(:,ix) + Ftilde*S(:,ix);
end %for i


%overall kinship matrices, concatenating all kin
allkin=cat(3,A,B,C,D,G,H,M,N,P,Q,R,S,T,V);

%combining older and younger categories
allkin2=cat(3,A,B,C,D,G,H,M+N,P+Q,R+S,T+V);

%output structure
out.allkin=allkin;
out.allkin2=allkin2;
out.pi=pi;
out.om=om;
out.lambda=lambda;
out.type_death=type_death;
out.U=U;
out.Utilde=Utilde;
out.F=F;
out.Ftilde=Ftilde;




