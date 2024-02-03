%function deconv_olaf_condsim
clear
rand('state',sum(100*clock));
prefix='medQ_long';
load([prefix '.mat']); % three vectors of equal length: time, in, out
% slope of the linear variogram
theta=.2;   % Big numbers approximate delta correlation, but hard to converge!
% standard deviation of epistemic error (initial guess)
sigma = .231;
% length of transfer-function vector (dt remains the same)
n_h=1024;
% number of realizations
nreal=10;  % Make 100 for nice final plots
% input signal
x=in;
% output signal
y=out;
t=time;
% time increment
dt=time(2)-time(1);
% construction of Jacobian (convolution matrix)
c=dt*x;
r=dt*zeros(1,n_h);
J=toeplitz(c,r);

theta_old=0;

while(abs(theta_old-theta)/theta>0.01)
% construction of generalized covariance matrix (here time in units of "time")
c=[n_h:-1:1]*dt*theta;
Q=toeplitz(c);
C=chol(Q);
iQ=inv(Q);

% vector of indices
ii=[1:n_h]';

% Best estimate
% initialization of constraints
 hL=[];
 nL=0;
 iter=0;
 while 1==1
       iter=iter+1;
       % construction of unconstrained matrix
       JRJ=J'*J/sigma^2;
       u=ones(n_h,1);
       umat=[JRJ+iQ,JRJ*u;u'*JRJ,u'*JRJ*u];
       urhs=[J'*y/sigma^2;u'*J'*y/sigma^2];
       % matrix related to the Lagrange multipliers
       Lmat=zeros(n_h+1,nL);
       Lrhs=zeros(nL,1);
       for j=1:nL
           Lmat(hL(j),j)=1;
           Lmat(n_h+1,j)=1;
           Lrhs(j)=0;
       end
       mat=[umat,Lmat;Lmat',zeros(nL)];
       rhs=[urhs;Lrhs];
       a=diag(mat);a(n_h+2:end)=1;
       warning off;imat=inv(diag(a.^-1)*mat)*diag(a.^-1);warning on;
       sol=imat*rhs;
       h_be=sol(1:n_h)+sol(n_h+1);
       h_be(hL)=0;
       nu=sol(n_h+2:end);
       sim=J*h_be;
       sigma=sqrt((y-sim)'*(y-sim)/(length(y)-n_h+nL-1));
       disp(sprintf(['iteration %i: sigma = %8.3g, number of Lagrange ' ...
                     'multipliers %i'],[iter,sigma,nL]));
       hLold=hL;
       % set of entries that need Lagrange multiplier
       hLadd=ii(h_be<0);
       % remove entries that don't need a Lagrange multiplier anymore
       hLrem=hL(nu>0);
       hL=hL(~ismember(hL,hLrem));
       hL=union(hL,hLadd);
       if (isempty(hL)), hL=1; elseif(hL(1)~=1), hL=[1;hL]; end
       nL=length(hL);
       if isempty(setdiff(hLold,hL)) & isempty(setdiff(hL,hLold)), break, end
end    

% initialize sum of h and sum of h squared
h_all  = zeros(n_h,nreal);

%figure(1);clf;
% loop over all realizations
for ireal=1:nreal
    % unconditional realiuzation
    h_uc=C'*randn(n_h,1);
    % measurement error
    me = sigma*randn(size(y));
    % initialization of constraints
    hL=[];
    nL=0;
    iter=0;
    while 1==1
        iter=iter+1;
        % construction of unconstrained matrix
        JRJ=J'*J/sigma^2;
        u=ones(n_h,1);
        umat=[JRJ+iQ,JRJ*u;u'*JRJ,u'*JRJ*u];
        urhs=[J'*(y+me)/sigma^2-JRJ*h_uc;u'*J'*(y+me)/sigma^2-u'*JRJ*h_uc];
        % matrix related to the Lagrange multipliers
        Lmat=zeros(n_h+1,nL);
        Lrhs=zeros(nL,1);
        for j=1:nL
            Lmat(hL(j),j)=1;
            Lmat(n_h+1,j)=1;
            Lrhs(j)=-h_uc(hL(j));
        end
        mat=[umat,Lmat;Lmat',zeros(nL)];
        rhs=[urhs;Lrhs];
        a=diag(mat);a(n_h+2:end)=1;
        warning off;imat=inv(diag(a.^-1)*mat)*diag(a.^-1);warning on;
        sol=imat*rhs;
        h=sol(1:n_h)+sol(n_h+1)+h_uc;
        h(hL)=0;
        nu=sol(n_h+2:end);
        sim=J*h;
        disp(sprintf('iteration %i: number of Lagrange multipliers %i',[iter,nL]));
       figure(1)
       set(gcf,'name',sprintf('Realization %i',ireal));
       subplot(3,1,1)
       handle=plot(t,y+me,'-r',t,sim,'-k');
       set(handle(1),'markersize',2);
       xlabel('t [hr]');
       legend('meas.','sim.');
       title(sprintf('output after iteration %i',iter));
       subplot(3,1,2)
       handle=plot([0:n_h-1]*dt,h,'k');
       title(sprintf('transfer function after iteration %i',iter));
       xlabel('\tau [hr]');
       drawnow
        hLold=hL;
        % set of entries that need Lagrange multiplier
        hLadd=ii(h<0);
        % remove entries that don't need a Lagrange multiplier anymore
        hLrem=hL(nu>0);
        hL=hL(~ismember(hL,hLrem));
        hL=union(hL,hLadd);
        if (isempty(hL)), hL=1; elseif(hL(1)~=1), hL=[1;hL]; end
        nL=length(hL);
        if isempty(setdiff(hLold,hL)) & isempty(setdiff(hL,hLold)), break, end
    end    
    % plot realization
    t_h=dt*[0:n_h-1]';
    h_all(:,ireal) = h;
    subplot(3,1,3)
    hand=plot(t_h,mean(h_all(:,1:ireal),2),'r');
    set(hand,'linewidth',1.5);
    xlabel('\tau [hr]');
    hold on
    plot(t_h,prctile(h_all(:,1:ireal),[10 90],2),'b');
    plot(t_h,[min(h_all(:,1:ireal),[],2),max(h_all(:,1:ireal),[],2)],'b:');
    legend('mean','10%','90%','min','max');
    hold off
end

theta_old=theta;
theta=exp(fminsearch(@(lntheta) sumprob(h_all,nreal,n_h,lntheta),log(theta)))
h=h_be;
save([prefix '_transfer_func_condreal.mat'],'t_h','h','h_all','theta','sigma');
end

%Some stats:
h_trim=h(1:200);
time=dt*(0:length(h_trim)-1)';
m_0=dt*sum(h_trim)
m_1=(dt/m_0)*sum(time.*h_trim)
var=(dt/m_0)*sum(((time-m_1).^2).*h_trim) 

%  Generate an Inverse Gaussian with same parameters

dx=10;              % m downstream of the upstream
v=dx/m_1;           % m/hr
Disp=var*v^3/dx;    % m^2/hr
(var/m_1)*v^3/dx/2

%Disp=650;

tplot=linspace(0,1,500); dt=tplot(2)-tplot(1);

InvG = (dx./sqrt(2*pi*Disp*tplot.^3)).*exp((dx-v*tplot).^2./(-2*Disp*tplot));
InvG(1)=0;


% plot final results
figure(2)
plot(t_h(1:end),h(1:end),'rx','MarkerSize',4);
xlabel('\tau [hr]');
ylabel('g [1/hr]');
hold on
plot(t_h,mean(h_all(:,1:ireal),2),'r');
plot(t_h,prctile(h_all(:,1:ireal),10,2),'b');
%plot(t_h,min(h_all(:,1:ireal),[],2),'c');
%cross=plot(t_h(5:10:end),h(5:10:end),'rx');set(cross,'markersize',14);
plot(t_h,prctile(h_all(:,1:ireal),90,2),'b');
%plot(t_h,max(h_all(:,1:ireal),[],2),'c');
plot(tplot,InvG,'b--')
axis([0 .5 0 40]);
%ll=legend('conditional mean','10%, 90%','min, max','best estimate');
txt ={['kernel mass: ' num2str(m_0)],['mean (hr):   ' num2str(m_1)],['Var (hr^2): ' num2str(var)]};
text(0.25,20,txt)
txt ={['v (m/hr): ' num2str(v)],['D (m^2/hr): ' num2str(Disp)]};
text(0.25,30,txt)

%set(ll,'xcolor','w','ycolor','w');
%set(gcf,'color','w','paperunits','centimeters','paperposition',[1 1 9 6],'inverthardcopy','off');
hold off


function lnpsum = sumprob(h_all,nreal,n_h,lntheta)
% probability of h according to prior statistics
theta=exp(lntheta);
lnpsum=0;
for ireal=1:nreal
    nnz=sum(h_all(:,ireal)>0);
    ng =size(h_all,1);
    nz = ng-nnz;
    lnp_i= -nnz/2*log(4*pi*theta)-(ng-1)/2*log(1);
    % loop over all zero entries
    im=1; while (h_all(im,ireal)>0), im=im+1;end
    ip=im;
    for jj=1:nz-1
        % determine next entry
        ip=ip+1;while (h_all(ip,ireal)>0), ip=ip+1;end
    end
    % loop pver all entries
    for jj=1:ng-1
        lnp_i = lnp_i - (h_all(jj+1,ireal)-h_all(jj,ireal))^2/(4*theta);
    end
    lnpsum=lnpsum-lnp_i;
end
end


