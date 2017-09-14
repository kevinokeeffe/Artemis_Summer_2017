set(0,'defaultaxeslinewidth',4)
set(0,'defaultlinelinewidth',4)
set(0,'defaultaxesfontsize',26)
set(0,'defaultlineMarkerSize',12)

%
D = 4:2:24;  
P = [160 340 480 650 750 830 890 920 940 950 960];

Pfit=fittype('a*(1-exp(-x.^2/b^2))','coefficients',{'a','b'});

F1=fit((D/2)',P',Pfit,'Startpoint',[145 3]);

cF1=coeffvalues(F1);
conf1=confint(F1);

W1FWHM=2*sqrt(log(2))*cF1(2)
W1do=2*sqrt(log(2))*conf1(1,2);
W1up=2*sqrt(log(2))*conf1(2,2);
err_up = W1up-W1FWHM
err_do = W1FWHM-W1do

D_l = linspace(0,25,1000);
plot(D/2,P,'o')
hold on
plot(D_l/2, feval(F1,D_l/2))
hold off
axis tight
xlabel('Aperture radius (mm)')
ylabel('Transmitted power (mW)')
set(gca,'linewidth',4)
%%
W_esqrd_diameter = 2*sqrt(2)*cF1(2)