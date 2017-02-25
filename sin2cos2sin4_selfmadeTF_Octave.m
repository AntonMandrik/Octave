clear all
nmtr=10;
l0=20;
A=2.2;
g=9.81;
den=[1 0.1 g./l0];
nom=[1];
ci=1;
wst=0.01;
wmax=10;
G1p=cell(1,7);

for wni=1:1:(nmtr*2+1)
Wneg(wni,wni)=wni-nmtr-1;
end

f11=0:0.1:6.3; 

for ppp=[2 4 6 8]
for wi=1:1:round((wmax./wst)) 
 
for fi=1:1:length(f11) 
    f=f11(fi);
   
    
Wsin2wf=zeros(nmtr*2+1,nmtr*2+1);
for ni=3:1:(nmtr*2+1)
Wsin2wf(ni-2,ni)=(A./2).*exp(j.*0.5.*pi-j.*f);
Wsin2wf(ni,ni-2)=-(A./2).*exp(j.*0.5.*pi+j.*f);
end
Wsin2wf=-g.*(l0.^-2).*Wsin2wf*j;

f1=f+0.5.*pi./2;
Wcos2wf=zeros(nmtr*2+1,nmtr*2+1);
for ni=3:1:(nmtr*2+1)
Wcos2wf(ni-2,ni)=(A./2).*exp(j.*0.5.*pi-j.*f1);
Wcos2wf(ni,ni-2)=-(A./2).*exp(j.*0.5.*pi+j.*f1);
end
Wcos2wf=4.*(l0.^-1).*wi.*wi.*Wcos2wf*j;
Wcos2wf=Wneg*Wcos2wf;

f2=2.*f;
Wsin4wf=zeros(nmtr*2+1,nmtr*2+1);
for ni=5:1:(nmtr*2+1)
Wsin4wf(ni-4,ni)=(A./2).*exp(j.*0.5.*pi-j.*f2);
Wsin4wf(ni,ni-4)=-(A./2).*exp(j.*0.5.*pi+j.*f2);
end
Wsin4wf=-2.*(l0.^-2).*wi*wi.*(A.^1).*Wsin4wf.*j;
Wsin4wf=Wneg*Wsin4wf;

for gi=0:1:nmtr
    Ww(nmtr+1-gi,nmtr+1-gi)=1/(-den(1)*wi*wi*gi*gi-den(2)*wi*gi*i+den(3));
    Ww(nmtr+1+gi,nmtr+1+gi)=1/(-den(1)*wi*wi*gi*gi+den(2)*wi*gi*i+den(3));
end

G=-Ww*(Wsin2wf+Wcos2wf+Wsin4wf);

mag(wi,ppp,fi)=abs(((G^ppp)(nmtr+2,nmtr+2)).^(1./ppp));
phase(wi,ppp,fi)=angle(((G^ppp)(nmtr+2,nmtr+2)).^(1./ppp));








 
end
end
   for wi=1:1:round((wmax./wst))
    Mf(wi,ppp)=f11(1);
    for fi=2:1:length(f11)
    
        if (mag(wi,ppp,1)<mag(wi,ppp,fi))
        mag(wi,ppp,1)=mag(wi,ppp,fi);
        Mf(wi,ppp)=f11(fi);
        end
    
    magp(ppp,wi)=mag(wi,ppp,1)
    
    end
    end

ppp
end
wi=1:1:round((1.4./wst));
wi=wi*wst;
save myfile.mat magp wi
