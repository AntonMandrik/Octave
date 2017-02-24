clear all
nmtr=4;

l0=20;
A=2.2;
g=9.81;
s=poly(0,'s');
den=[1 0.1 g./l0];
nom=[1];
nomw=poly(nom,'s','c');
denw=poly(den,'s','c');
j=sqrt(-1);
pi=3.1415;


W=syslin('c',nomw,denw);
ci=1;
wst=0.5;
G1p=cell(1,7);

//ñîçäàåì ìàòðèöó íà äèàãîíàëè êîòîðîé ñòîÿò ïåðåäàòî÷íûå ôóíêöèè îò nw
//îñòàëüíûå ýëåìåíòû ìàòðèöû - íîëþ.
clear den1;
for wn=1:1:nmtr+1
den1(wn,3)=den(1).*(wn-1).*(wn-1);
den1(wn,2)=-den(2).*(wn-1);
den1(wn,1)=den(3);
den1c(wn)=poly(den1(wn,:),'s','c');
end


clear den2;
for wn=1:1:nmtr+1
den2(wn,3)=den(1).*(wn-1).*(wn-1);
den2(wn,2)=den(2).*(wn-1);
den2(wn,1)=den(3);
den2c(wn)=poly(den2(wn,:),'s','c');
end

for Wwi=0:1:nmtr
Ww(nmtr-Wwi+1,nmtr-Wwi+1)=syslin('c',[1],den1c(Wwi+1));
Ww(nmtr+Wwi+1,nmtr+Wwi+1)=syslin('c',[1],den2c(Wwi+1));
end
//%ñîçäàåì ìàòðèöó ñ ýëåìåíòàìè -w äî ñðåäíåé ñòðîêè è w ïîñëå ñðåäíåé ñòðîêè
for i=1:1:(nmtr*2+1)
Wneg(i,i)=i-nmtr-1;
end
//ñîçäàåì ìàòðèöó ñ ýëåìåíòàìè -w äî ñðåäíåé ñòðîêè è w ïîñëå ñðåäíåé ñòðîêè

f11=0:0.5:6.3; 

ppp=2;
wi=1;

fi=1;
   f=f11(fi);
   
    
Wsin2wf=zeros(nmtr*2+1,nmtr*2+1);
for ni=3:1:(nmtr*2+1)
Wsin2wf(ni-2,ni)=(A./2).*exp(j.*0.5.*pi-j.*f);
Wsin2wf(ni,ni-2)=-(A./2).*exp(j.*0.5.*pi+j.*f);
end
Wsin2wf=-g.*(l0.^-2).*Wsin2wf;

f1=f+0.5.*pi./2;
Wcos2wf=zeros(nmtr*2+1,nmtr*2+1);
for ni=3:1:(nmtr*2+1)
Wcos2wf(ni-2,ni)=(A./2).*exp(j.*0.5.*pi-j.*f1);
Wcos2wf(ni,ni-2)=-(A./2).*exp(j.*0.5.*pi+j.*f1);
end
Wcos2wf=4.*(l0.^-1).*(s/j).*s.*Wcos2wf;
Wsin2wf=Wneg*Wsin2wf;

f2=2.*f;
Wsin4wf=zeros(nmtr*2+1,nmtr*2+1);
for ni=5:1:(nmtr*2+1)
Wsin4wf(ni-4,ni)=(A./2).*exp(j.*0.5.*pi-j.*f2);
Wsin4wf(ni,ni-4)=-(A./2).*exp(j.*0.5.*pi+j.*f2);
end
Wsin4wf=-2.*(l0.^-2).*(s/j).*s.*(A.^2).*Wsin4wf;
Wsin4wf=Wneg*Wsin4wf;


