%%%�����! ���� ��� ����� ������ ��� ������ � ��������������� ����������
%%%�������, �� ������� ����� ��� ����������� ����, ���� �� ����������� �������, � ������ �� ������� � ��� ���� ������ �����������,
%%% ��� �������� ������ ���������� ����, ��� ���������� ������ ��������������� ��������
clear all
nmtr=6;
s=tf('s')
l0=20;
A=2.2;
g=9.81;
den=[1 0.1 g./l0];
nom=[1];
W=tf(nom,den);
ci=1;
wst=0.03;
G1p=cell(1,7);

%������� ������� �� ��������� ������� ����� ������������ ������� �� nw
%��������� �������� ������� - ����.
clear den1;
den1(1,3)=den(3);
for wn=2:1:nmtr+1
den1(wn,1)=den(1).*(wn-1).*(wn-1);
den1(wn,2)=-den(2).*(wn-1);
den1(wn,3)=den(3);
end

clear den2;
den2(1,3)=den(3);
for wn=2:1:nmtr+1
den2(wn,1)=den(1).*(wn-1).*(wn-1);
den2(wn,2)=den(2).*(wn-1);
den2(wn,3)=den(3);
end

for Wwi=0:1:nmtr
Ww(nmtr-Wwi+1,nmtr-Wwi+1)=tf([1],den1(Wwi+1,:));
Ww(nmtr+Wwi+1,nmtr+Wwi+1)=tf([1],den2(Wwi+1,:));
end
%������� ������� �� ��������� ������� ����� ������������ ������� �� nw
%��������� �������� ������� ����� ����.

%������� ������� � ���������� -w �� ������� ������ � w ����� ������� ������
for i=1:1:(nmtr*2+1)
Wneg(i,i)=i-nmtr-1;
end
%������� ������� � ���������� -w �� ������� ������ � w ����� ������� ������

f11=0:0.05:6.3; 

for ppp=[2]
for wi=1:1:round((1.4./wst)) 
 
for fi=1:1:length(f11) 
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

G=-Ww*(Wsin2wf+Wcos2wf+Wsin4wf);
for gi=1:1:(nmtr.*2+1)
    for gk=1:1:(nmtr.*2+1)
    
    [magG,phaseG]=bode((G(gi,gk)),wi.*wst);
    G1(gi,gk)=magG.*exp(j.*phaseG);
    end
end


G1p{1,ppp}=G1^ppp;
mag(wi,ppp,fi)=abs(G1p{1,ppp}(nmtr+2,nmtr+2).^(1./ppp));
phase(wi,ppp,fi)=angle(G1p{1,ppp}(nmtr+2,nmtr+2).^(1./ppp));
 
end
end
   for wi=1:1:round((1.4./wst))
    Mf(wi,ppp)=f11(1);
    for fi=2:1:length(f11)
    
        if (mag(wi,ppp,1)<mag(wi,ppp,fi))
        mag(wi,ppp,1)=mag(wi,ppp,fi);
        Mf(wi,ppp)=f11(fi);
        end
    
    end
    end

ppp
end

%/������ ������������ ������� �� �������� ����� �������




%������� �� ����� ���������� ������� ������������ ������� �� �������� �����
%�������
hold off
ci=1;
wi=1:1:round(1.4./wst);
for ppp=[2]
    
ci=ci+(6.28./10);
subplot(2,1,1)
hold on
plot(wi.*wst,mag(:,ppp,1),'color',[250.*0.8.*(1.1+sin(ci))./1.7 250.*0.8.*(1.1+sin(ci+4.*pi./3))./1.7 250.*0.8.*(1.1+sin(ci+2.*pi./3))./1.7]/255);

subplot(2,1,2)
hold on
plot(wi.*wst,Mf(:,ppp),'color',[250.*0.8.*(1.1+sin(ci))./1.7 250.*0.8.*(1.1+sin(ci+4.*pi./3))./1.7 250.*0.8.*(1.1+sin(ci+2.*pi./3))./1.7]/255);
end
legend('2')
hold off
%������� �� ����� ���������� ������� ������������ ������� �� �������� �����
%�������



%���� ���������
ppp=2

for wi=1:1:round((1.4./wst)); 
Ad=0.1;
A1(wi,ppp,1)=A;
f=Mf(wi,ppp);
mag1(wi,ppp,1)=mag(wi,ppp,1);
while (abs((1-mag1(wi,ppp,1)))>0.05)

Wsin2wf=zeros(nmtr*2+1,nmtr*2+1);
for ni=3:1:(nmtr*2+1)
Wsin2wf(ni-2,ni)=(A1(wi,ppp,1)./2).*exp(j.*0.5.*pi-j.*f);
Wsin2wf(ni,ni-2)=-(A1(wi,ppp,1)./2).*exp(j.*0.5.*pi+j.*f);
end
Wsin2wf=-g.*(l0.^-2).*Wsin2wf;

f1=f+0.5.*pi./2;
Wcos2wf=zeros(nmtr*2+1,nmtr*2+1);
for ni=3:1:(nmtr*2+1)
Wcos2wf(ni-2,ni)=(A1(wi,ppp,1)./2).*((ni-2-nmtr-1)).*exp(j.*0.5.*pi-j.*f1);
Wcos2wf(ni,ni-2)=-(A1(wi,ppp,1)./2).*((ni-nmtr-1)).*exp(j.*0.5.*pi+j.*f1);
end
Wcos2wf=4.*(l0.^-1).*(s/j).*s.*Wcos2wf;





G=-Ww*(Wsin2wf+Wcos2wf);
Gp=cell(1,7);
Gp{1,1}=G;
    for pppn=2:1:ppp
    Gp{1,pppn}=Gp{1,1}*Gp{1,pppn-1};
    end

[mag1(wi,ppp,1),phase1(wi,ppp,1)]=bode((Gp{1,ppp}(nmtr+2,nmtr+2)),wi.*wst);
mag1(wi,ppp,1)=abs((mag1(wi,ppp,1).^(1./ppp)));  
Ad=(1./mag1(wi,ppp,1));
A1(wi,ppp,1)=A1(wi,ppp,1).*Ad;

        
end 
end
  %  -���� ���������






%hold on
%plot(wi.*wst,Mf(:,1))