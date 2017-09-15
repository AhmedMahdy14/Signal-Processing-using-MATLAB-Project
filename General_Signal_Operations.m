
clear ;close all;clc;
Fs = input('What is the Sample rate of the signal? '); %sample rate=#samples per one second
ti=input('What is the Start  of time basis? ');
te=input('What is the  end of time basis? ');
t=linspace(ti,te,(te-ti)*Fs);
num=input('What is the  Number of the break points? ');
positions=[];
xyz=zeros(1,Fs*(te-ti));
y_index=1;
for ii=1:num %get positions of breakpoints
    e=sprintf('What is the  Number of  %d positions=  \n',ii);
    disp(e);
    g=input('');
    positions(1,ii)=g;
end
%this is a structure to save properties of signals from user 17:24
p.n(1).a='Amplitude';       %DC signal   
p.n(2).a='slope';        %Ramp signal
p.n(2).b='intercept';    %Ramp signal
p.n(3).a='Amplitude';    %Exponential signal
p.n(3).b='exponent';     %Exponential signal
p.n(4).a='Amplitude';    %Sinusoidal signal
p.n(4).b='frequency';    %Sinusoidal signal
p.n(4).c='phase';        %Sinusoidal signal
%displaying all signals that user choose from it
r=sprintf('According to %d break points you should choose the specications for %d signals \n',num,num+1);
disp(r);
disp('1-DC signal')
fprintf('-%s\n',p.n(1).a);
disp('2-Ramp signal')
fprintf('-%s\n-%s\n',p.n(2).a,p.n(2).b);
disp('3-Exponential signal')
fprintf('-%s\n-%s\n',p.n(3).a,p.n(3).b);
disp('4-Sinusoidal signal')
fprintf('-%s\n-%s\n-%s\n',p.n(4).a,p.n(4).b,p.n(4).c);
index=[];
%choose the signals according to its index
for jj=1:num+1
    m=sprintf('Choose the index of %d signal \n',jj);
    disp(m);
    x=input('=');
    index(1,jj)=x; 
end
%get values of properties of signal that user choosed previously
l=[]; v=[];y=[];i=[];o=[];I=[];J=[];O=[];
for zz=1:length(index)
    switch index(zz)
        case 1
            disp('DC signal')
            yu=input('Amplitude=');
            l=[l yu];
            [p.n(1,index(zz)).a]=l;
        case 2
            disp('Ramp signal')
            ty=input('slope=');
            v=[v ty];
            [p.n(1,index(zz)).a]=v;
            re=input('intercept=');
            y=[y re];
            [p.n(1,index(zz)).b]=y;
        case 3
            disp('Exponential signal')
            fd=input('Amplitude=');
            i=[i fd];
            [p.n(1,index(zz)).a]=i;
            kj=input('exponent=');
            o=[o kj];
            [p.n(1,index(zz)).b]=o;
        otherwise
            disp('Sinusoidal signal')
            tg=input('Amplitude=');
            I=[I tg];
            [p.n(1,index(zz)).a]=I;
            lop=input('frequency=');
            O=[O lop];
            [p.n(1,index(zz)).b]=O;          
            mk=input('phase=');
            J=[J mk];
           [p.n(1,index(zz)).c]=J;
  end
end
%plot signals that user choosed in time and frequency domain
positions=[ti positions te];
aoi=0;boi=0;coi=0;doi=0;
for pp=1:length(index) 
    time=positions(pp+1)-positions(pp);
    switch index(pp)
        case 1
            aoi=aoi+1; 
            t1=linspace(positions(pp),positions(pp+1),time*Fs);
            f1=linspace(-Fs/2,Fs/2,time*Fs);
            y1=p.n(1,1).a(aoi)*ones(1,time*Fs);
            yf1=abs(fftshift(fft(y1)));
            subplot(2,1,1);
            stem(t1,y1);
            hold on;
            subplot(2,1,2);
            stem(f1,yf1);
            hold on;
			xyz(y_index:(y_index-1+(time)*Fs))=y1;
			y_index=y_index+(time)*Fs;
        case 2
            boi=boi+1;
            t2=linspace(positions(pp),positions(pp+1),(time)*Fs);
            f2=linspace(-Fs/2,Fs/2,time*Fs);
            y2=p.n(1,2).a(boi)*linspace(0,time,time*Fs)+p.n(1,2).b(boi); 
            yf2=abs(fftshift(fft(y2)));
            subplot(2,1,1);
            stem(t2,y2);
            hold on;
            subplot(2,1,2);
            stem(f2,yf2);
            hold on;
			xyz(y_index:(y_index-1+(time)*Fs))=y2;
			y_index=y_index+(time)*Fs;
        case 3
            coi=coi+1;
            t3=linspace(positions(pp),positions(pp+1),(time)*Fs);
            f3=linspace(-Fs/2,Fs/2,time*Fs);
            y3=p.n(1,3).a(coi)*exp(p.n(1,3).b(coi)*linspace(0,time,time*Fs)); 
            yf3=abs(fftshift(fft(y3)));
            subplot(2,1,1);
            stem(t3,y3);
            hold on;
            subplot(2,1,2);
            stem(f3,yf3);
            hold on;
			xyz(y_index:(y_index-1+(time)*Fs))=y3;
			y_index=y_index+(time)*Fs;
        otherwise
            doi=doi+1;
            t4=linspace(positions(pp),positions(pp+1),(time)*Fs);
            f4=linspace(-Fs/2,Fs/2,time*Fs);
            y4=p.n(1,4).a(doi)*cos(2*pi*p.n(1,4).b(doi)*linspace(0,time,time*Fs) + p.n(1,4).c(doi)*(pi/180));
            yf4=abs(fftshift(fft(y4)));
            subplot(2,1,1); 
            stem(t4,y4);
            hold on;
            subplot(2,1,2);
            stem(f4,yf4);
            hold on;
			xyz(y_index:(y_index-1+(time)*Fs))=y4;
			y_index=y_index+(time)*Fs;
    end
 end 
subplot(2,1,1);
xlim([ti te]);
ylabel('plotting in time domain','color','k','FontSize',13);
xlabel('n','color','k','FontSize',13);
subplot(2,1,2);
ylabel('plotting in frequency domain','color','k','FontSize',13);
xlabel('f','color','k','FontSize',13);   
