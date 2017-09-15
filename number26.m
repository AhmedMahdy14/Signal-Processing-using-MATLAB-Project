%Ahmed Saed Ahmed Elmahdy 26
%Omar Mahmoud Mahmoud Amin 141
%Amr Ahmed Mohamed Mostafa 143
%Mohamed Ali Maher Mohamed 192
%Mohamed Ali Mohamed Ali 193
%Mohamed Kamed Ali Abd Elal 194
close all
clear all
clc

%% Part I: Sound Signal

[long, fs] = wavread('long_track') ;%reading short track
[short, fs1] = wavread('short_track'); % reading long track
sound (long,fs); %playing the long signal
sound(short,fs1); %playing the short signal 
long_upsampled = upsample(long,2); %upsampling the long with douple
short_downsampled = downsample(short,3); %downsampling the short with 1/3
wavwrite(long_upsampled,2*fs,'long1'); % saving the long wave signal 
wavwrite(short_downsampled,fs1/3,'short1'); % saving the short wave signal
sound(long_upsampled,2*fs); %playing the long_upsampled signal
sound(short_downsampled,fs1/3); %playing the short_downsampled signal 
%start of step 4 (fading)
t_long=linspace(0,12.5388,12.5388*fs); %time base of long signal
long_row= long'; % rotate the long signal to be row matrix
figure
plot(t_long,long_row); %plot the long signal 
 xlim([0 12.5388])
%the time of long song exactly= max index of long/fs=138240/11025= 12.5388
t=linspace(0,2,(2-0)*fs); %time base of fade in signal
fade_scale_in =t*.5;  % create fade
t_new=linspace(10.5388,12.5388,2*fs);
fade_scale_out = 6.2694-.5*t_new; % create fade
one_signal = ones(1,(10.5388-2)*fs); % the mid one signal
fade_signal=[fade_scale_in , one_signal , fade_scale_out] ; %the fading filter signal 
figure
plot(t_long,fade_signal); %plotting the total fade filter it self
output_signal= fade_signal.*long_row; % the output faded signal
figure 
plot(t_long,output_signal); % plotting the output faded signal
xlim([0 12.5388])
sound(output_signal,fs); %playing the o/p faded signal
sound (long_row); %playing the original signal
% starting the 5th step
%the time of short song exactly = max index of short/fs1=94208/11025=8.5449
long_trimmed = long_row(1:8.5449*fs+1); %timming the long signal to be the same of short signal
t_5=linspace(0,8.5449,8.5449*fs+1); %time base of long signal
short_row = short'; %rotating the short signal to be row matrix
figure
plot(t_5,short_row);
figure
plot(t_5,long_trimmed);
freq=.5;
sin_filter=sin(2*pi*freq*t_5); % the sine filter signal
figure
plot(t_5,sin_filter) ; %plotting the filter signal
for qq=1:length(sin_filter) %this loop to switch between long_trimmed in positive sign and short_row in negative sign
   if sin_filter(qq)>0
    output_alt_signal(1,qq)=long_trimmed(qq);
   else
    output_alt_signal(1,qq)=short_row(qq);   
   end
end
figure
plot(t_5,output_alt_signal); %plotting the output signal
sound(output_alt_signal) ; %play the output alternative signal

%% Part II:2.2.1 General Signal Operations

% clear ;close all;clc;
Fs = input('What is the Sample rate of the signal? '); %sample rate=#samples per one second
ti=input('What is the Start  of time basis? ');
te=input('What is the  end of time basis? ');
t_tot=linspace(ti,te,(te-ti)*Fs);
f_tot=linspace(-Fs/2,Fs/2,(te-ti)*Fs);
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
            subplot(2,2,1);
            stem(t1,y1);
            hold on;
            subplot(2,2,3);
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
            subplot(2,2,1);
            stem(t2,y2);
            hold on;
            subplot(2,2,3);
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
            subplot(2,2,1);
            stem(t3,y3);
            hold on;
            subplot(2,2,3);
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
            subplot(2,2,1); 
            stem(t4,y4);
            hold on;
            subplot(2,2,3);
            stem(f4,yf4);
            hold on;
			xyz(y_index:(y_index-1+(time)*Fs))=y4;
			y_index=y_index+(time)*Fs;
    end
 end 
subplot(2,2,1);
xlim([ti te]);
ylabel('plotting in time domain','color','k','FontSize',13);
xlabel('n','color','k','FontSize',13);
subplot(2,2,3);
ylabel('plotting in frequency domain','color','k','FontSize',13);
xlabel('f','color','k','FontSize',13);   
%2.2.2 Create a system:

output=input('get the output of the system via:\n1)impulse response\n2)transfer function\n');
%the user choosed the impulse response
if output==1	
	%here we should repeat all steps in the previuos code..to get h
	%from user ..then conv(input(from previously cod),h)
    %the next code is the same code from 64:215
Fs2 = input('What is the Sample rate of the impulse? '); %sample rate=#samples per one second
ti2=input('What is the Start  of time basis? ');
te2=input('What is the  end of time basis? ');
t_tot2=linspace(ti2+ti,te2+te,((te-ti)*Fs+(te2-ti2)*Fs2)-1 );
f_tot2=linspace(-Fs2/2,Fs2/2,((te-ti)*Fs+(te2-ti2)*Fs2)-1);
num2=input('What is the  Number of the break points? ');
positions2=[];
impulse=zeros(1,Fs2*(te2-ti2));
y_index2=1;
for ii2=1:num2 %get positions of breakpoints
    e2=sprintf('What is the  Number of  %d positions=  \n',ii2);
    disp(e2);
    g2=input('');
    positions2(1,ii2)=g2;
end
%this is a structure to save properties of signals from user 240:247
f.n(1).a='Amplitude';       %DC signal   
f.n(2).a='slope';        %Ramp signal
f.n(2).b='intercept';    %Ramp signal
f.n(3).a='Amplitude';    %Exponential signal
f.n(3).b='exponent';     %Exponential signal
f.n(4).a='Amplitude';    %Sinusoidal signal
f.n(4).b='frequency';    %Sinusoidal signal
f.n(4).c='phase';        %Sinusoidal signal
%displaying all signals that user choose from it
r2=sprintf('According to %d break points you should choose the specications for %d signals \n',num2,num2+1);
disp(r2);
disp('1-DC signal')
fprintf('-%s\n',f.n(1).a);
disp('2-Ramp signal')
fprintf('-%s\n-%s\n',f.n(2).a,f.n(2).b);
disp('3-Exponential signal')
fprintf('-%s\n-%s\n',f.n(3).a,f.n(3).b);
disp('4-Sinusoidal signal')
fprintf('-%s\n-%s\n-%s\n',f.n(4).a,f.n(4).b,f.n(4).c);
index2=[];
%choose the signals according to its index
for jj2=1:num2+1
    m2=sprintf('Choose the index of %d signal \n',jj2);
    disp(m2);
    x2=input('=');
    index2(1,jj2)=x2; 
end
%get values of properties of signal that user choosed previously
l2=[]; v2=[];y2=[];i2=[];o2=[];I2=[];J2=[];O2=[];
for zz2=1:length(index2)
    switch index2(zz2)
        case 1
            disp('DC signal')
            yu2=input('Amplitude=');
            l2=[l2 yu2];
            [f.n(1,index2(zz2)).a]=l2;
        case 2
            disp('Ramp signal')
            ty2=input('slope=');
            v2=[v2 ty2];
            [f.n(1,index2(zz2)).a]=v2;
            re2=input('intercept=');
            y2=[y2 re2];
            [p.n(1,index2(zz2)).b]=y2;
        case 3
            disp('Exponential signal')
            fd2=input('Amplitude=');
            i2=[i2 fd2];
            [f.n(1,index2(zz2)).a]=i2;
            kj2=input('exponent=');
            o2=[o2 kj2];
            [f.n(1,index2(zz2)).b]=o2;
        otherwise
            disp('Sinusoidal signal')
            tg2=input('Amplitude=');
            I2=[I2 tg2];
            [f.n(1,index2(zz2)).a]=I2;
            lop2=input('frequency=');
            O2=[O2 lop2];
            [f.n(1,index2(zz2)).b]=O2;          
            mk2=input('phase=');
            J2=[J2 mk2];
           [f.n(1,index2(zz2)).c]=J2;
  end
end
%plot signals that user choosed in time and frequency domain
positions2=[ti2 positions2 te2];
aoi2=0;boi2=0;coi2=0;doi2=0;
for pp2=1:length(index2) 
    time2=positions2(pp2+1)-positions2(pp2);
    switch index2(pp2)
        case 1
            aoi2=aoi2+1; 
            t12=linspace(positions2(pp2),positions2(pp2+1),time2*Fs2);
            f12=linspace(-Fs2/2,Fs2/2,time2*Fs2);
            y12=f.n(1,1).a(aoi2)*ones(1,time2*Fs2);
            yf12=abs(fftshift(fft(y12)));
            subplot(2,2,2);
            stem(t12,y12);
            hold on;
            subplot(2,2,4);
            stem(f12,yf12);
            hold on;
			impulse(y_index2:(y_index2-1+(time2)*Fs2))=y12;
			y_index2=y_index2+(time2)*Fs2;
        case 2
            boi2=boi2+1;
            t22=linspace(positions2(pp2),positions2(pp2+1),time2*Fs2);
            f22=linspace(-Fs2/2,Fs2/2,time2*Fs2);
            y22=f.n(1,2).a(boi2)*linspace(0,time2,time2*Fs2)+f.n(1,2).b(boi2); 
            yf22=abs(fftshift(fft(y22)));
            subplot(2,2,2);
            stem(t22,y22);
            hold on;
            subplot(2,2,4);
            stem(f22,yf22);
            hold on;
			impulse(y_index2:(y_index2-1+(time2)*Fs2))=y22;
			y_index2=y_index2+(time2)*Fs2;
        case 3
            coi2=coi2+1;
            t32=linspace(positions2(pp2),positions2(pp2+1),time2*Fs2);
            f32=linspace(-Fs2/2,Fs2/2,time2*Fs2);
            y32=f.n(1,3).a(coi2)*exp(f.n(1,3).b(coi2)*linspace(0,time2,time2*Fs2)); 
            yf32=abs(fftshift(fft(y32)));
            subplot(2,2,2);
            stem(t32,y32);
            hold on;
            subplot(2,2,4);
            stem(f32,yf32);
            hold on;
			impulse(y_index2:(y_index2-1+(time2)*Fs2))=y32;
			y_index2=y_index2+(time2)*Fs2;
        otherwise
            doi2=doi2+1;
            t42=linspace(positions2(pp2),positions2(pp2+1),time2*Fs2);
            f42=linspace(-Fs2/2,Fs2/2,time2*Fs2);
            y42=f.n(1,4).a(doi2)*cos(2*pi*f.n(1,4).b(doi2)*linspace(0,time2,time2*Fs2) + f.n(1,4).c(doi2)*(pi/180));
            yf42=abs(fftshift(fft(y42)));
            subplot(2,2,2); 
            stem(t42,y42);
            hold on;
            subplot(2,2,4);
            stem(f42,yf42);
            hold on;
			impulse(y_index2:(y_index2-1+(time2)*Fs2))=y42;
			y_index2=y_index2+(time2)*Fs2;
    end
end 
subplot(2,2,2) 
xlim([ti2 te2]);
ylabel('plotting in time domain','color','k','FontSize',12);
xlabel('n','color','k','FontSize',12);
subplot(2,2,4)
ylabel('plotting in frequency domain','color','k','FontSize',12);
xlabel('f','color','k','FontSize',12);
output_of_the_system_time=conv(xyz,impulse);
output_of_the_system_freq= abs(fftshift(fft(output_of_the_system_time)));
figure
stem(t_tot2,output_of_the_system_time);
ylabel('output of the system time','color','k','FontSize',12);
xlabel('n','color','k','FontSize',12);
figure
stem(f_tot2, output_of_the_system_freq);
ylabel('output of the system freq','color','k','FontSize',12);
xlabel('f','color','k','FontSize',12);
%the user choosed the transfer function:
elseif output==2
	num=input('What is the numerator of the transfer function? '); %enter numerator value by user
	den=input('What is the denominator of the transfer function? '); %enter denominator value by user
	figure
    zplane(num,den) %plot z-plane
    H=tf(num,den);
	[p,z]= pzmap(H); %find the values of poles and zeros of tf
    bb=0;nb=0;s=[];
    for t=1:length(p)
		s(t)=abs(p(t));
	end
    for vv=1:length(s)
       if s(vv)>1
           disp('unstable');
           break
       elseif s(vv)==1
            bb=bb+1;
           if (bb==1)&&(vv==length(s))
            disp('marginally stable') ;
           elseif bb>1
               disp('unstable');
               break
           end
       elseif s(vv)<1
            nb=nb+1;
           if (nb==length(s))&&(vv==length(s))
            disp('stable');  
           end
       end
    end
signal_time=filter(num,den,xyz); %calculate the output of system
signal_freq=abs(fftshift(fft(signal_time)));
figure
stem(t_tot,signal_time);
figure
stem(f_tot,signal_freq);
end
