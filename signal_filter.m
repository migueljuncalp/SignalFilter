
%--------------------------------------------------------------------------
%  Practica de Filtrado de sañales 
%  Objetivo: filtrar una señal formada por suma de señales sinusoidales
%--------------------------------------------------------------------------


figure(1);


fs=48000; % Frecuencia de muestreo utilizada en el DVD en Hz.
Ts=1/fs; %Periodo de muestreo: inverso de la frec. muestreo en Hz
t=Ts:Ts:T;%Vector eje tiempo
T=1; %Duración de la señal en segundos
L=length(t); %Longitud del vector t. Comprobar que sea igual a fs*T

%%-----Señales------------------------

% Señal coseno_1
f1=440; %Frecuencia de la señal senoidal
A1=0.1; %Amplitud de la señal senoidal
xdet1=A1*cos(2*pi*f1*t); %Generación de x(t)

% Señal coseno_2
f2=550; %Frecuencia de la señal senoidal
A2=0.1; %Amplitud de la señal senoidal
xdet2=A2*cos(2*pi*f2*t); %Generación de x(t)

% Señal coseno_3
f3=660; %Frecuencia de la señal senoidal
A3=0.1; %Amplitud de la señal senoidal
xdet3=A3*cos(2*pi*f3*t); %Generación de x(t)

% Señal coseno_4
f4=880; %Frecuencia de la señal senoidal
A4=0.1; %Amplitud de la señal senoidal
xdet4=A4*cos(2*pi*f4*t); %Generación de x(t)

%-------------------------------------------

%señal sin ruido
xsum = xdet1+xdet2+xdet3+xdet4;

%señal con ruido
R=randn(1,L); %Generacion de ruido blanco gaussiano
xsumr = xsum + R ; %sumamos ruido blanco gausiano a la señal suma 


subplot (5, 1, 1);
plot(t,xdet1);
axis([0 0.1 -0.5 0.5]);

subplot (5, 1, 2);
plot(t,xdet2);
axis([0 0.1 -0.5 0.5]);

subplot (6, 1, 3);
plot(t,xdet3);
axis([0 0.1 -0.5 0.5]);

subplot (6, 1, 4);
plot(t,xdet4);
axis([0 0.1 -0.5 0.5]);

subplot (6, 1, 5);
plot(t,xsum);
axis([0 0.1 -0.5 0.5]);


subplot (6, 1, 6);
plot(t,xsumr);
axis([0 0.1 -0.5 0.5]);




%--------------------------------------------------
% Transformada de Fourier
%--------------------------------------------------

figure(2)

subplot(6,1,1);
%transformada de Fourier señal coseno 1
Xdef=fftshift(fft(xdet1)); % Transformada de Fourier de xdet
f= -fs/2+fs/L:fs/L:fs/2; %Declaración del eje de frecuencias
plot(f,abs(Xdef)/max(abs(Xdef))); %Dibujo del módulo de Xdef
axis([-1000 1000 0 1.1]); %Ajuste de los ejes de la figura

subplot(6,1,2)
%transformada de Fourier señal coseno 2
Xdef=fftshift(fft(xdet2)); % Transformada de Fourier de xdet
f= -fs/2+fs/L:fs/L:fs/2; %Declaración del eje de frecuencias
plot(f,abs(Xdef)/max(abs(Xdef))); %Dibujo del módulo de Xdef
axis([-1000 1000 0 1.1]); %Ajuste de los ejes de la figura

subplot(6,1,3);
%transformada de Fourier señal coseno 3
Xdef=fftshift(fft(xdet3)); % Transformada de Fourier de xdet
f= -fs/2+fs/L:fs/L:fs/2; %Declaración del eje de frecuencias
plot(f,abs(Xdef)/max(abs(Xdef))); %Dibujo del módulo de Xdef
axis([-1000 1000 0 1.1]); %Ajuste de los ejes de la figura

subplot(6,1,4);
%transformada de Fourier señal coseno 3
Xdef=fftshift(fft(xdet3)); % Transformada de Fourier de xdet
f= -fs/2+fs/L:fs/L:fs/2; %Declaración del eje de frecuencias
plot(f,abs(Xdef)/max(abs(Xdef))); %Dibujo del módulo de Xdef
axis([-1000 1000 0 1.1]); %Ajuste de los ejes de la figura

subplot(6,1,5);
%transformada de Fourier señal coseno suma
Xdef=fftshift(fft(xsum)); % Transformada de Fourier de xdet
f= -fs/2+fs/L:fs/L:fs/2; %Declaración del eje de frecuencias
plot(f,abs(Xdef)/max(abs(Xdef))); %Dibujo del módulo de Xdef
axis([-1000 1000 0 1.1]); %Ajuste de los ejes de la figura

subplot(6,1,6);
%transformada de Fourier señal coseno suma
Xdef=fftshift(fft(xsumr)); % Transformada de Fourier de xdet
f= -fs/2+fs/L:fs/L:fs/2; %Declaración del eje de frecuencias
plot(f,abs(Xdef)/max(abs(Xdef))); %Dibujo del módulo de Xdef
axis([-1000 1000 0 1.1]); %Ajuste de los ejes de la figura

%---------------------------------------------------------------
% Filtros 

% Filtro de pase bajo
flow=500; %Ejemplo de frecuencia de corte (Hz)
Hdef=zeros(1,L/2);
Hdef(1:flow/fs*L)=ones(1,flow/fs*L);
Hdef=[fliplr(Hdef) Hdef];

subplot(5,1,4);
hold on;  % sobre escribimos en el subplot
plot(f,abs(Hdef),"r");
axis([-1000 1000 0 1.1]);

figure(3);
subplot(2,1,1);
Ydef=Xdef.*Hdef;
plot(f,abs(Ydef)/max(abs(Ydef)));
axis([-1000 1000 0 1.1]);

subplot(2,1,2);
ydet=ifft(fftshift(Ydef));
ydet=real(ydet);
plot(t,ydet);
axis([0 0.1 -2 2]);


% Filtro de pase alto
fhigh=800; %Ejemplo de frecuencia de corte (Hz)
Hdef=ones(1,L/2);
Hdef(1:fhigh/fs*L)=zeros(1,fhigh/fs*L);
Hdef=[fliplr(Hdef) Hdef];




