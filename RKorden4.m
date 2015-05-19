%MÉTODO DE RUNGE-KUTTA DE ORDEN CUATRO
%Ejemplo de clase (t^2y + y)dy=(ty^2 - t)dt 
% - Introduzca la ecuación diferencial        : 'Dy=(t*(y^2)-t)/((t^2)*y + y)'
% - Introduzca la condición y(a)=b            : 'y(0.25)=1.25'
% - Introduzca la función de trabajo          : (t*(y^2)-t)/((t^2)*y + y)
% - Introduzca la condición inicial           : 1.25
% - Introduzca el valor de a                  : 0.25
% - Introduzca el valor de b                  : 0.30
% - Introduzca el tamaño de paso h            : 0.01

fprintf('\n');
clear all
clc
fprintf('            -------------------------------------\n')
fprintf('            MÉTODO DE RUNGE-KUTTA DE ORDEN CUATRO\n')
fprintf('            -------------------------------------\n')
fprintf('Las variables a utilizar son "t" y "y"\n');
syms t y 
d=input(' - Introduzca la ecuación diferencial   : ');
n=input(' - Introduzca la condición y(a)=b       : ');
f1=input(' - Introduzca la función de trabajo     : ');
ya=input(' - Introduzca la condición inicial      : ');
a=input(' - Introduzca el valor de a             : ');
b=input(' - Introduzca el valor de b             : ');
h=input(' - Introduzca el tamaño de paso h       : ');

fprintf('\n\n');

fprintf(' - La solución de la ecuación diferencial es : \n\n\n');

m = dsolve(d,n,'t');
pretty(m);

fprintf('\n\n\n');

f=f1;
w(1)=ya;
ti=a;
n=(b-a)/h;
q(1)=a;
v(1)=a;
d=0;
c=0;
g=0;
e=1;

fprintf('- w0 = %1.8f ',ya);
fprintf('\n\n');

for i=0:(n)
   fprintf('---------------');		
   fprintf('\n\n');
   fprintf('- Iteración: %1.0f\n',e);
   fprintf('---------------');
   fprintf('\n\n');
   k1=h*subs(f,{t,y},{ti,w(i+1)});
   fprintf('- K1 = h * f(t%1.0f,w%1.0f)',i,i);
   fprintf('\n');
   fprintf('- K1 = %1.5f * f(%1.5f,w%1.0f)',h,ti,i);
   fprintf('\n');
   fprintf('- K1 = %2.15f',double(k1))
   fprintf('\n\n');
   k2=h*subs(f,{t,y},{(ti+(h/2)),(w(i+1)+(k1/2))});
   fprintf('- K2 = h * f(t%1.0f + h/2 , w%1.0f + K1/2)',i,i);
   fprintf('\n');
   fprintf('- K2 = %1.5f * f(%1.5f + %1.5f , w%1.0f + %1.15f)',h,ti,h/2,i,(k1)/2);
   fprintf('\n');
   fprintf('- K2 = %2.15f',double(k2))
   fprintf('\n\n');
   k3=h*subs(f,{t,y},{(ti+(h/2)),(w(i+1)+(k2/2))});
   fprintf('- K3 = h * f(t%1.0f + h/2 , w%1.0f + K2/2)',i,i);
   fprintf('\n');
   fprintf('- K3 = %1.5f * f(%1.5f + %1.5f , w%1.0f + %1.15f)',h,ti,h/2,i,(k2)/2);
   fprintf('\n');
   fprintf('- K3 = %2.15f',double(k3))
   fprintf('\n\n');
   k4=h*subs(f,{t,y},{(ti+h),(w(i+1)+k3)});
   fprintf('- K4 = h * f(t%1.0f + h , w%1.0f + K3)',i,i);
   fprintf('\n');
   fprintf('- K4 = %1.5f * f(%1.5f + %1.5f , w%1.0f + %1.15f)',h,ti,h,i,k3);
   fprintf('\n');
   fprintf('- K4 = %2.15f',double(k4))
   fprintf('\n\n');
   w(i+2)=w(i+1)+((1/6)*(k1+(2*k2)+(2*k3)+k4));
   fprintf('- w%1.0f = w%1.0f + (1/6)*(K1+2K2+2K3+K4)',i+1,i)
   fprintf('\n');
   fprintf('- w%1.0f = w%1.0f + (1/6)*(%1.15f+2*%1.15f+2*%1.15f+%1.15f)',i+1,i,k1,k2,k3,k4)
   fprintf('\n');
   fprintf('- w%1.0f = %2.15f',i+1,double(w(i+2)))
   fprintf('\n\n');
   ti=ti+h;
   e=e+1;
end

fprintf('\n\n');
%Este for obtiene y guarda todos los valores de t
%También se utiliza para evaluar la ecuación diferencial
for p=a:h:b
   d=1+d;
   q(d)=p;
   v(d)=subs(m,p);
end

%Este for se usa para contabilizar las iteraciones
for s=c:1:(d-1)
   g=1+g;
   k(g)=(g-1);
end

k3=k(end);

%Presentación de los datos 
fprintf('              i             ti                wi+1                y(t)');   
fprintf('\n\n');

for k1=0:k3
   k2=k1+1;

   fprintf('              %1.0f        %10.15f        %10.15f         %10.15f',k(k2),q(k2),w(k2),v(k2));
   fprintf('\n');                                      
end

fprintf('\n');   