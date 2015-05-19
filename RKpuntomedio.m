%MÉTODO DEL PUNTO MEDIO
%Ejemplo de clase  (y/(t-1)dt+(log(t-1)+2y)dy=0
% - Introduzca la ecuación diferencial        : 'Dy=-y/((x-1)*(log(x-1)+2*y))'
% - Introduzca la condición y(a)=b            : 'y(2)=4'
% - Introduzca la función de trabajo          : -y/((x-1)*(log(x-1)+2*y))
% - Introduzca la condición inicial           : 4
% - Introduzca el valor de a                  : 2
% - Introduzca el valor de b                  : 2.2
% - Introduzca el tamaño de paso h            : 0.05

fprintf('\n');
clear all
clc
fprintf('             ----------------------\n')
fprintf('             MÉTODO Runge Kutta de Punto Medio\n')
fprintf('             ----------------------\n')
fprintf('Las variables a utilizar son "x" y "y"\n');
syms x y
d=input(' - Introduzca la ecuación diferencial     : ');
n=input(' - Introduzca la condición y(a)=b         : ');
f1=input(' - Introduzca la función de trabajo       : ');
ya=input(' - Introduzca la condición inicial        : ');
a=input(' - Introduzca el valor de a               : ');
b=input(' - Introduzca el valor de b               : ');
h=input(' - Introduzca el tamaño de paso h         : ');

fprintf('\n\n');

fprintf(' - La solución de la ecuación diferencial es : \n\n\n');

m = dsolve(d,n,'x');
pretty(m);

fprintf('\n\n\n');

%Condiciones para el funcionamiento de los lazos FOR
f=f1;
w(1)=ya;
i=0;
t(1)=a;
v(1)=a;
d=0;
c=0;
g=0;

%Este for obtiene y guarda todos los valores de t
%También se utiliza para evaluar la ecuación diferencial
for p=a:h:b
   d=1+d;
   t(d)=p;
   v(d)=subs(m,p);
end

%Este for se usa para contabilizar las iteraciones
for s=c:1:(d-1)
   g=1+g;
   k(g)=(g-1);
end

k3=k(end);
   
%Este for obtiene los valores aproximados de solución
fprintf('      ----------------------------------------------------------------------------');
fprintf('\n');
fprintf('                               FÓRMULAS DE CADA ITERACIÓN');
fprintf('\n');
fprintf('      ----------------------------------------------------------------------------');
fprintf('\n\n');
fprintf('       - w0 = %1.5f ',ya);
fprintf('\n');

for j=a:h:(b-h)
   i=1+i;
   w(i+1)=w(i)+(h*(subs(f,{x,y},{(t(i)+h/2),w(i)+((h/2)*(((subs(f,{x,y},{t(i),w(i)})))))})));   
   fprintf('\n');
   fprintf('       - w%1.0f = w%1.0f + h f(t%1.0f + h/2,w%1.0f + h/2*f(t%1.0f,w%1.0f))',i,i-1,i-1,i-1,i-1,i-1);
   fprintf('\n');
   fprintf('       - w%1.0f = w%1.0f + %1.5f f(%1.9f + %1.5f,w%1.0f + %1.5f*f(%1.9f,w%1.0f))',i,i-1,h,t(i),(h/2),i-1,(h/2),t(i),i-1);
   fprintf('\n');
   fprintf('       - w%1.0f = %10.15f',i,w(i+1));
   fprintf('\n');
end

fprintf('\n');
fprintf('      ----------------------------------------------------------------------------');
fprintf('\n');

%Presentación de los datos 

fprintf('\n\n');

fprintf('     i          ti               wi               y(t)');   

fprintf('\n');

for k1=0:k3
   k2=k1+1;

   fprintf('     %1.0f | %10.5f | %10.15f | %10.15f',k(k2),t(k2),w(k2),v(k2));
   fprintf('\n');                                      
end

fprintf('\n');   


