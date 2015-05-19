%MÉTODO MODIFICADO DE EULER
%Ejemplo de clase t^2dy -tydy=y^2dt
% - Introduzca la ecuación diferencial        : 'Dy=(ty+y^2)/(t^2)'
% - Introduzca la condición y(a)=b            : 'y(1)=2'
% - Introduzca la función de trabajo          : (ty+y^2)/(t^2)
% - Introduzca la condición inicial           : 2
% - Introduzca el valor de a                  : 1
% - Introduzca el valor de b                  : 1.5
% - Introduzca el tamaño de paso h            : 0.1

fprintf('\n');
clear all
clc
fprintf('             ----------------------\n')
fprintf('             Metodo Runge Kutta de HEUN\n')
fprintf('             ----------------------\n')
fprintf('Las variables a utilizar son "x" y "y"\n');
syms x y
d=input(' - Introduzca la ecuación diferencial : ');
n=input(' - Introduzca la condición y(a)=b     : ');
f1=input(' - Introduzca la función de trabajo   : ');
ya=input(' - Introduzca la condición inicial    : ');
a=input(' - Introduzca el valor de a           : ');
b=input(' - Introduzca el valor de b           : ');
h=input(' - Introduzca el tamaño de paso h     : ');

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
   v(d)=double(subs(m,p));
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
      i=i+1;
      w(i+1)=w(i)+(h/4)*(subs(f,{x,y},{t(i),w(i)}))+(0.75*h)*(subs(f,{x,y},{t(i)+(2*h/3),(w(i)+(2*h/3)*subs(f,{x,y},{t(i),w(i)}))}));
%   w(i+1)=w(i) + (h/2)*(subs(f,{x,y},{t(i),w(i)})+ subs(f,{x,y},{t(i+1),(w(i)+(h*(subs(f,{x,y},{t(i),w(i)}))))}));   
   fprintf('\n');
   fprintf('       - w%1.0f = w%1.0f + h/4*f(t%1.0f,w%1.0f)+(3h/4)*f(t%1.0f + (2h/3), w%1.0f + (2h/3)*f(t%1.0f,w%1.0f))',i,i-1,i-1,i-1,i-1,i-1,i-1,i-1);
%   fprintf('       - w%1.0f = w%1.0f + h/2 [ f(t%1.0f,w%1.0f) + f(t%1.0f,w%1.0f + h*f(t%1.0f,w%1.0f))]',i,i-1,i-1,i-1,i,i-1,i-1,i-1);
   fprintf('\n');
   fprintf('       - w%1.0f = %10.15f',i,w(i+1))
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