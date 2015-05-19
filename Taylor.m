%MÉTODO DE TAYLOR DE ORDEN SUPERIOR

%- Indique el orden a utilizar               : 4
% - Introduzca la ecuación diferencial        : 'Dy=y-(x^2)+1'
% - Introduzca la condición y(a)=b            : 'y(0)=0.5'
% - Se pidió aproximación de valor 1=(Si)     : 1
% - Introduzca el valor a aproximar           : 0.75

% - Función: 1

% - Introduzca la función de trabajo          : y-(x^2)+1

% - Función: 2

% - Introduzca la función de trabajo          : y-(x^2)-(2*x)+1

% - Función: 3

% - Introduzca la función de trabajo          : y-(x^2)-(2*x)-1

% - Función: 4

% - Introduzca la función de trabajo          : y-(x^2)-(2*x)-1

% - Introduzca la condición inicial           : 0.5
% - Introduzca el valor de a                  : 0
% - Introduzca el valor de b                  : 1
% - Introduzca el tamaño de paso h            : 0.1
warning('off','MATLAB:dispacher:InexactCaseMatch')
clear all
clc
fprintf('----------------------------------\n')
fprintf('MÉTODO DE TAYLOR DE ORDEN SUPERIOR\n')
fprintf('----------------------------------\n')
fprintf('\n');
fprintf('las variables a utilizar son x ^ y\n')

syms x
syms y
eb=0;
ix=0;
e=input('Indique el orden a utilizar: ');
d=input('Introduzca la ecuación diferencial: ');
n=input('Introduzca la condición y(a)=b: ');
z=input('Se pidió aproximación de valor 1=(Si): ');
if z==1
	kz=input('Introduzca el valor a aproximar: ');
end

fprintf('\n')

for ea=1:1:e
   ix=ix+1;
   eb=eb+1;
   fprintf(' - Función: %1.0f\n\',eb);
   f1(ix)=input('Introduzca la función de trabajo: ');
   fprintf('\n')
end
   
ya=input('Introduzca la condición inicial: ');
a=input('Introduzca el valor de a: ');
b=input('Introduzca el valor de b: ');
h=input('Introduzca el tamaño de paso h: ');
fprintf('\n\n');

fprintf(' - La solución de la ecuación diferencial es : \n\n\n');

m = dsolve(d,n,'x');
pretty(m);

fprintf('\n\n');

%Condiciones para el funcionamiento de los lazos FOR

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

for sp=a:h:(b-h)
   af=sprintf(' w%1.0f = w%1.0f',i+1,i);
   fprintf('\n');
   ag=sprintf(' w%1.0f = %1.5f',i+1,w(i+1));
   fprintf('\n');
   i=1+i;
   kc=0;
   df=0;
   qc=0;
   sd=1;
   o=0;
   while sd <= e
      kc=kc+1;
      df=df+1;
      o=o+1;
      qc=qc+(((h^(df))/(factorial(df)))*(subs((f1(o)),{x,y},{t(i),w(i)})));
      ad=sprintf(' + h^%1.0f/%1.0f f%1.0f (t%1.0f, w%1.0f)',df,factorial(df),df-1,i-1,i-1);
      ah=sprintf(' + %1.8f^%1.0f/%1.0f f%1.0f (%1.8f, %1.8f)',h,df,factorial(df),df-1,t(i),w(i));
      af=[af,ad];
      ag=[ag,ah];
      sd=sd+1;         
   end
   disp(af)
   w(i+1)=w(i)+qc; 
   disp(ag)
   af=' ';
   ag=' ';
end

fprintf('\n\n');

fprintf('i \t\t ti \t\t\t wi \t\t\t\t y(t) \t\t\t\t Error\n');  


for k1=0:k3
   k2=k1+1;
   fprintf('%1.0f \t %10.5f \t %10.15f \t %10.15f \t %e',k(k2),t(k2),w(k2),v(k2),abs((v(k2)-w(k2))));
   fprintf('\n');                                      
end

fprintf('\n');

if z == 1
   
fprintf('\n');

bb=t;
im=1;

while bb(im)<kz
   im=im+1;
end

xi=[t(im-1) t(im)];

aa=t;
af=1;
ae=aa(af);
ac=xi(1);

%Cálculo de los valores de f(x)
while ae < ac
   af=af+1;
   ae=aa(af);
end

af;
ag=af+1;

%Valores de f(x) 
yl=[w(af) w(ag)];

%Cálculo de los valores de la derivada
ah=f1(1);
ai=subs(ah,x,xi(1));
aj=subs(ai,y,w(af));

ak=subs(ah,x,xi(2));
al=subs(ak,y,w(ag));

am=[aj al];
%Valores de la derivada

yderiv=am;

n=2*length(xi);
D=zeros(n,n);
b=1;

format long

for a=1:2:n
    
   z(a)=xi(b);
   z(a+1)=xi(b);
   y2(a)=yl(b);
   y2(a+1)=yl(b);  
   b=b+1;
   
end

D(:,1)=y2';

x=0;

j=0;
cont=1;
while j<n-1
   
	x=x+1;
   j=j+1;
   for i=x:n-1
      if (z(i+1)-z(i+1-j))==0
      	D(i+1,j+1)=yderiv(cont);   
         cont=cont+1;   
      else
         D(i+1,j+1)=(D(i+1,j)-D(i,j))/(z(i+1)-z(i+1-j));
      end
   end
end

D

s=1;
c=0;
q=0;
p=0;
for i=2:n
   for j=i
      s=1;
      for m=1:1
         q=q+1;
         for h=1:(i-1)          
            s=s*(kz-z(i-q));
            q=q-1;
          end
         p=p+1;
         q=q+p;
      end
   end
   c=c+(D(i,j)*s);
end

for i=1:1
   for j=i
      c=c+(D(i,j));
   end
end 
fprintf('\n');


fprintf('ESCRIBIR LA ECUACION DE APROXIMACION DE NEVILLE A MANO\n');
fprintf('p(t) = a0 + a1(t - ti1) + a2(t - ti1)^2 + a3(t - ti1)^2*(t - ti2)\n\n');
fprintf('El resultado de y(%10.5f) es: %10.15f\n',kz,c);
fprintf('\n');
end