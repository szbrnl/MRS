function MRS(f_a, f_b, f_c, f_f, uL, uR, n)
% obliczamy odleg³oœæ pomiêdzy kolejnymi punktami
h=1/n;
% tworzymy tablice o rozmiarze (n-1)x3
macierz = zeros(n-1,3);
% tworzymy wektor funkcji o rozmiarze (n-1)
wektor_pr = zeros(n-1,1);
% tworzymy wektor do przechowywania wyników
  wynik = zeros(n-1,1);


function [p1] = P1(k)
  p1 = (f_a((k+1)*h)-f_a((k-1)*h))/(4*h^2) + (f_a(k*h)/(h^2)) + (f_b(k*h)/(2*h));  
end;

function [p2] = P2(k)
  p2 = ((-2*f_a(k*h))/(h^2)+f_c(k*h));
end;

function [p3] = P3(k)
  p3 = (-f_a((k+1)*h)+f_a((k-1)*h))/(4*h^2) + f_a(k*h)/(h^2) - f_b(k*h)/(2*h);
end;

% wype³nianie macierzy 
macierz(1,3) = P1(1);
macierz(1,2) = P2(1);
for k = 2:n-2
    macierz(k,1) = P3(k);
    macierz(k,2) = P2(k);
    macierz(k,3) = P1(k);
end
macierz(n-1,2) = P2(n-1);
macierz(n-1,1) = P3(n-1);

wektor_pr(1,1) = f_f(1) - uL*P3(1);
wektor_pr(n-1,1) = f_f(n-1) - uR*P1(n-1);

% wype³niamy wektor prawej strony 
for k = 2:n-2
    wektor_pr(k,1) = f_f(k*h);
end




% eliminacja gausa dla macierzy (n-1)x3 i wektora prawej strony
wektor_pr(1,1) = wektor_pr(1,1)/macierz(1,2);
macierz(1,:) = macierz(1,:)/macierz(1,2);
for k = 2:n-1
    dzielnik = macierz(k,1)/macierz(k-1,2);
    macierz(k,1) = macierz(k,1) - dzielnik * macierz(k-1,2);
    macierz(k,2) = macierz(k,2) - dzielnik * macierz(k-1,3);
    wektor_pr(k,1) = wektor_pr(k,1) - dzielnik * wektor_pr(k-1,1);
    wektor_pr(k,1) = wektor_pr(k,1)/macierz(k,2);
    macierz(k,:) = macierz(k,:)/macierz(k,2);
end
% podstawienie odwrotne po eliminacji gausa
wynik(n-1,1) = wektor_pr(n-1,1);
for k = n-2:-1:1
    wynik(k,1) = wektor_pr(k,1) - macierz(k,3)*wynik(k+1,1);
end
% tworzymy wektor punktów 
punkty = [0:1/n:1];
% oraz wartosci u(xn)
wartosci = zeros(1,n+1);
for k = 2:n
    wartosci(k)=wynik(k-1);
end

%zmiana wartoœci na brzegach 

wartosci(1) = uL;
wartosci(n+1) = uR;


p = 3;
x=punkty;
y=wartosci;
% x, y - wektory wspó³rzednych
%p = stopien wielomianu 
A(1:p+1,1:p+1)=0;
b(1:p+1)=0;
%petla po punktach
for k=1:n
    % petla po wierszach
    for i=1:p+1
        %petla po kolumnach   
         for j=1:p+1
            A(i,j)=A(i,j)+x(k)^(i+j-2);
        end
        b(i)=b(i)+y(k)*x(k)^(i-1);
    end
end
a=A\b';
%rysowanie wielomianu aproksymujacego
dokl = 10000;
punkt = [0:1/dokl:1];
%petla po punktach dla ktorych rysujemy wartosci wielomianu
for i=1:dokl+1
    %petla po wspolczynnikach wielomianu
    wielomian(i)=0;
    for j=1:p+1
        wielomian(i)=wielomian(i)+a(j)*punkt(i)^(j-1);
    end
end
plot(punkt,wielomian);
hold on
plot(x,y,'rx');
hold off
end;