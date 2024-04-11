% Kod izole ve etkileşim altında olan anten verilerinden kompanzasyon matrisini elde 
eder. Burada Khan S. tezinde alınmış yöntem uygulanır.
% izole durumda antenlerin aldığı sinyallerin karmaşık değerleri 
bdrtudr1 = [s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 1.alıcı izole 
durumda verici dikey durumda S21(gerçel)
bditudi1= [s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 1.alıcı izole 
durumda verici dikey durumda S21(sanal)
byrtuyr1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 1.alıcı izole 
durumda verici yatay durumda S21(gerçel)
byituyi1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 1.alıcı izole 
durumda verici yatay durumda S21(sanal)
idrtydr1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 2.alıcı izole 
durumda verici dikey durumda S21(gerçel)
iditydi1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];% 2.alıcı izole 
durumda verici dikey durumda S21(sanal)
iyrtyyr1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...; % 2.alıcı izole 
durumda verici yatay durumda S21(gerçel)
iyityyi1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 2.alıcı izole 
durumda verici yatay durumda S21(sanal)
udrtudr2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 3.alıcı izole 
durumda verici dikey durumda S21(gerçel)
uditudi2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 3.alıcı izole 
durumda verici dikey durumda S21(sanal)
uyrtuyr2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 3.alıcı izole 
durumda verici dikey durumda S21(gerçel)
uyituyi2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 3.alıcı izole 
durumda verici dikey durumda S21(sanal)
ddrtydr2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 4.alıcı izole 
durumda verici dikey durumda S21(gerçel)
dditydi2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 4.alıcı izole 
durumda verici dikey durumda S21(sanal)
dyrtyyr2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 4.alıcı izole 
durumda verici yatay durumda S21(gerçel)
dyityyi2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 4.alıcı izole 
durumda verici yatay durumda S21(sanal)
“151
bdrddur1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...]; % 1.alıcı etkileşimli
durumda verici dikey durumda S21(gerçel)
bdiddui1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];% 1.alıcı etkileşimli
durumda verici dikey durumda S21(sanal)
idrdydr1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];% 2.alıcı etkileşimli
durumda verici dikey durumda S21(gerçel)
ididydi1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];% 2.alıcı etkileşimli
durumda verici dikey durumda S21(sanal)
byrduyr1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];1.alıcı etkileşimli
durumda verici yatay durumda S21(gerçel)
byiduyi1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];1.alıcı etkileşimli
durumda verici yatay durumda S21(sanal)
iyrdyyr1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];2.alıcı etkileşimli
durumda verici yatay durumda S21(gerçel)
iyidyyi1=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];2.alıcı etkileşimli
durumda verici yatay durumda S21(sanal)
udrddur2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];3.alıcı etkileşimli
durumda verici dikey durumda S21(gerçel)
udiddui2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];3.alıcı etkileşimli
durumda verici dikey durumda S21(sanal)
ddrdydr2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];4.alıcı etkileşimli
durumda verici dikey durumda S21(gerçel)
ddidydi2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];4.alıcı etkileşimli
durumda verici dikey durumda S21(sanal)
uyrduyr2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];3.alıcı etkileşimli
durumda verici yatay durumda S21(gerçel)
uyiduyi2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];3.alıcı etkileşimli
durumda verici yatay durumda S21(sanal)
dyrdyyr2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];4.alıcı etkileşimli
durumda verici yatay durumda S21(gerçel)
dyidyyi2=[s21(θ1) s21(θ2).s21(θ3)s21(θ4)s21(θ5)s21(θ6)s21(θ7). . . . . . . . . . . ...];4.alıcı etkileşimli
durumda verici yatay durumda S21(sanal)
% vericinin yatay olarak çalıştığı durumda etkileşimli S21 ile izole durumda S21’in farkı 
dfiy=dyidyyi2-dyityyi2;
dfry=dyrdyyr2-dyrtyyr2;
ufiy=uyiduyi2-uyituyi2;
“152
ufry=uyrduyr2-uyrtuyr2;
ifiy=iyidyyi1-iyityyi1;
ifry=iyrdyyr1-iyrtyyr1;
bfiy=byiduyi1-byituyi1;
bfry=byrduyr1-byrtuyr1;
% vericinin dikey olarak çalıştığı durumda etkileşimli S21 ile izole durumda S21’in farkı 
dfid=ddidydi2-dditydi2;
dfrd=ddrdydr2-ddrtydr2;
ufid=udiddui2-uditudi2;
ufrd=udrddur2-udrtudr2;
ifid=ididydi1-iditydi1;
ifrd=idrdydr1-idrtydr1;
bfid=bdiddui1-bditudi1;
bfrd=bdrddur1-bdrtudr1;
% Etkileşimli akım değerleri verici dikey durumda (Etkileşim halde voltaj değerinin yüke 
50 ohma oranı )
ibid=(bdiddui1/50);
ibrd=(bdrddur1/50);
iiid=(ididydi1/50);
iird=(idrdydr1/50);
iuid=(udiddui2/50);
iurd=(udrddur2/50);
idid=(ddidydi2/50);
idrd=(ddrdydr2/50);
% Etkileşimli akım değerleri verici yatay durumda (Etkileşim halde voltaj değerinin yüke 
50 ohma oranı )
ibiy=(byiduyi1/50);
ibry=(byrduyr1/50);
iiiy=(iyidyyi1/50);
iiry=(iyrdyyr1/50);
“153
iuiy=(uyiduyi2/50);
iury=(uyrduyr2/50);
idiy=(dyidyyi2/50);
idry=(dyrdyyr2/50);
Ay(4,3,12)=0;% yatay bileşen akım matrisleri (12 açı için, 4x3 lük akım matrisleri
by(4,1,12)=0; % yatay bileşen voltaj farkı matrisleri (12 açı için, 4x1 lik voltaj farkı 
matrisleri)
Ad(4,3,12)=0; % dikey bileşen akım matrisleri
bd(4,1,12)=0; % dikey bileşen voltaj farkı matrisleri (12 açı için, 4x1 lik voltaj farkı 
matrisleri)
for k=1:1:12;
Ay(:,:,k)=[iiry(k)+i*iiiy(k) iury(k)+i*iuiy(k) idry(k)+i*idiy(k); 
ibiy(k)+i*ibry(k)+iury(k)+i*iuiy(k) idry(k)+i*idiy(k) 0; 
iiry(k)+i*iiiy(k)+idry(k)+i*idiy(k) ibiy(k)+i*ibry(k) 0; iury(k)+i*iuiy(k) 
iiry(k)+i*iiiy(k) ibiy(k)+i*ibry(k)];
by(:,:,k)=[bfry(k)+i*bfiy(k);ifry(k)+i*ifiy(k); ufry(k)+i*ufiy(k) ; dfry(k)+i*dfry(k) ]
end
for k=1:1:12;
Ad(:,:,k)=[iird(k)+i*iiid(k) iurd(k)+i*iuid(k) idrd(k)+i*idid(k); 
ibid(k)+i*ibrd(k)+iurd(k)+i*iuid(k) idrd(k)+i*idid(k) 0; 
iird(k)+i*iiid(k)+idrd(k)+i*idid(k) ibid(k)+i*ibrd(k) 0; iurd(k)+i*iuid(k) 
iird(k)+i*iiid(k) ibid(k)+i*ibrd(k)];
bd(:,:,k)=[bfrd(k)+i*bfid(k);ifrd(k)+i*ifid(k); ufrd(k)+i*ufid(k) ; dfrd(k)+i*dfrd(k)]
end
xy(3,1,12)=0; % Ax=b için çözüm olan x değişkenleri
xd(3,1,12)=0;
for k=1:1:12
 xy(:,:,k)=inv(Ay(:,:,k)'*Ay(:,:,k))*(Ay(:,:,k)')*by(:,:,k) % tersi(yatay akım 
değerlerinin transposu * yatay akım değerleri)*(yatay voltaj farkı değerleri) 
 xd(:,:,k)=inv(Ad(:,:,k)'*Ad(:,:,k))*(Ad(:,:,k)')*bd(:,:,k) % tersi(dikey akım 
değerlerinin transposu * dikey akım değerleri)*(dikey voltaj farkı değerleri)
end
“154
Cy(4,4,12)=0; % yatay bileşen için kompanzasyon matrisi
for k=1:1:12;
Cy(:,:,k)=[1 -xy(1,1,k)/50 -xy(2,1,k)/50 -xy(3,1,k)/50 ;
-xy(1,1,k)/50 1 -xy(1,1,k)/50 -xy(2,1,k)/50 ; 
-xy(2,1,k)/50 -xy(1,1,k)/50 1 -xy(1,1,k)/50;
xy(3,1,k)/50 -xy(2,1,k)/50 -xy(1,1,k)/50 1 ];
end
Cd(4,4,12)=0; % dikey bileşen için kompanzasyon matrisi
for k=1:1:12;
Cd(:,:,k)=[1 -xd(1,1,k)/50 -xd(2,1,k)/50 -xd(3,1,k)/50 ;
-xd(1,1,k)/50 1 -xd(1,1,k)/50 -xd(2,1,k)/50 ;
-xd(2,1,k)/50 -xd(1,1,k)/50 1 -xd(1,1,k)/50;
-xd(3,1,k)/50 -xd(2,1,k)/50 -xd(1,1,k)/50 1 ];
end
uz=[];
uk=[uzz ;uzz; uzz ;uzz];
for k=1:1:12; % ışıma deseninden 12 adet açı alınmıştır.
ady(:,k)=Cy(:,:,k)*uk(:,k);% kompanze edilmiş dizi elemanları yatay ışıma deseni.
add(:,k)=Cd(:,:,k)*uk(:,k);% kompanze edilmiş dizi elemanları dikey ışıma bileşenleri.
diziy(:,k)=ady(1,k)+ady(2,k)+ady(3,k)+ady(4,k); % yatay ışıma desenlerinin toplamı.
dizid(:,k)=add(1,k)+add(2,k)+add(3,k)+add(4,k);% dikey ışıma desenlerinin toplamı.
end
dizi(:,:)=(diziy+dizid); % kompanze edilmiş dizi ışın deseni
the=1:1:12;
plot(the,dizi)
“155
EK 3: Yön Bulan ve kör nokta atayan küresel konumlama sistemi alıcısı Matlab Kodları
% Kod ilk önce antenlerden gelen verileri bir satırlık matrisler halinde s1 değişkeni 
olarak kaydediyor. Bu veriler tek bir matris içerisine konuluyor ve kompanzasyon matrisi 
ile çarpılıp verilerin kuplaj olmadan alındığı veriler haline getiriliyor. Konduktan sonra 
korelasyon matrisi hesaplanıyor ve tekil değer ayrışması ile MDL algoritması çalıştırılıyor 
elde edilen kaynak sayısı yön bulmak için ESPRIT algoritmasında değişken olarak 
kullanılıyor. Burada bulunan yönler kör nokta atanması için kör nokta atayan dizi çarpanı 
bulan algoritmaya geliyor. Dizi çarpanı elde edildikten sonra gelen sinyaller modifiye 
edilip alıcı algoritması ile uydular bulunmaya çalışılıyor. En sonda ESPRIT algoritmasına 
alternative olarak MUSIC ve CAPON algoritmaları verilmiştir. 
 % Anten dizisinden gelen sinyaller ile karıştırıcı sinyal sayısının MDL algoritması 
ile çıkarılması
M=4; % Dizideki anten sayısı
s1=[]; %1. anten veriler 
s2=[]; %2. anten veriler
s3=[]; %3. anten veriler
s4=[]; %4. anten veriler
% Anten verilerinin dalgacık dönüşümü sayesinde gürültü değerinin düşürülmesi 
[c1,l1]=wavedec(s1,7,'db16'); % dalgacık dönüşümüyle sinyalin ayrıştırılması
r1=wthresh(c1,'h',22); % Ayrıştırılan sinyalin eşik değerinden geçirilmesi
 
r2=waverec(r1,l1,'db16'); % Sinyalin Tekrar oluşturulması
[d1,f1]=wavedec(s1,7,'db16'); % dalgacık dönüşümüyle sinyalin ayrıştırılması
t1=wthresh(d1,'h',100); % Ayrıştırılan sinyalin düşük sinyal gürültü oranı(SNR) 
uydu sinyalleri için düşük eşik değerinden geçirilmesi
 
t2=waverec(t1,f1,'db16'); % Sinyalin Tekrar oluşturulması
[c2,l2]=wavedec(s2,7,'db16');
r3=wthresh(c2,'h',5);
r4=waverec(r3,l2,'db16');
[d2,f2]=wavedec(s2,7,'db16'); % dalgacık dönüşümüyle sinyalin ayrıştırılması
t3=wthresh(d2,'h',100); % Ayrıştırılan sinyalin düşük sinyal gürültü oranı(SNR) 
uydu sinyalleri için düşük eşik değerinden geçirilmesi
 
t4=waverec(t3,f2,'db16'); % Sinyalin Tekrar oluşturulması
[c3,l3]=wavedec(s3,7,'db16');
r5=wthresh(c3,'h',22);
r6=waverec(r5,l1,'db16');
[d3,f3]=wavedec(s3,7,'db16'); % dalgacık dönüşümüyle sinyalin ayrıştırılması
t5=wthresh(d3,'h',5); % Ayrıştırılan sinyalin düşük sinyal gürültü oranı(SNR) 
uydu sinyalleri için düşük eşik değerinden geçirilmesi
 
“156
t6=waverec(t5,f3,'db16'); % Sinyalin Tekrar oluşturulması
[c4,l4]=wavedec(s4,7,'db16');
r7=wthresh(c4,'h',22);
r8=waverec(r7,l4,'db16');
[d4,f4]=wavedec(s4,7,'db16'); % dalgacık dönüşümüyle sinyalin ayrıştırılması
t7=wthresh(d4,'h',5); % Ayrıştırılan sinyalin düşük sinyal gürültü oranı(SNR) 
uydu sinyalleri için düşük eşik değerinden geçirilmesi
 
t8=waverec(t7,f4,'db16'); % Sinyalin Tekrar oluşturulması
Xm=[r2+t2;r4+t4;r6+t6;r8+t8]; % Dalgacık dönüşümü ile gürültü değeri düşürülen anten 
verilerinin matris haline getirilmesi
C()*Xm = X ; % Anten verileri kompanzasyon matrisi ile çarpılarak kompanze edilir. 
R=(X*X')/(N+1); % Kovaryans matrisinin hesaplanması
[U,D,V]=svd(R); % tekil değer ayrıştırması 
e=diag(D);
% MDL algoritmasıyla kaynak sayısının belirlenmesi 
for k = 0 : M-1
la = e(k+1:M);
lam=la.^(1/(M-k));
MDL(k+1)=-(M-k)*(N+1)*log10(prod(lam)/(sum(e(k+1:M))/(M-k)))+0.5*k*(2*M-k)*log10((N+1)); 
% MDL formülü
[min1,index]=min(MDL);
index_MDL=index-1;
% Karıştırıcı sayısı belirlendikten sonra sisteme gelen sinyallerin açısının ESPRIT 
algoritması ile bulunması
% Dairesel dizide ESPRIT algoritması kullanılamadığı için dairesel diziden sanal doğrusal 
diziye dönüşüm matrisi elde edilmesi
for h=0:M 
p=1:M+1 
F(h+1,p)=exp((j*2*pi*(h-M/2)*(p-1))/M);
end 
F1=(1/sqrt(M))*F;
% Diagonal matrix 
for i=0:M 
“157
Jo(:,i+1)=1/(sqrt(M)*j^(i-M/2)*besseli(i-M/2,1.6*2*pi));
end 
J=diag(Jo);
T=F1’*J; % Dairesel diziden sanal doğrusal diziye dönüşüm matrisi
X=T*X; % Dönüşüm matrisi ile sinyal matrislerinin modifiye edilmesi 
Rxx=X*X'; 
[S V D]=svd(X*X');
V1=V(1:index_MDL,1:index_MDL); % kaynak sayısı kadar vektörün alınması
S1=S(:,1:index_MDL); 
D1=D(:,1:index_MDL);
Q=D1*V1*S1.'; 
Q1=Q(1:h-1,:);% ilk elemandan sondan bir önceki elemana kadar olan matris
Q2=Q(2:h,:); % İkinci elemandan sona kadar olan matris
P=pinv(Q1)*Q2; 
[V2 D2]= eig(P);
z=(sort(diag(D2)));
for O=1:index_MDL
Y(O)=aci(z(O,1));
acitahmini(O)=abs((Y(O))/(pi/360));
end
% Açılar belirlendikten sonra bu açılara kör nokta atayan anten dizi çarpanı hesaplanması
N = 4;
c = 300000000;
fc=1575420000;
lambda = c/fc;
d = 0.85*lambda;
r = N*d/(2*pi);
angp = -180:180;
for O = 1:1:index_MDL
[w,pos] = diffbfweights(N,r/lambda,angle_estim(O),'ArrayGeometry','UCA'); 
“158
% belirlenmiş açıya kör nokta atayan dizi faktörü hesabı 
 bp(O) = arrayfactor(pos,angp,w);
 bp=bp*bp(O);% her bir kör nokta için bulunmuş dizi çarpanlarının çarpılarak tüm kör 
noktaları kapsayan dizi çarpanı bulunması 
end
Xns=bp.*X; % Dizi matrisinin dizi çarpanı ile çarpılması
% antenden alınan verilerin dizi çarpanı ile yenilendikten sonra diziden alınan toplam 
sinyal( uyduları görüntülemek için alıcıya giden veriler) 
S=Xns(1)+Xns(2)+Xns(3)+Xns(4); % Dizi elemanlarının toplanarak dizi ışıma görüntüsünün 
bulunması
% Kör nokta atayan dizi çarpanı ile diziden elde edilen sinyaller yenilendikten sonra 
alıcı algoritması ile gördüğü uyduların belirlenmesi
% GPS alıcı algoritması
orneklemefrekansi= "" ; % Örnekleme frekansı
C/A-koduzunlugu=1023; % Kaba edinim (C/A) kod uzunluğu 
C/A-kodfrekansi=1.023 MHz;
toplamafrekansbandi=500e3 ;
kodornekleri = round(orneklemefrekansi / (C/A-kodfrekansi / C/A-koduzunlugu ));
% 1 ms boyutlarında sinyal verilerinin elde edilmesi
alinansinyal=S;
sinyal1 = alinansinyal(1 : kodornekleri);
sinyal2 = alinansinyal(kodornekleri+1 : 2*kodornekleri);
sinyal0DC = alinansinyal - mean(alinansinyal); 
% örnekleme periyodu
ts = 1 /orneklemefrekansi;
% Yerel taşyıcı sinyalin faz noktalarının belirlenmesi 
faznoktalari = (0 : (kodornekleri-1)) * 2 * pi * ts;
% Toplama bandı için frekans noktalarının belirlenmesi (500Hz steps)
frekanssayilari= round(toplamafrekansbandi * 2) + 1;
% Kaba edinim kodlarının üretilmesi ve sinyal örnekleme frekansında örneklenmesi 
kodornekleri = round(orneklemefrekansi /(C/A-kodfrekansi /C/A-koduzunlugu));
C/A-kodtablosu = zeros(32, kodornekleri);
“159
%--- Zaman sabitlerinin belirlenmesi --------------------------------------------------
ts = 1/orneklemefrekansi; % Örnekleme periyodu saniye
tc = 1/(C/A-kodfrekansi); % Kaba edinim kodu periyodu
%=== Tüm uydu numaraları için
for PRN = 1:32
 %--- PRN lar için kaba edinim kodlarının hesaplanması
g2s = [ 5, 6, 7, 8, 17, 18, 139, 140, 141, 251,252, 254, 255, 256, 257, 258, 469, 470, 
471, 472,473, 474, 509, 512, 513, 514, 515, 516, 859, 860,861, 862];
%--- Verilen PRN için sağa kaydırma ----------------------------
g2kaydirma = g2s(PRN);
%--- G1 kodunun üretilmesi -----------------------------------------------------
g1 = zeros(1, 1023);
reg = -1*ones(1, 10);
for i=1:1023
 g1(i) = reg(10);
 bithafiza= reg(3)*reg(10);
reg(2:10) = reg(1:9);
reg(1) = bithafiza;
end
%--- G2 kodun üretilmesi -----------------------------------------------------
g2 = zeros(1, 1023);
reg = -1*ones(1, 10);
for i=1:1023
g2(i) = reg(10);
bithafiza = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10);
reg(2:10) = reg(1:9);
reg(1) = bithafiza;
end
%--- G2 kodunun kaydırılması --------------------------------------------------------
%g2 kodunun önceki ve sonraki verilerinin yer değiştirilerek birleştirilmesi
g2 = [g2(1023-g2kaydirma+1 :1023), g2(1 : 1023-g2kaydirma)];
%--- kaba edinim kodunun elde edilmesi
C/A-kod = -(g1 .* g2);
 
koddeger = ceil((ts * (1:kodornekleri)) / tc);
“160
 
koddeger(end) = 1023;
 
C/A-kodtablosu(PRN, :) = C/A-kod(koddeger);
 
end % her bir uydu PRN için PRN = 1:32
sonuclar = zeros(frekanssayilari, kodornekleri);
frekanslar = zeros(1, frekanssayilari);
tasiyicifrekansi = zeros(1, 32);
kodfazi = zeros(1, 32);
buyuklukmetrigi = zeros(1, 32);
fprintf('(');
% tüm PRN lar için test edilmesi acqSatelliteList=1:1:32;
for PRN = 1:1:32
%% sinyalin korele edilmesi
 %--- kaba edinim kodlarının ayrık fourier dönüşümünün yapılması---
C/A-kodfrekans = conj(fft(C/A-kodtablosu(PRN, :)));
%--- frekans bandı boyunca korele edilmesi
 for frekansnumaralari = 1:frekanssayilari
%--- taşıyıcı sinyalin üretilmesi (0.5kHz step) -----------
frekanslar(frekansnumaralari) = AF - (toplamafrekansbandi/2) * 1000 + 0.5e3 * 
(frekansnumaralari - 1);
 %--- yerel cos ve sin sinyallerinin üretilmesi -------------------------------
tasiyicisinus = sin(frekanslar(frekansnumaralari) * faznoktalari);
 tasiyicikosinus = cos(frekanslar(frekansnumaralari) * faznoktalari);
 %--- sinyalden taşıyıcı sinyalin ayrıştırılması -----------------------------
I1= tasiyicisinus .* sinyal1;
Q1= tasiyicikosinus .* sinyal1;
 I2= tasiyicisinus .* sinyal2;
 Q2= tasiyicikosinus .* sinyal2;
 %--- taban sinyalin frekans düzlemine alınması--------------
 IQfrekans1 = fft(I1 + j*Q1);
 IQfrekans2 = fft(I2 + j*Q2);
 %--- frekans düzleminde çarpım
 konvolusyonIQ1 = IQfrekans1 .* C/A-kodfrekans;
 konvolusyonIQ2 = IQfrekans2 .* C/A-kodfrekans;
“161
 %--- ters ayrık zaman fourier dönüşümünün yapılması ------------
 toplamasonuc1 = abs(ifft(konvolusyonIQ1)) .^ 2;
 toplamasonuc2 = abs(ifft(konvolusyonIQ2)) .^ 2;
 
 %--- yüksek gücün seçilip kaydedilmesi
 if (max(toplamasonuc1) > max(toplamasonuc2))
 sonuclar(frekansnumaralari, :) = toplamasonuc1;
 else
 sonuclar(frekansnumaralari, :) = toplamasonuc2;
 end
 
 end % frekanslar = 1: frekansnumaralari
 
 
 [enbuyuk frekansnumaralari] = max(max(sonuclar, [], 2));
 
 [enbuyuk kodfazi] = max(max(sonuclar));
 
 kodcipornekleri = round(orneklemefrekansi / kodfrekansi);
 aralik1 = kodfazi - kodcipornekleri;
 aralik2 = kodfazi + kodcipornekleri;
 
 if aralik1 < 2
 kodfaziaraligi = aralik2:(kodornekleri + aralik1);
 
 elseif aralik2 >= kodornekleri
 kodfazaraligi =(aralik2-kodornekleri):aralik1;
 else
 kodfazaraligi = [1:aralik1,aralik2:kodornekleri];
 end
 enbuyuk2 = max(sonuclar(frekansnumaralari, kodfazaraligi));
 %--- sonucun yazılması -----------------------------------------------------
 
buyukmetrik(PRN) = enbuyuk/enbuyuk2;
 
 if (enbuyuk/enbuyuk2) > esikdegeri
 fprintf('%02d ', PRN);
 
%% CAPON algoritması ile dairesel dizi için Geliş Yönü hesaplama 
clear all
close all
%Olcay Yiğit tarafından 18.10.2022 tarinde düzrenlenmiştir.
% Ege Universitesi Elektrik Elektronik Mühendisliği 
%----------------------------------------------------------------------------------------
--------------------------
“162
%----------------------------------------------------------------------------------------
--------------------------
% d elemanlar arası uzaklık
% p sinyal wayısı
% w açısal frekans
% A yönelim vektörü
% S: simule edilmiş sinyal
% Gelen sinyal X
% R : kovaryans matrisi
% NN tahmin edilen gürültü kümesi
%------------ Sistem parametreleri ve geometri ---------------------------
r = 0.18; % platform yarıçapı (m)
M = 4; % eleman sayısı
c = 3e8; % ışık hızı (m/s)
gamma = 2*pi/M*(0:M-1); % elemanlar arası açı
t = 0.1; % gözlenen zaman aralığı (sec)
Tcamp = 1/12000; % örnek aralığı
N = round(t/Tcamp); % örnek sayısı
%-------------------- Capon algoritması ----------------------------------
IR = inv(X*X'/M); % X: sinyal matrisi IR: kovaryans matrisi
for kk = 1:length(phi)
a = exp(1i*2*pi/c*f(k)*r*cos((phi(kk)-gamma)*pi/180)).'; % Yönelim vektörü
CAPON(k,kk) = 1/real(a'*IR*a);
end
end
mesh(phi,0:length(f)-1,abs(CAPON))
grid on;
xlabel('yatay (aci)');
ylabel('frekans (Hz)');
zlabel('algoritma Capon');
title('CAPON');
% MUSIC algoritması ile dairesel dizi için Geliş Yönü hesaplama 
clear all;
close all;
clc;
tic;
% Olcay Yiğit tarafından 18.10.2022 tarinde düzrenlenmiştir.
% Ege Universitesi Elektrik Elektronik Mühendisliği 
%----------------------------------------------------------------------------------------
--------------------------
r=0.13; % Yarıçap (m)
N=4; % Eleman sayısı
d=2*r*sin(pi/N); % Dizideki elemanlar arası boşluk
s=1; % Kaynak sinyal sayısı
noise_var=0;
“163
gamma=(2*pi/N)*(0:N-1); % iki alıcı arasındaki açı
fc=1575000e3; % Taşıyıcı frekans
c=3e8; % ışık hızı (m/s)
lambda=c/fc; % dalga boyu
R=X*X'; % X: alınan sinyallerden oluşan matris
[Q,D]=eig(R); % Kovaryans matrisin özdeğer ayrışması
[D,I]=sort(diag(D),1,'descend'); % en büyük s tane özdeğerin bulunması
Q=Q(:,I); % özvektörlerin sıralanması
Qs=Q(:,1:s); % Sinyal özvektörlerin elde edilmesi
Qn=Q(:,s+1:N); % Gürültü özvektörlerin elde edilmesi
theta=0:180;
phi=0:1:359;
p_MUSIC=zeros(length(theta),length(phi));
for ii=1:length(theta)
for iii=1:length(phi)
zeta=2*pi/lambda*r*sin(theta(ii)*pi/180);
A=exp(1i*zeta*cos((phi(iii)-gamma)*pi/180)).'; % Yönelim vektörü
p_MUSIC(ii,iii)=(1/(A'*(Qn*Qn')*A));
end
end
mesh(phi,theta,abs(p_MUSIC))
grid on;xlabel('\phi');ylabel('\theta');zlabel('PMUSIC');title('UCA MUSIC');
toc;