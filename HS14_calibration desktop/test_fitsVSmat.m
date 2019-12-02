path='empty_0012.fits';
s=dir(path)
fprintf('\n fitsinfo\n')
info=fitsinfo(path)
info.PrimaryData
a=fitsread(path);
class(a)
a=a+0.01;
save('a_matfile','a')
fitswrite(a,'a_fitsfile.fits')

b=fitsread('a_fitsfile.fits');
%%
c=load('a_matfile');
c=c.a;