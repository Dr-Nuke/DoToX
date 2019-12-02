function []= f_BoFig2PDF(h,fname)

set(h,'PaperOrientation','landscape');

set(h,'PaperPosition', [1 1 28 19]);


savefig(h,fname)
print(h, '-dpdf',fname);
set(h,'PaperOrientation','portrait');
print(h, '-dpng',fname);

end

