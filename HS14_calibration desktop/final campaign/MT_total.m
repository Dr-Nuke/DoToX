%% MT script: Method test


%% get some data

MT.raw.startplanes=350;
MT.raw.nplanes=5;
MT.raw.d=mt.loadData(T,F,MT);

%% reconstruct by copy paste the old formula

MT.rec1=mt.OldRec(T,F,MT);

MT.rec2=mt.NewRec1(T,F,MT);

mt.Slide2Recons(MT,2)

%%
mt.SOQuestion(F,recon)