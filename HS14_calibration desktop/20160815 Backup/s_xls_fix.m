clear all

load('cl3xls.mat');
load('d2oxls.mat');
load('empxls.mat');


for i=1:size(empty_xls,1)
  for j=[1,2,3,6,7,8,9,10], %4 and 5 are stupid time stings
      empty_log(i,j)=empty_xls{i,j};
  end
end



for i=1:size(d2o_xls,1)
  for j=[1,2,3,6,7,8,9,10], %4 and 5 are stupid time stings
      d2o_log(i,j)=d2o_xls{i,j};
  end
end


for i=1:size(cl3_xls,1)
  for j=[1,2,3,6,7,8,9,10], %4 and 5 are stupid time stings
      cl3_log(i,j)=cl3_xls{i,j};
  end
end




save('cl3_log','cl3_log');
save('d2o_log','d2o_log');
save('empty_log','empty_log');
