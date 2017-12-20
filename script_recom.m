
overall_mae=0;
overall_rmse=0;
user_num =  943;
movie_num = 1682;
tran = 1;
data_type = 'ml-100k';
stderr = 2;
%% parameter=[T1,T2,SDPrank,SDPiter,mu,stepsize,TOL,Tikhonov_lambda_w,Tikhonov_lambda_bias_solver,cg_iter,cg_eps]
parameter = [200,3,10,40,32000,5,1e-5,2,0.8,5,0.1];
for i=1:5	
	train_name  = strcat('./data/',data_type,'/u',num2str(i),'.base'); 
	test_name   = strcat('./data/',data_type,'/u',num2str(i),'.test'); 
	tmp_train   = load(train_name);
	tmp_test    = load(test_name);
	data_train  = sparse(tmp_train(:,1),tmp_train(:,2),tmp_train(:,3),user_num,movie_num);
	data_test   = sparse(tmp_test(:,1),tmp_test(:,2),tmp_test(:,3),user_num,movie_num);
	fprintf(stderr,'Training fold %d\n',i);
	[RMSE,MAE,~,~]=LFrecom(data_train,data_test,parameter,data_type,i,tran);
	fprintf(stderr,'fold-%d: RMSE=%f MAE=%f\n',i,RMSE,MAE);
	re = fopen(strcat('./',data_type,'/result'),'a+');
    fprintf(re,'data_source-%s,tran-%d,fold-%d,T1=%d,T2=%d,SDPrank=%f,SDPiter=%f,mu=%f,stepsize=%f,lambda1=%f,lambda2=%f\n : RMSE=%f MAE=%f\n',train_name,tran,i,parameter(1)...
    	 ,parameter(2), parameter(3), parameter(4), parameter(5), parameter(6), parameter(8), parameter(9),RMSE,MAE);
	overall_rmse = overall_rmse+RMSE;
	overall_mae  = overall_mae+MAE;
end

fprintf(stderr,'Averaged: RMSE=%f MAE=%f\n',overall_rmse/5,overall_mae/5);
re = fopen(strcat('./',data_type,'/result'),'a+');
fprintf(re,'data_source-%s,tran-%d,T1=%d,T2=%d,SDPrank=%f,SDPiter=%f,mu=%f,stepsize=%f,lambda1=%f,lambda2=%f\n Averaged: RMSE=%f MAE=%f\n',data_type,tran,parameter(1)...
    	 ,parameter(2), parameter(3), parameter(4), parameter(5), parameter(6), parameter(8), parameter(9),overall_rmse/5,overall_mae/5);



% Later use  Z or W for clustering?
% m1-100K user:943, movie:1682