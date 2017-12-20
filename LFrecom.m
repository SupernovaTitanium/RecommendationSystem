
function  [RMSE,MAE,Z,W] = LFrecom(R_train,R_test,parameter,data_type,fold_ind,tran)

stderr=2;
if(tran == 1)
	R_train = R_train';
	R_test = R_test';
end
fprintf(stderr,'T1=%d\n',parameter(1));
fprintf(stderr,'T2=%d\n',parameter(2));
fprintf(stderr,'SDP_rank=%d\n',parameter(3));
fprintf(stderr,'SDP_iter=%d\n',parameter(4));
fprintf(stderr,'mu=%d\n',parameter(5));
fprintf(stderr,'stepsize=%d\n',parameter(6));
fprintf(stderr,'TOL=%d\n',parameter(7));
fprintf(stderr,'Tikhonov_lambda_w=%d\n',parameter(8));
fprintf(stderr,'Tikhonov_lambda_bias=%d\n',parameter(9));
fprintf(stderr,'CG iter=%d\n',parameter(10));
fprintf(stderr,'CG eps=%f\n',parameter(11));




T1 = parameter(1);
T2 = parameter(2);
SDP_rank = parameter(3);
SDP_iter = parameter(4);
mu = parameter(5);
stepsize = parameter(6);
TOL = parameter(7);
lambda = parameter(8);
lambda2 = parameter(9);
cgiter = parameter(10);
cgeps = parameter(11);


n = size(R_train,1);
d = size(R_train,2);

Z = [];
c = [];


opts.SYM=true;

[R_train_ub,u_bias,i_bias] = bias_solver(R_train,lambda2,n,d);


name = strcat(data_type); 
if(7 ~= exist(name,'dir'))
	mkdir(name);
end
cd(name);
if(7 ~= exist(strcat('fold-',num2str(fold_ind)),'dir') & (tran ~=1))
	mkdir(strcat('trans-fold-',num2str(fold_ind)));
end
if(7 ~= exist(strcat('fold-',num2str(fold_ind)),'dir') & (tran ==1))
	mkdir(strcat('trans-fold-',num2str(fold_ind)));
end
cd ..
outname = strcat('SDP_',num2str(SDP_rank),'_',num2str(SDP_iter),'_lambda_',num2str(lambda),'_stepsize_',num2str(stepsize),'_mu_',num2str(mu)); 

if(tran==1)
	out = fopen(strcat(name,'/trans-fold-',num2str(fold_ind),'/',outname),'w');
else
	out = fopen(strcat(name,'/fold-',num2str(fold_ind),'/',outname),'w');
end

out1 = fopen('Timing1','w');

fprintf(out,'T1=%d\n',parameter(1));
fprintf(out,'T2=%d\n',parameter(2));
fprintf(out,'SDP_rank=%d\n',parameter(3));
fprintf(out,'SDP_iter%d\n',parameter(4));
fprintf(out,'mu=%d\n',parameter(5));
fprintf(out,'stepsize=%d\n',parameter(6));
fprintf(out,'TOL=%d\n',parameter(7));
fprintf(out,'Tikhonov_lambda_w=%d\n',parameter(8));
fprintf(out,'Tikhonov_lambda_bias=%d\n',parameter(9));
bRMSE = inf;
bMAE = inf;
bt=1;
bK0 = 0;






for t = 1:T1


	fprintf(stderr,'On outer iteration:%d\n',t);
	% compute A and use it to solve maxcut for Z empty and not empty
	if length(c)==0
		A = setA(R_train_ub,zeros(n,1),zeros(1,d));
		z = MixMaxCutSparseAAT(A,SDP_rank,SDP_iter);
		% Z = [Z z];
		% c = [c;0];
		% Zc = Z * diag(sqrt(c));
		% [W_o,~] = w_solver(R0_sp,R1_sp,Zc,Zc',zeros(size(Zc,2),n),lambda);
		% 'maxcut done'
	else
		Zc = Z * diag(sqrt(c));	
		% fprintf('At the start of loop\n');	

		W_o = wlinsolve(Z,c,lambda,R_train_ub,opts,out1,cgiter,cgeps);
	
		A = setA(R_train_ub,Zc,W_o);
		z = MixMaxCutSparseAAT(A,SDP_rank,SDP_iter);
		% Z = [Z z];
		% c = [c;0];
		% Zc = Z * diag(sqrt(c));
		% [W_o,~] = w_solver(R0_sp,R1_sp,Zc,Zc',zeros(size(Zc,2),n),lambda);
		% 'maxcut done'
	end

	Z = [Z z];
	c = [c;0];

	% 'maxcut done'

	



	
	%fully corrective by prox-GD


	k = length(c);
	h = diag_hessian(Z);

	eta = stepsize/(max(h)*k); %step size

	for t2 = 1:T2
		% t2
		% pause(0.2);
		% c
		% eta
		% mu
		grad_c = gradient_c_fast(R_train_ub,Z,c,lambda,opts,out1,cgiter,cgeps);
		% if(length(c)==3)
		% 	save('t1');
		% 	exit(0);				
		% end

		c = prox( c - eta*grad_c, eta*mu);

	end
	'prox done'


	tic;
	%shrink c and Z for j:cj=0
	Z = Z(:,c'>TOL);
	c = c(c>TOL);
	[obj,~,~] = LFrecomObj(R_train_ub,lambda,Z,c,opts,out1,cgiter,cgeps);
	obj
	
	%dump info
	if mod(t,100)==0
			%M = compute_M(c,Z);
			%obj = boolLassoObj(R,lambda,M,c, beta);

			% match = zeros(1,length(c));
			% for k = 1:size(Z0,2)
			% 		match_k = f(Z,Z0(:,k),n*tol_rate);
			% 		match(match_k>0)=k;
			% end
			% P = [match;c'];
			% P
			[~,ind] = sort(c,'descend');
			for K0 = min(5,length(c)):5:length(c)
				Z2 = Z(:,ind(1:min(end,K0)));
				c2 = c(ind(1:min(end,K0))); 
				[obj,R_guess,W] = LFrecomObj(R_train_ub,lambda,Z2,c2,opts,out1,cgiter,cgeps);
				% [obj,R_guess,W] = LFLassoObj(R0_sp,R1_sp,lambda,Z,c);
				R_guess = R_guess + repmat(u_bias,1,d) + repmat(i_bias',n,1);
				[RMSE,MAE] = statistic(R_train,R_guess);
				fprintf(stderr,'Training: t=%f k=%d RMSE=%f MAE=%f obj=%f\n',t,K0,RMSE,MAE,obj);	
				fprintf(out,'Training: t=%f k=%d RMSE=%f MAE=%f obj=%f\n',t,K0,RMSE,MAE,obj);	

				[RMSE,MAE] = statistic(R_test,R_guess);	

				% [acc,rec,prec,F1,F1_2,auc] = all_statistic(R_true,R_test,R_guess,Z,1,threshold);
				% [~,~,~,~,~,g_auc] = all_statistic(R_true,R_test,Rt_guess,Z,1,threshold);
				fprintf(stderr,'Testing: t=%f k=%d RMSE=%f MAE=%f obj=%f\n',t,K0,RMSE,MAE,obj);	
				fprintf(out,'Testing: t=%f k=%d RMSE=%f MAE=%f obj=%f\n',t,K0,RMSE,MAE,obj);
				if( RMSE< bRMSE)
					bRMSE = RMSE;
					bMAE = MAE;
					bt  = t;
					bK0 = K0;
				end	
				fprintf(stderr,'Best Testing(RMSE): t=%f k=%d RMSE=%f MAE=%f\n',bt,bK0,bRMSE,bMAE);
				fprintf(out,'Best Testing(RMSE): t=%f k=%d RMSE=%f MAE=%f\n',bt,bK0,bRMSE,bMAE);
				
				% save('t1');
				% exit(0);	
				% pause(2);
			end
	end
	
end	
end

function M = compute_M(c,Z)
	
	[n,k] = size(Z);
	M = zeros(n,n);
	for j = 1:k
		M = M + c(j)*Z(:,j)*Z(:,j)';
	end
end


function grad_M = gradient_M(R,M, beta)
	
	grad_M = -R .* beta .* exp(-beta*M) ./ (1-exp(-beta*M)+1e-2) + (1-R)*beta;
end


function grad_c = gradient_c(R,Z,c, beta)
	
	k = length(c);
	grad_c = zeros(k,1);
	%for j=1:k
	%	ZTzj = Z'*Z(:,j); %k by 1
	%	grad_c(j) = -Z(:,j)'*R*Z(:,j) + 0.5*c'*(ZTzj.^2) + 0.5*c(j)*ZTzj(j).^2;
	%end
	grad_M = gradient_M(R,compute_M(c,Z), beta);
	for j=1:k
		grad_c(j) = Z(:,j)'*grad_M*Z(:,j);
	end
end

function grad_c = gradient_c_fast(R_train,Z,c,lambda,opts,out1,cgiter,cgeps)

	k = length(c);
	grad_c = zeros(k,1);
	Zc = Z * diag(sqrt(c));
	W_o = wlinsolve(Z,c,lambda,R_train,opts,out1,cgiter,cgeps);
	A =  setA(R_train,Zc,W_o);


	for j=1:k
		grad_c(j) = -Z(:,j)'*A*A'*Z(:,j);
	end
	% if norm((Zc'*Zc+lambda)*W_o-(Zc'*R_train))>1e-7
	% 	save('t1');
	% 	exit(0);
	% end
end

function h = diag_hessian(Z)
	k = size(Z,2);
	h = zeros(k,1);
	for i = 1:k
		h(i) = (Z(:,i)'*Z(:,i)).^2;
	end
end

function c2 = prox( c, lambda )
	
	c2 = c;
	c2(c<=lambda)=0;
	c2(c>lambda) = c(c>lambda)-lambda;
end
function A = setA(R_train,Zc,W_o)
	% -grad_L(ZW)

	[I1,J1,V1] = find(R_train);
	E_1 = sum(Zc(I1,:).*(W_o(:,J1)'),2);
	temp1 = V1-E_1;
    % sum(temp1>=0)
    % sum(temp1<0)
  	n = size(Zc,1);
	d = size(W_o,2);
	S1 = sparse(I1,J1,temp1,n,d);
	A = -S1;
	% save('t2');
	% E = Z*W;
	% R = full(R1_sp);
	% A = -(-R.*exp(-E)+1-R)./(1+exp(-E));
end
function [obj,R_guess,W_o] =  LFrecomObj(R_train,lambda,Z,c,opts,out1,cgiter,cgeps)
	Zc = Z * diag(sqrt(c));
	W_o = wlinsolve(Z,c,lambda,R_train,opts,out1,cgiter,cgeps);
	% fprintf('At the calculation of obj\n');
	% pause(1);
	R_guess = Zc*W_o;	
	[I1,J1,V1] = find(R_train);
	E_1 = sum(Zc(I1,:).*(W_o(:,J1)'),2);
	temp1 = V1-E_1;
	obj = 0.5*norm(temp1)^2;
end

function [RMSE,MAE] = statistic(R_test,R_guess)	
	[I1,J1,V1] = find(R_test);
	ind = find(R_test);
	MAE = mean(abs(V1-R_guess(ind)));
	RMSE =sqrt(mean((V1 -R_guess(ind)).^2));
end
function [W]=wlinsolve(Z,c,lambda,R_train,opts,out,cgiter,cgeps)

	k = size(Z,2);
	d = size(R_train,2);
	Zc = Z * diag(sqrt(c));	
	tic;
	W = w_solver_recom(R_train,Zc,Zc',lambda,cgeps,cgiter);
	a=toc;
	% save('here');
	% exit(0);
	

	% C++ solver ,give 

	%Z dense
	%R_train :sparse format
	%
	% for i=1:d
		% fprintf(2,'%d\n',i);	
		% tic;	
		% r = full(R_train(:,i));

		% ind = find(r);
		% Zr = Zc(ind,:);
		% r = r(ind);
		% toc;
		% tic;
		% tic;
		% Q1=Zr'*Zr+lambda*eye(k);
		% Q2=Zr'*r;
		% a=toc;
		% fprintf(out,'Zr^t*Zr+lambda*eye(k) at k=%d: %f\n',k,a);
		
		% [W(:,i),~] = linsolve(Q1,Q2,opts);				
		% a1=toc;
		% sa1=sa1+a1;
		% fprintf(out,'linsolve at k=%d: %f\n',k,a);
		% f =  @(x)  mvprod(x,Zr,lambda,out);		
		% tic;
		% [W(:,i),~,relres,riter] = pcg(f,Zr'*r,[],5);		
		% a2=toc;
		% sa2=sa2+a2;
		% sa3=sa3+relres;
		% sa4=sa4+riter;
		% fprintf(out,'PCG at k=%d: %f\n',k,a);
		% norm(W(:,i)-temp)
		% fprintf(2,'%d %d %f\n',i,length(c),msg);
		% r = full(R_train(:,15));ind = find(r);Zr = Zc(ind,:);r = r(ind);
	% end
	 fprintf(out,'CG at k=%d: take %f secs total\n',k,a);
	% fprintf(out,'PCG at k=%d: %f rel-residual=%f iter=%f\n',k,sa2/d,sa3/d,sa4/d);
end
function [g] = mvprod(x,Zr,lambda,out)
	% tic;
	g = Zr'*(Zr*x)+lambda*x;	
	% a=toc;
	% fprintf(out,'matrix vector product %f\n',a);
end
function [R_train_ub,u_bias,i_bias] = bias_solver(R_train,lambda2,n,d)

	rate_num = nnz(R_train);
	[I1,J1,V1] = find(R_train);
	v1 = ones(rate_num*2,1);
	X = sparse([[1:rate_num]';[1:rate_num]'],[I1;n+J1],v1,rate_num,n+d);
	y = (X'*X + lambda2*eye(n+d))\(X'*V1);
	u_bias = y(1:n);
	i_bias = y(n+1:n+d); 
	R_train_ub = sparse(I1,J1,V1-u_bias(I1)-i_bias(J1),n,d);
end