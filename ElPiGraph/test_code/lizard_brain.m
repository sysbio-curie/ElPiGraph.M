function [data,irx,v,u,s,np,ed] = lizard_brain(number_of_branches_minus_2,dimension,varargin)

fprintf('===========================================================\n');
fprintf('LIZARD BRAIN FUNCTION FROM test_code FOLDER OF ELPIGRAPH!!!\n');
fprintf('===========================================================\n');

rand('twister',10)
epsilon = 0.005;
add_noise = 0.02;
min_branch_points = 50;
k_forknngraph = 8;
make_knn_graph = 0;
noise_type = 'laplacian';
make_pca_plot = 1;

for i=1:2:length(varargin)
   if strcmpi(varargin{i},'epsilon')
	epsilon = varargin{i+1};
   elseif strcmpi(varargin{i},'add_noise')
	add_noise = varargin{i+1};
   elseif strcmpi(varargin{i},'min_branch_points')
	min_branch_points = varargin{i+1};
   elseif strcmpi(varargin{i},'k_forknngraph')
	k_forknngraph = varargin{i+1};
   elseif strcmpi(varargin{i},'make_knn_graph')
	make_knn_graph = varargin{i+1};	
   elseif strcmpi(varargin{i},'noise_type')
	noise_type = varargin{i+1};
   elseif strcmpi(varargin{i},'make_pca_plot')
	make_pca_plot = varargin{i+1};
   elseif strcmpi(varargin{i},'random_seed')
	rand('twister',varargin{i+1});
   end
end

irx = [];

x0 = zeros(1,dimension);
i1 = 1;
i2 = 2;
branch = [];
while size(branch,1)<min_branch_points
x0(i1) = rand();
x0(i2) = rand();
branch = make_branch(x0,i1,i2,epsilon);
end

data = branch;
irx(1:size(data,1)) = 1;
ed = [];
for i=1:size(data,1)-1
    ed(i,1) = i;
    ed(i,2) = i+1;
end


k = 0;
while k<=number_of_branches_minus_2
n = floor(size(branch,1)/2);
%x0 = branch(n,:);
new_start = floor(rand()*size(data,1)+1);
x0 = data(new_start,:);

i1 = floor(rand()*dimension+1);
i2 = floor(rand()*dimension+1);
while(i2==i1)
    i2 = floor(rand()*dimension+1);
end
%disp(sprintf('Dim (%i,%i)',i1,i2));

[newbranch] = make_branch(x0,i1,i2,epsilon);
n1 = size(data,1);
n2 = size(newbranch,1);

if n2>min_branch_points-1
data(n1+1:n1+n2,:) = newbranch(:,:);
irx(n1+1:n1+n2) = k+2;
branch = newbranch;
k = k+1;

p = size(ed,1)+1;
 ed(p,1) = new_start;
 ed(p,2) = n1+1;
 for i=n1+1:n1+n2-1
     p = p+1;
     ed(p,1) = i;
     ed(p,2) = i+1;
 end

end

% plot(branch(:,1),branch(:,2),'ko'); hold on;
%  plot([x0(:,1) x0(:,1)+v1(:,1)/20],[x0(:,2) x0(:,2)+v1(:,2)/20],'b-');
%  plot([x0(:,1) x0(:,1)+v2(:,1)/20],[x0(:,2) x0(:,2)+v2(:,2)/20],'b-');
end

np = data;

if add_noise>0
    if strcmp(noise_type,'uniform')
	    data = data + rand(size(data,1),size(data,2))*add_noise;
    elseif strcmp(noise_type,'laplacian')
	    data = data + randl(size(data,1),size(data,2))*add_noise;
    end
end

[v,u,s] = pca(data);
if make_pca_plot
colors = [1 0 0; 0 1 0; 0 0 1; 0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0; 0.5 0.5 0.5; 0 0 0.5; 0 0.5 0; 0.5 0 0; 0 0.25 0.5; 0 0.5 0.25; 0.25 0 0.5; 0.25 0.5 0; 0.5 0 0.25; 0.5 0.25 0];
colors = colors(randperm(length(colors)),:);
number_of_classes = max(irx);
for i=1:number_of_classes
    j = rem(i,length(colors));
    plot(u(irx==i,1),u(irx==i,2),'ko','MarkerSize',2,'Color',colors(j,:)); hold on; drawnow;
end
xlabel(sprintf('PC1 %2.2f%%',s(1)/sum(s)*100));
ylabel(sprintf('PC2 %2.2f%%',s(2)/sum(s)*100));
axis equal;

mn = mean(data);
npc = np - repmat(mn,size(data,1),1);
npcu = (v'*npc')';
%plot(npcu(:,1),npcu(:,2),'ko','MarkerSize',3); hold on; drawnow;
adj_mat = zeros(size(np,1),size(np,1));
for i=1:size(ed,1)
    j = rem(irx(ed(i,1)),length(colors));
    plot([npcu(ed(i,1),1) npcu(ed(i,2),1)],[npcu(ed(i,1),2) npcu(ed(i,2),2)],'r-','LineWidth',3,'Color',colors(j,:)); hold on; 
    adj_mat(ed(i,1),ed(i,2))=1;
    adj_mat(ed(i,2),ed(i,1))=1;
end
conn = sum(adj_mat);
plot(npcu(conn>2,1),npcu(conn>2,2),'ro','MarkerSize',10,'LineWidth',3); hold on; drawnow;
end

if make_knn_graph

knngraph = knnsearch(data,data,'k',k_forknngraph);
fid = fopen('knn1.sif','w');
for i=1:size(knngraph,1)
    for k=2:size(knngraph,2)
        fprintf(fid,'%i\tna\t%i\n',knngraph(i,1),knngraph(i,k));
        %plot([u(knngraph(i,1),1) u(knngraph(i,k),1)],[u(knngraph(i,1),2) u(knngraph(i,k),2)],'b-','MarkerSize',2); hold on;
    end
    %drawnow;
end
fclose(fid);

end

end

function [branch] = make_branch(x0,i1,i2,epsilon)
dimension = size(x0,2);
v1 = zeros(1,dimension);
v2 = zeros(1,dimension);
v1(i1) = rand()-0.5; 
v1(i2) = rand()-0.5; 
v1 = v1/norm(v1);
v2(i1) = -v1(i2); 
v2(i2) = v1(i1); 

[branch] = parabolic_branch(x0,v1,v2,epsilon);
end

function [x] = parabolic_branch(x0,v1,v2,epsilon)
x = [];
t = epsilon/1000;
i = 1;
irx1 = find(v1~=0);
irx2 = find(v2~=0);
    while 1
      xn = x0+t*v1+t*t*v2;
      if (max(xn(irx1))<1)&&(min(xn(irx2))>0)
          x(i,:) = xn;
          i = i+1;
          t = t+epsilon;
      else
          break;
      end
    end
end
