%-- Programming by Navid Teymourian and Habib Izadkhah--
%-- Department of computer Science, Faculty of Mathematical science,
%   University of Tabriz, Tabriz, Iran.

%\\ The algorithm inputs are in the matrix of zero and one.
%\\ The output of the algorithm consists of clustered nodes.
%\\ Results display includes elapsed time, number of clusters, and TurboMQ
clc;
clear;
close all;
tic
%\\ Write the path of the nodes file.
A1=load('C:\Name file.txt');
b1=size(A1,1);
C=zeros(b1,b1);
for i4=1:b1
    for j4=1:b1
        if (A1(i4,j4)==1)
            C(i4,j4)=1;
            C(j4,i4)=1;
        end
    end
end
p2=1;
E=[];
for i4=1:b1
        if (C(i4,:)==0)
            E(1,p2)=i4;
            p2=p2+1;
        end
end
D=C;
i=0;
B1=E;   % \\ B1= List of nodes that have no relation
for kx1=1:numel(B1)
  C(E(kx1)-i,:)=[];
  C(:,E(kx1)-i)=[];
  i=i+1;
end
p3=0;
for i6=1:b1
   [row,colu]=find(B1==i6);
   co=numel(colu);
   if co==0
       p3=p3+1;
       F(1,p3)=i6;
   end
end
A=C;
%--------------------------------------------------------
b=size(A,1);
B=zeros(b,b);
sra=sum(A,2);
for i1=1:b
    for j1=1:b
        if (A(i1,j1)==1)
            B(i1,j1)=sra(j1);
        end
    end
end
%\\ srC=Add row 
%\\ cont= Counter The number of communications
srC=sum(B,2);
for i2=1:b
    as=0;
    cont=0;
    for j2=1:b
        if (A(i2,j2)~=0)
           as =srC(j2)+as;
           cont=cont+1;
        end
    end
    %\\ seg=Node communications 
    seg(i2,1)=srC(i2)/as;
end
%\\ so=Sort seg in descending order
 so=sort(seg,'descend');
 i3=1;
 while i3<=b
     [row,column]=find(seg==so(i3,1));
        %\\ nrw= Number of nodes equal
        nrw=numel(row);
        if nrw==1
            Ta(i3,1)=row(1,1);
            i3=i3+1;
        else
            for j3=1:nrw    
                %\\ Ta=Nodes ready for clustering
                Ta(i3,1)=row(j3,1);
                i3=i3+1;
            end
        end
 end
for j=1:b
   t=A(j,:);
   [rw,coln]=find(t==1);
    t=B(j,coln);
    E=min(t);
    minT(j,1)=E;      
end
ert=0;
p=0;
Clus{1,1}=[];
BT=1;
TaB=[];
for i=1:b
    BFT=0;
    ermd=[];
    termd3=[];
    termd2=[];
    termd1=[];
    termd=[];
     celltonum = cell2mat(Clus);
     s=Ta(i,:);
     %\\ Test the existence of a number in the cluster
     r=ismember(s,celltonum);
     if r==0
        x=A(s,:);
        [row,columnp]=find(x==1);
        y=B(s,columnp);
        my=min(y);
        Br=B(s,:);
        [row2,column2]=find(Br==my);
        nrc=numel(column2); 
        %\\ column2=The number of y nodes that the x node wants to cluster with
        BFT=0;
          for k=1:nrc
              ert=0;
                col=column2(1,k);
                if minT(col,1)==B(col,s) 
                      rr=ismember(col,celltonum); 
                      if rr==0
                           BFT=1;
                           ermd(k,1)=col;  
                           ermd(k,2)=1;   
                           ermd(k,3)=sra(col,1)*srC(col,1); 
                           ermd(k,4)=0; 
                           ermd(k,5)=1;
                           ermd(k,6)=k; 
                      else
                        BFT=1;
                          ermd(k,1)=col;
                          number=numel(Clus);
                          for kx=1:number   
                              na=ismember(col,Clus{kx}); 
                                 if na==true
                                    TF=0;
                                    break
                                 end
                          end
                          rclus=Clus{kx};
                          numrclus=numel(rclus); 
                          for z1=1:numrclus
                              if A(s,rclus(1,z1))==1
                                   ert=ert+1;
                              end
                          end
                          ermd(k,2)=ert;
                          sum1=0;
                          sum2=0;
                          for z2=1:numrclus
                              sum1=srC(rclus(1,z2),1)+sum1;
                              sum2=sra(rclus(1,z2),1)+sum2;
                          end
                          ermd(k,3)=(sum1/numrclus)*(sum2/numrclus);
                          ermd(k,4)=1;  
                          ermd(k,5)=numrclus;
                          ermd(k,6)=k;
                          ermd(k,7)=kx; 
                      end   
                end      
          end
          if BFT==0
               %\\ TaB= Table nodes B
               TaB(BT,1)=s; 
               BT=1+BT;    
          end
          if BFT==1
                m1=max(ermd(:,2));
                [r1,c1]=find(ermd(:,2)==m1);
                for h=1:numel(r1)
                    termd(h,:)=ermd(r1(h),:);
                end
                m2=min(termd(:,3));
                [r2,c2]=find(termd(:,3)==m2);
                for i=1:numel(r2)
                    termd1(i,:)=termd(r2(i),:);
                end
                m3=min(termd1(:,5));
                [r3,c3]=find(termd1(:,5)==m3);
                for i=1:numel(r3)
                    termd2(i,:)=termd1(r3(i),:);
                end
                m4=min(termd2(:,6));
                [r4,c4]=find(termd2(:,6)==m4);
                for i=1:numel(r4)
                    termd3(i,:)=termd2(r4(i),:);
                end
          else
            termd3(1,4)=1;  
          end
     else
         termd3(1,4)=1; 
     end 
        if termd3(1,4)==0
            p=p+1;
            Clus{1,p}=union(termd3(1,1),s);
        else
            if BFT==1
                Clus{termd3(1,7)}=union(Clus{termd3(1,7)},s);
            end
        end  
end 
%-------------------------------- Separate tables B and C  
at=0;
if at==0
 g=0;
 g1=0;
 TaB3=[];
 T2=[];
 nTaB=numel(TaB);
 if nTaB~=0
    cur=0;
   for i=1:nTaB
     celltonum = cell2mat(Clus);
     s=TaB(i,:);    
        Br=B(s,:);
        [row2,column2]=find(Br~=0);
        nrc=numel(column2); 
          for k=1:nrc              
              col=column2(1,k);
                      rr=ismember(col,celltonum); 
                      if rr==1
                         cur=1;
                      end
          end
          if cur==1
              g=g+1;
              TaB3(g,1)=s;
          else
              g1=g1+1;
              T2(g1,1)=s;
          end
          cur=0;
   end
 end
end
b3=numel(TaB3);
so1=sort(TaB3,'descend');
i6=1;
while i6<=b3
    fi=find(TaB3==so1(i6,1));
    for j5=1:numel(fi)
         so1(i6,2)=fi(1,j5);
         i6=i6+1;
    end
end
%----------------------------------Matrix Calculation B
atff1=0;
if atff1==0
    if numel(TaB3)~=0
      TaB=so1(:,1);
    end
nTaB=numel(TaB);
if (nTaB~=0)
for i=1:nTaB
    BFT=0;
    ermd=[];
    termd3=[];
    termd2=[];
    termd1=[];
    termd=[];
     celltonum = cell2mat(Clus);
     s=TaB(i,:);
     r=ismember(s,celltonum);
     if r==0       
        Br=B(s,:);
        [row2,column2]=find(Br~=0);
        nrc=numel(column2); 
        BFT=0;
          for k=1:nrc    
              ert=0;
              col=column2(1,k);
                      rr=ismember(col,celltonum); 
                      if rr==0
                           BFT=1;
                           ermd(k,1)=col;  
                           ermd(k,2)=1;  
                           ermd(k,3)=sra(col,1)*srC(col,1);   
                           ermd(k,4)=0;  
                           ermd(k,5)=1; 
                           ermd(k,6)=k; 
                      else
                        BFT=1;   
                          ermd(k,1)=col;
                          number=numel(Clus);
                          for kx=1:number   
                              na=ismember(col,Clus{kx}); 
                                 if na==true
                                    TF=0;
                                    break
                                 end
                          end
                          rclus=Clus{kx};
                          numrclus=numel(rclus);
                          for z1=1:numrclus
                              if A(s,rclus(1,z1))==1
                                   ert=ert+1;
                              end
                          end
                          ermd(k,2)=ert;
                          sum1=0;
                          sum2=0;
                          for z2=1:numrclus
                              sum1=srC(rclus(1,z2),1)+sum1;
                              sum2=sra(rclus(1,z2),1)+sum2;
                          end
                          ermd(k,3)=(sum1/numrclus)*(sum2/numrclus);
                          ermd(k,4)=1; 
                          ermd(k,5)=numrclus; 
                          ermd(k,6)=k;
                          ermd(k,7)=kx; 
                      end       
          end
          if BFT==1
                m1=max(ermd(:,2));
                [r1,c1]=find(ermd(:,2)==m1);
                for h=1:numel(r1)
                    termd(h,:)=ermd(r1(h),:);
                end
                m2=min(termd(:,3));
                [r2,c2]=find(termd(:,3)==m2);
                for i=1:numel(r2)
                    termd1(i,:)=termd(r2(i),:);
                end
                m3=min(termd1(:,5));
                [r3,c3]=find(termd1(:,5)==m3);
                for i=1:numel(r3)
                    termd2(i,:)=termd1(r3(i),:);
                end
                m4=min(termd2(:,2));
                [r4,c4]=find(termd2(:,2)==m4);
                for i=1:numel(r4)
                    termd3(i,:)=termd2(r4(i),:);
                end
                if termd(:,2)>1
                    termd3=termd;
                end
          else
            termd3(1,4)=1;
          end    
     else
         termd3(1,4)=1; 
     end    
     if numel(termd3(:,1))==1
        if termd3(1,4)==0
            p=p+1;
            Clus{1,p}=union(termd3(1,1),s);
        else
            if BFT==1
                Clus{termd3(1,7)}=union(Clus{termd3(1,7)},s);
            end
        end
     else
         B3=termd3(:,1);
         if termd3(1,7)~=0  
           for i6=1:numel(B3)
             Clus{termd3(1,7)}=union(Clus{termd3(i6,7)},Clus{termd3(1,7)});
           end
           Clus{termd3(1,7)}=union(Clus{termd3(1,7)},s);
         else
             p=p+1;
             Clus{1,p}=[];
           for i8=1:numel(B3) 
             Clus{1,p}=union(termd3(i8,1),Clus{1,p}); 
           end
            Clus{1,p}=union(Clus{1,p},s); 
         end
         if termd3(1,7)~=termd3(2,7)
           for i7=2:numel(B3)
             Clus{termd3(i,7)}=[];
           end
         end
     end
end
end
celltonum1 = cell2mat(Clus);
end
%-----------------------------------End of Matrix Calculation B

%-----------------------------------Matrix Calculation C
atff=0;
if atff==0
    TaB=T2;
nTaB=numel(TaB);
if (nTaB~=0)
for i=1:nTaB
    BFT=0;
    ermd=[];
    termd3=[];
    termd2=[];
    termd1=[];
    termd=[];
     celltonum = cell2mat(Clus);
     s=TaB(i,:);
     r=ismember(s,celltonum);
     if r==0       
        Br=B(s,:);
        [row2,column2]=find(Br~=0);
        nrc=numel(column2); 
        BFT=0;
          for k=1:nrc
              ert=0;
              col=column2(1,k);
                      rr=ismember(col,celltonum); 
                      if rr==0
                           BFT=1;
                           ermd(k,1)=col;  
                           ermd(k,2)=1;   
                           ermd(k,3)=sra(col,1)*srC(col,1);  
                           ermd(k,4)=0; 
                           ermd(k,5)=1;
                           ermd(k,6)=k; 
                      else
                          BFT=1;
                          ermd(k,1)=col;
                          number=numel(Clus);
                          for kx=1:number   
                              na=ismember(col,Clus{kx}); 
                                 if na==true
                                    TF=0;
                                    break
                                 end
                          end
                          rclus=Clus{kx};
                          numrclus=numel(rclus);
                          for z1=1:numrclus
                              if A(s,rclus(1,z1))==1
                                   ert=ert+1;
                              end
                          end
                          ermd(k,2)=ert;
                          sum1=0;
                          sum2=0;
                          for z2=1:numrclus
                              sum1=srC(rclus(1,z2),1)+sum1;
                              sum2=sra(rclus(1,z2),1)+sum2;
                          end
                          ermd(k,3)=(sum1/numrclus)*(sum2/numrclus);
                          ermd(k,4)=1;
                          ermd(k,5)=numrclus; 
                          ermd(k,6)=k;
                          ermd(k,7)=kx; 
                      end          
          end
          if BFT==1
                m1=max(ermd(:,2));
                [r1,c1]=find(ermd(:,2)==m1);
                for h=1:numel(r1)
                    termd(h,:)=ermd(r1(h),:);
                end
                m2=min(termd(:,3));
                [r2,c2]=find(termd(:,3)==m2);
                for i=1:numel(r2)
                    termd1(i,:)=termd(r2(i),:);
                end
                m3=min(termd1(:,5));
                [r3,c3]=find(termd1(:,5)==m3);
                for i=1:numel(r3)
                    termd2(i,:)=termd1(r3(i),:);
                end
                m4=min(termd2(:,2));
                [r4,c4]=find(termd2(:,2)==m4);
                for i=1:numel(r4)
                    termd3(i,:)=termd2(r4(i),:);
                end
                if termd(:,2)>1
                    termd3=termd;
                end
          else
            termd3(1,4)=1; 
          end    
     else
         termd3(1,4)=1; 
     end    
     if numel(termd3(:,1))==1
        if termd3(1,4)==0
            p=p+1;
            Clus{1,p}=union(termd3(1,1),s);
        else
            if BFT==1
                Clus{termd3(1,7)}=union(Clus{termd3(1,7)},s);
            end
        end
     else
          B3=termd3(:,1);
         if termd3(1,7)~=0  
           for i6=1:numel(B3)
             Clus{termd3(1,7)}=union(Clus{termd3(i6,7)},Clus{termd3(1,7)});
           end
           Clus{termd3(1,7)}=union(Clus{termd3(1,7)},s);
         else
             p=p+1;
             Clus{1,p}=[];
           for i8=1:numel(B3) 
             Clus{1,p}=union(termd3(i8,1),Clus{1,p}); 
           end
            Clus{1,p}=union(Clus{1,p},s); 
         end
         if termd3(1,7)~=termd3(2,7)
           for i7=2:numel(B3)
             Clus{termd3(i,7)}=[];
           end
         end 
     end  
end
end
celltonum1 = cell2mat(Clus);
end
toc
%-----------------------------------End of the algorithm

%--------------------------------------------------TurboMQ
for i=1:b
    for j=1:b
        if (A(i,j)==1)
            B(i,j)=j;
        end
    end
end
T={};
for i=1:b
    T{i,1}=i;
end
for i=1:b
    for j=1:b
        if B(i,j)~=0
            T{i,1}=[T{i,1} B(i,j)];
        end
    end
end
Items=[];
for i=1:numel(T)
    Items=union(Items,T{i});
end
Items=reshape(Items,1,[]);
Count=zeros(size(Items));
for i=1:numel(T)
    Count=Count+ismember(Items,T{i});
end
[Count,SortOrder]=sort(Count,'descend');
Items=Items(SortOrder);
count44=[];
    for i=1:numel(Clus)
        count33=[];
        for j=1:numel(Clus{i})
            count22=[];
            for t=1:numel(Clus)
                count22(t)=0;
            end
            a=numel(T{Clus{i}(j)})-1;
            for t=1:numel(Items)
                if ismember(Clus{i}(j),T{t}) && t~=Clus{i}(j)
                    a=a+1;
                end
            end
            for t=1:numel(Clus)
                if t~=i
                    for k=1:numel(Clus{t})
                        if ismember(Clus{i}(j),T{Clus{t}(k)})
                            count22(t)=count22(t)+1;
                        end
                        if ismember(Clus{t}(k),T{Clus{i}(j)})
                            count22(t)=count22(t)+1;
                        end
                    end
                end
            end
            Max=0;
            for t=1:numel(Clus)
                if Max<=count22(t)
                    Max=count22(t);
                    SC=t;
                end
            end
            b=0;
            for t=1:numel(Clus{i})
                if ismember(Clus{i}(j),T{Clus{i}(t)}) && t~=j
                    b=b+1;
                end
                if ismember(Clus{i}(t),T{Clus{i}(j)}) && t~=j
                    b=b+1;
                end
            end
            if b~=0 || count22(SC)~=0
                c=(b/a);
                if b<count22(SC)
                    c=(count22(SC)/a);
                end
            else
                c=0;
            end
            if b-count22(SC)~=0
                dd=0;
                ddd=0;
                if b~=0
                    dd=b/a;
                end
                if count22(SC)~=0
                    ddd=count22(SC)/a;
                end
                count33(j)=((dd-ddd)/c);
            else
                count33(j)=0;
            end     
        end
        a=0;
        b=0;
        for j=1:numel(Clus{i})
            a=a+count33(j);
        end
        b=numel(Clus{i});
        if a~=0
            count44(i)=(a/b);
        else
            count44(i)=0;
        end
    end
    a=0;
    b=0;
    for j=1:numel(Clus)
        a=a+count44(j);
    end
    b=numel(Clus);
    if a~=0
        Zrib=(a/b);
    else
        Zrib=0;
    end
    nB1=numel(B1);
    oi=0;
    for e=1:b
        numberc=numel(Clus{e});
        if numberc==0
            oi=oi+1;
        end
    end
    Zrib=(Zrib*numel(Clus))/(numel(Clus)-oi);
    if numel(B1)~=0
      Zrib=(Zrib*numel(Clus))/(numel(Clus)+1);
    end
 %
count11=[];
for i=1:numel(Clus)
    count11(i)=0;
end
for i=1:numel(Clus)
    count2=[];
    for j=1:numel(Clus{i})
        count2(j)=0;
    end
    count3=[];
    for j=1:numel(Clus{i})
        count3(j)=(numel(T{Clus{i}(j)})-1);
    end
    for j=1:numel(Clus{i})
        for t=1:numel(Clus{i})
            if t~=j
                if ismember(Clus{i}(j),T{Clus{i}(t)})
                    count2(t)=count2(t)+1;
                end
            end
        end
    end
    for j=1:numel(Clus{i})
        count3(j)=count3(j)-count2(j);
    end
    a=0;
    b=0;
    for j=1:numel(Items)
        if ~ismember(Items(j),Clus{i})
            for t=1:numel(Clus{i})
                if ismember(Clus{i}(t),T{Items(j)})
                    count3(t)=count3(t)+1;
                end
            end
        end
    end
    for j=1:numel(Clus{i})
        a=a+count2(j);
        b=b+count3(j);
    end
    if a~=0
        count11(i)=((2*a)/((2*a)+(b)));
    end   
end
TurboMQ=0;
for i=1:numel(Clus)
    TurboMQ=TurboMQ+count11(i);
end
%+==============================
if numel(A)~= numel(A1)
    nu=numel(Clus);
    for i2=1:nu
       Clus1{i2}=[];
       for j=1:numel(F)
          if ismember(j,Clus{i2})
             Clus1{i2}=union(F(1,j),Clus1{i2});
          end 
       end
    end
    Clus=Clus1; 
end
%+==============================
disp('TurboMQ=');
disp(TurboMQ);
disp('Number of Cluster=');
nB1=numel(B1);
if nB1==1
    Clus{1,1}=union(B1(1,1),Clus{1,1});
else
    if nB1~=0
      p=p+1;
      Clus{1,p}=[];
      for i=1:nB1
        Clus{1,p}=union(B1(1,i),Clus{1,p});
      end
    end
end
disp(numel(Clus)-oi);
nu2=numel(Clus);
tmojo=[];
for j6=1:nu2
    ku=Clus{j6};
    nku=numel(ku);
    for k6=1:nku
        tmojo(ku(1,k6),1)=j6;
    end
end